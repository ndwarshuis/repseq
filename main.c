#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEADER_PREFIX '>'

void print_entry(char* chr, int p, int n, char* unit) {
  printf("%s\t%i\t%i\tunit=%s\n", chr, p - n, p, unit);
}

typedef struct {
  int size;
  int* elements;
} Ring;

typedef struct {
  int d;
  Ring* matches;
} SubMatch;

typedef struct {
  Ring* last_bases;
  SubMatch** submatches;
  int submatch_index;
  int n_submatches;
  int length;
  char* chr;
  char* unit_buffer;
} SeqState;

Ring* init_ring (const int size) {
  Ring* ring;

  ring = malloc(sizeof(*ring));
  ring->size = size;
  ring->elements = malloc(size * sizeof(int));

  return ring;
}

void free_ring (Ring* ring) {
  free(ring->elements);
  free(ring);
}

void write_ring (const Ring* ring, const int i, const int c) {
  ring->elements[i % ring->size] = c;
}

int read_ring (const Ring* const ring, const int i) {
  return ring->elements[i % ring->size];
}

int int_compare (const void* a, const void* b) {
  const int* a0 = (const int*)a;
  const int* b0 = (const int*)b;

  return (a0 - b0);
}

SubMatch* init_submatch (int d, int rep) {
  SubMatch* sm;

  sm = malloc(sizeof(*sm));

  sm->d = d;
  sm->matches = init_ring(d * ((rep / d) - 1));

  return sm;
}

void free_submatch (SubMatch* sm) {
  free_ring(sm->matches);
  free(sm);
}

int* expand_array (int size, int* result) {
  int* tmp;

  tmp = realloc(result, size * sizeof(*result));

  if (tmp == NULL) {
    printf("Error expanding array");
    exit(-1);
  } else {
    return tmp;
  }
}

int* find_divisors (int x, int* k) {
  int* candidates;

  int i = 2;
  int is_mid;
  int u;

  *k = 1;

  candidates = malloc(*k * sizeof(int));

  candidates[0] = 1;

  u = floor(sqrt(x));

  while (i <= u) {
    if (x % i == 0) {
      is_mid = x / i == i;
      *k = *k + 1 + !is_mid;
      candidates = expand_array(*k, candidates);
      candidates[*k - 1] = i;
      if (!is_mid) {
        candidates[*k - 2] = x / i;
      }
    }
    i++;
  }

  qsort(candidates, *k, sizeof(int), int_compare);

  return candidates;
}

int* find_exclusive_divisors (int x, int* k) {
  int i, j, d;
  int* candidates;
  int* result;

  *k = 1;

  result = malloc(*k * sizeof(int));

  candidates = find_divisors(x, &d);

  for (i = 0; i < d - 1; i++) {
    j = i + 1;
    while (candidates[j] % candidates[i] > 0 && j < d) {
      j++;
    }
    if (j == d) {
      (*k)++;
      result = expand_array(*k, result);
      result[*k - 2] = candidates[j];
    }
  }

  result[*k - 1] = candidates[d - 1];

  free(candidates);

  return result;
}

void init_submatches (SeqState* st, int rep) {
  SubMatch** sms;

  int* xdivisors = NULL;
  int i, k;

  xdivisors = find_exclusive_divisors(rep, &k);

  sms = malloc(k * sizeof(*sms));

  for (i = 0; i < k; i++) {
    sms[i] = init_submatch(xdivisors[i], rep);
  }

  st->submatches = sms;
  st->n_submatches = k;
}

void print_entryN (SeqState* st, int p, int n) {
  int r, o, j;

  if (n >= st->length) {
    r = st->last_bases->size;
    o = (p - n) % r;

    for (j = 0; j < r; j++) {
      st->unit_buffer[j] = read_ring(st->last_bases, (j + o));
    }

    print_entry(st->chr, p, n, st->unit_buffer);
  }
}

int all_match (SubMatch* sm) {
  for (int i = 0; i < sm->matches->size; i++) {
    if (!sm->matches->elements[i]) {
      return 0;
    }
  }

  return 1;
}

int invalid_repeat (SeqState* st) {
  SubMatch** d;

  for (d = st->submatches; d < st->submatches + st->n_submatches; d++) {
    if (all_match(*d)) {
      return 1;
    }
  }
  return 0;
}

void update_match (Ring* last_bases, SubMatch* sm, const int p) {
  int c0;
  int c1;

  c0 = read_ring(last_bases, p);
  c1 = read_ring(last_bases, (p - sm->d));

  write_ring(sm->matches, p, c0 == c1);
}

void update_matches (SeqState* st, const int p) {
  for (int j = 0; j < st->submatch_index; j++) {
    update_match(st->last_bases, st->submatches[j], p);
  }
}

SeqState* init_seq_state (char* chr, int rep, int len) {
  SeqState* st;

  st = malloc(sizeof(*st));

  st->last_bases = init_ring(rep);
  st->submatch_index = 0;
  st->length = len;
  st->chr = chr;
  st->unit_buffer = malloc((rep + 1) * sizeof(char));
  st->unit_buffer[rep] = '\0';

  init_submatches(st, rep);

  return st;
}

void free_seq_state (SeqState* st) {
  for (int j = 0; j < st->n_submatches; j++) {
    free_submatch(st->submatches[j]);
  }

  free_ring(st->last_bases);

  free(st->submatches);
  free(st->unit_buffer);

  free(st);
}

void update_submatch_index (SeqState* st, int n) {
  if (st->submatch_index < st->n_submatches && st->submatches[st->submatch_index]->d == n + 1) {
    st->submatch_index++;
  }
}

int next_n (SeqState* st, int p, int c) {
  write_ring(st->last_bases, p, c);
  update_matches(st, p);
  return st->last_bases->size - invalid_repeat(st);
}

void scan_seqN(FILE* fp, char* chr, const int rep, const int len) {
  SeqState* st;

  int c, c0;
  int p = 0;
  int n = 0;

  st = init_seq_state(chr, rep, len);

  while (1) {
    c = fgetc(fp);

    // ignore newlines
    if (c != '\n') {

      if (c == EOF || c == HEADER_PREFIX) {
        print_entryN(st, p, n);
        break;

      } else if (c == 'N') {
        print_entryN(st, p, n);
        n = 0;
        st->submatch_index = 0;

      } else if (n < rep - 1) {
        write_ring(st->last_bases, p, c);
        update_matches(st, p);
        update_submatch_index(st, n);
        n++;

      } else if (n == rep - 1) {
        n = next_n(st, p, c);

      } else {
        c0 = read_ring(st->last_bases, p);
        if (c0 == c) {
          n++;
          update_matches(st, p);
        } else {
          print_entryN(st, p, n);
          n = next_n(st, p, c);
        }
      }
      p++;
    }
  }

  free_seq_state(st);
}

void scan_seq1(FILE* fp, char* chr, int len) {
  char lastBase[2] = {'N', '\0'};
  int p = 0;
  int n = 1;
  int c;

  while (1) {
    c = fgetc(fp);

    // ignore newlines
    if (c != '\n') {

      if (c == EOF || c == HEADER_PREFIX) {
        if (n >= len) {
          print_entry(chr, p, n, lastBase);
        }
        break;

      } else if (c == 'N') {
        n = 1;
        lastBase[0] = c;

      } else if (c == lastBase[0]) {
        n++;

      } else {
        if (n >= len) {
          print_entry(chr, p, n, lastBase);
        }
        n = 1;
        lastBase[0] = c;
      }
      p++;
    }
  }
}

int seek_char(FILE* fp, char t) {
  int c;

  do {
    c = fgetc(fp);
  } while (EOF != c && t != c);

  return c;
}

int parse_header(FILE* fp, char* chr) {
  int c;

  c = fscanf(fp, "%31s", chr);

  if (c != 1 && c != EOF) {
    printf("Error when parsing chromosome header");
    exit(-1);
  }

  return seek_char(fp, '\n');
}

int read_fasta(FILE* fp, int rep, int len) {
  char chr[32];
  
  if (fp == NULL) {
    perror("Error in opening file");
    return (-1);
  }

  // find the first header (which should be the first character, but genomics
  // is weird...
  seek_char(fp, HEADER_PREFIX);


  while (parse_header(fp, chr) != EOF) {
    fprintf(stderr, "Parsing chromosome %s\n", chr);

    if (rep == 1) {
      scan_seq1(fp, chr, len);
    } else {
      scan_seqN(fp, chr, rep, len);
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  if (argc == 4) {
    int rep;
    int len;

    rep = atoi(argv[1]);
    len = atoi(argv[2]);

    read_fasta(fopen(argv[3], "r"), rep, len);
  } else {
    printf("Usage: REPS LENGTH INFILE\n");
  }
}
