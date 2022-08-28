#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEADER_PREFIX '>'

#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

void print_entry(char* chr, int start, int end, char* unit) {
  printf("%s\t%i\t%i\tunit=%s\n", chr, start, end, unit);
}

typedef struct {
  int n;
  int* elements;
} Ring;

typedef struct {
  int d;
  int n_empty;
  int n_required;
  Ring* matches;
} Divisor;

typedef struct {
  Ring* last_bases;
  Divisor** divisors;
  int div_index;
  int n_divisors;
  int length;
  char* chr;
  char* unit_buffer;
  int is_even;
} SeqState;

Ring* init_ring (const int n) {
  Ring* r;

  r = malloc(sizeof(*r));
  r->n = n;
  r->elements = malloc(n * sizeof(int));

  return r;
}

void free_ring (Ring* rng) {
  free(rng->elements);
  free(rng);
}

void write_ring (const Ring* r, const int i, const int c) {
  r->elements[i % r->n] = c;
}

int read_ring (const Ring* const r, const int i) {
  return r->elements[i % r->n];
}

/* divisors */

int int_compare (const void* a, const void* b) {
  const int* a0 = (const int*)a;
  const int* b0 = (const int*)b;

  return (a0 - b0);
}

Divisor* init_divisor (int d, int r) {
  Divisor* div;

  int n; 
  int empty;

  n = d * ((r / d) - 1);
  empty = d % n;

  div = malloc(sizeof(*div));

  div->d = d;
  div->n_required = n;
  div->n_empty = empty;
  div->matches = init_ring(empty == 0 ? r / 2 : r);

  return div;
}

void free_divisor (Divisor* div) {
  free_ring(div->matches);
  free(div);
}

/* int* sieve_of_eratosthenes (const int n) { */
/*   /\* ASSUME n > 1 *\/ */
/*   int i, j, r; */
/*   int* candidates; */

/*   candidates = malloc(n * sizeof(*candidates)); */

/*   for (i = 0; i < n; i++) { */
/*     candidates[i] = 1; */
/*   } */

/*   r = sqrt(n); */

/*   for (i = 2; i <= r; i++) { */
/*     if (candidates[i]) { */
/*       for(j = i * i; j <= n; j += i) { */
/*         candidates[j] = 0; */
/*       } */
/*     } */
/*   } */

/*   return candidates; */
/* } */

int* expand_array (int n, int* result) {
  int* tmp;

  tmp = realloc(result, n * sizeof(*result));

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

void init_subpatterns (SeqState* st, int r) {
  Divisor** subpatterns;

  int* xdivisors = NULL;
  int i, k;

  xdivisors = find_exclusive_divisors(r, &k);

  subpatterns = malloc(k * sizeof(*subpatterns));

  fprintf(stderr, "%i\n", k);

  for (i = 0; i < k; i++) {
    subpatterns[i] = init_divisor(xdivisors[i], r);
  }

  st->divisors = subpatterns;
  st->n_divisors = k;
}

void print_entryN (SeqState* st, int i, int n) {
  int r;
  int o;
  int j;

  if (n >= st->length) {
    r = st->last_bases->n;
    o = (i - n) % r;

    for (j = 0; j < r; j++) {
      st->unit_buffer[j] = read_ring(st->last_bases, (j + o));
    }

    print_entry(st->chr, i - n, i, st->unit_buffer);
  }
}

int all_true (Divisor* div, const int i) {
  int j;

  for (j = 0; j < div->n_required; j++) {
    if (!read_ring(div->matches, i - j)) {
      return 0;
    }
  }

  return 1;
}

int invalid_repeat (SeqState* st, const int i) {
  Divisor** d;

  for (d = st->divisors; d < st->divisors + st->n_divisors; d++) {
    if (all_true(*d, i)) {
      return 1;
    }
  }
  return 0;

  /* int j; */
  /* for (j = st->n_divisors - 1; j >= 0; j--) { */
  /*   if (!all_true(st->divisors[j], i)) { */
  /*     return 0; */
  /*   } */
  /* } */
  /* return 1; */
}

void update_match (Ring* last_bases, Divisor* div, const int i) {
  int c0;
  int c1;

  c0 = read_ring(last_bases, i);
  c1 = read_ring(last_bases, (i - div->d));

  write_ring(div->matches, i, c0 == c1);
}

void update_matches (SeqState* st, const int i) {
  for (int j = 0; j < st->div_index; j++) {
    update_match(st->last_bases, st->divisors[j], i);
  }
}

SeqState* init_seq_state (char* chr, int r, int l) {
  SeqState* st;

  st = malloc(sizeof(*st));

  st->last_bases = init_ring(r);
  st->div_index = 0;
  st->length = l;
  st->chr = chr;
  st->unit_buffer = malloc((r + 1) * sizeof(char));
  st->unit_buffer[r] = '\0';
  st->is_even = !(r % 2);

  init_subpatterns(st, r);

  return st;
}

void free_seq_state (SeqState* st) {
  for (int j = 0; j < st->n_divisors; j++) {
    free_divisor(st->divisors[j]);
  }

  free_ring(st->last_bases);

  free(st->divisors);
  free(st->unit_buffer);

  free(st);
}

void update_blind (SeqState* st, int i, int n, int r) {
  int shift;
  int j;
  int k;

  shift = (n + 1) % r;

  if (shift > 1) {
    for (j = 0; j < st->n_divisors - st->is_even; j++) {
      for (k = max(1, shift - st->divisors[j]->n_empty); k < shift; k++) {
        update_match(st->last_bases, st->divisors[j], i - k);
      }
    }
  }
}

void update_divisor_index (SeqState* st, int n) {
  if (st->div_index < st->n_divisors && st->divisors[st->div_index]->d == n + 1) {
    st->div_index++;
  }
}

void scan_seqN(FILE* fp, char* chr, const int r, const int l) {
  SeqState* st;

  int c;
  int c0;
  int i = 0;
  int n = 0;

  st = init_seq_state(chr, r, l);

  while (1) {
    c = fgetc(fp);

    // ignore newlines
    if (c != '\n') {

      if (c == EOF || c == HEADER_PREFIX) {
        print_entryN(st, i, n);
        break;

      } else if (c == 'N') {
        print_entryN(st, i, n);
        n = 0;
        st->div_index = 0;

      } else if (n < r - 1) {
        write_ring(st->last_bases, i, c);
        update_matches(st, i);
        update_divisor_index(st, n);
        n++;

      } else if (n == r - 1) {
        write_ring(st->last_bases, i, c);
        update_matches(st, i);
        n = r - invalid_repeat(st, i);

      } else {
        c0 = read_ring(st->last_bases, i);
        if (c0 == c) {
          n++;
        } else {
          print_entryN(st, i, n);
          write_ring(st->last_bases, i, c);
          update_matches(st, i);
          // update all divisors with blind arrays
          update_blind(st, i, n, r);
          n = r - invalid_repeat(st, i);
        }
      }
      i++;
    }
  }

  free_seq_state(st);
}

void scan_seq1(FILE* fp, char* chr, int l) {
  char lastBase[2] = {'N', '\0'};
  int i = 0;
  int n = 1;
  int c;

  while (1) {
    c = fgetc(fp);

    // ignore newlines
    if (c != '\n') {

      if (c == EOF || c == HEADER_PREFIX) {
        if (n >= l) {
          print_entry(chr, i - n, i, lastBase);
        }
        break;

      } else if (c == 'N') {
        n = 1;
        lastBase[0] = c;

      } else if (c == lastBase[0]) {
        n++;

      } else {
        if (n >= l) {
          print_entry(chr, i - n, i, lastBase);
        }
        n = 1;
        lastBase[0] = c;
      }
      i++;
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

int read_fasta(FILE* fp, int r, int l) {
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

    if (r == 1) {
      scan_seq1(fp, chr, l);
    } else {
      scan_seqN(fp, chr, r, l);
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  if ( argc == 4 ) {
    int r;
    int l;

    r = atoi(argv[1]);
    l = atoi(argv[2]);

    read_fasta(fopen(argv[3], "r"), r, l);
  } else {
    printf("Usage: REPS LENGTH INFILE\n");
  }
}
