#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEADER_PREFIX '>'

/* helper functions */

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


int compare (const void* a, const void* b) {
  const Divisor* a0 = (const Divisor*)a;
  const Divisor* b0 = (const Divisor*)b;

  return (a0->d - b0->d);
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

void find_divisors(SeqState* st, int r) {
  Divisor** divisors;

  int i = 2;
  int n = 1;
  int is_mid;

  double u;

  divisors = malloc(n * sizeof(*divisors));

  divisors[0] = init_divisor(1, r);

  u = floor(sqrt(r));

  while (i <= u) {
    if (r % i == 0) {
      is_mid = r / i == i;
      n = n + 1 + !is_mid;
      // lets just pretend realloc won't return NULL here
      divisors = realloc(divisors, n * sizeof(*divisors));
      divisors[n - 1] = init_divisor(i, r);
      if (!is_mid) {
        divisors[n - 2] = init_divisor(r / i, r);
      }
    }
    i++;
  }

  qsort(divisors, n, sizeof(*divisors), compare);

  st->divisors = divisors;
  st->n_divisors = n;
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

int valid_repeat (SeqState* st, const int i) {
  Divisor** d;

  for (d = st->divisors; d < st->divisors + st->n_divisors; d++) {
    if (all_true(*d, i)) {
      return 0;
    }
  }
  return 1;
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

  find_divisors(st, r);

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

void scan_seqN(FILE* fp, char* chr, const int r, const int l) {
  SeqState* st;

  int c;
  int c0;
  int i = 0;
  int n = 0;
  int j;
  int k;
  int shift;

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
        if (st->div_index < st->n_divisors && st->divisors[st->div_index]->d == n + 1) {
          st->div_index++;
        }
        n++;

      } else if (n == r - 1) {
        write_ring(st->last_bases, i, c);
        update_matches(st, i);
        n = r - !valid_repeat(st, i);

      } else {
        c0 = read_ring(st->last_bases, i);
        if (c0 == c) {
          n++;
        } else {
          print_entryN(st, i, n);
          write_ring(st->last_bases, i, c);
          update_matches(st, i);

          // update all divisors with blind arrays
          shift = (n + 1) % r;
          if (shift > 1) {
            for (j = 0; j < st->n_divisors - st->is_even; j++) {
              for (k = max(1, shift - st->divisors[j]->n_empty); k < shift; k++) {
                update_match(st->last_bases, st->divisors[j], i - k);
              }
            }
          }

          n =  r - !valid_repeat(st, i);
        }
      }
      i++;
    }
  }

  free_seq_state(st);
}

void scan_seq1(FILE* fp, char* chr, int len) {
  char lastBase[2] = {'N', '\0'};
  long i = 0;
  long n = 1;
  int c;

  do {
    c = fgetc(fp);
    if (c == '\n') {
      // ignore newlines
    } else if (c == 'N') {
      i++;
      n = 1;
      lastBase[0] = c;
    } else if (c == lastBase[0]) {
      i++;
      n++;
    } else {
      if (n >= len) {
        print_entry(chr, i - n, i, lastBase);
      }
      i++;
      n = 1;
      lastBase[0] = c;
    }
  } while (c != EOF && c != HEADER_PREFIX);
}

int seek_char(FILE* fp, char t) {
  int c;

  do {
    c = fgetc(fp);
  } while (EOF != c && t != c);

  return c;
}

int parse_header(FILE* fp, char* chr) {
  int c = fscanf(fp, "%31s", chr);
  if (1 != c && EOF != c) {
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

  while (EOF != parse_header(fp, chr)) {
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
    int r = atoi(argv[1]);
    int l = atoi(argv[2]);
    read_fasta(fopen(argv[3], "r"), r, l);
  } else {
    printf("Usage: REPS LENGTH INFILE\n");
  }
}
