#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEADER_PREFIX '>'

/* #define min(a,b)             \ */
/* ({                           \ */
/*     __typeof__ (a) _a = (a); \ */
/*     __typeof__ (b) _b = (b); \ */
/*     _a < _b ? _a : _b;       \ */
/* }) */

#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

void print_entry(char* chr, int start, int end, char* unit) {
  printf("%s\t%i\t%i\tunit=%s\n", chr, start, end, unit);
}

struct ring {
  int n;
  int* elements;
};

struct ring* init_ring (const int n) {
  struct ring* r;

  r = malloc(sizeof(*r));
  r->n = n;
  r->elements = malloc(n * sizeof(int));

  return r;
}

void write_ring (const struct ring* r, const int i, const int c) {
  r->elements[i % r->n] = c;
}

int read_ring (const struct ring* const r, const int i) {
  return r->elements[i % r->n];
}

struct divisor {
  int d;
  int n_empty;
  int n_required;
  struct ring* matches;
};

struct seq_state {
  struct ring* last_bases;
  // TODO make this an array of pointers and not literal divisors?
  struct divisor* divisors;
  int div_index;
  int n_divisors;
  int length;
  char* chr;
  char* unit_buffer;
};

int compare (const void* a, const void* b) {
  const struct divisor* a0 = (const struct divisor*)a;
  const struct divisor* b0 = (const struct divisor*)b;

  return (a0->d - b0->d);
}

struct divisor init_divisor(int d, int r) {
  struct divisor div;

  int n; 
  int empty;

  n = d * ((r / d) - 1);
  empty = d % n;

  div.d = d;
  div.n_required = n;
  div.n_empty = empty;
  div.matches = init_ring(empty == 0 ? r / 2 : r);

  return div;
}

void find_divisors(struct seq_state* st, int r) {
  struct divisor* divisors;

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

  /* for (int j = 0; j < n; j++) { */
  /*   fprintf(stderr, "%i\n", divisors[j].n_required); */
  /* } */

  st->divisors = divisors;
  st->n_divisors = n;
}

void print_entryN (struct seq_state* st, int i, int n) {
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

int all_true (struct divisor* div, const int i) {
  int j;

  for (j = 0; j < div->n_required; j++) {
    if (!read_ring(div->matches, i - j)) {
      return 0;
    }
  }

  return 1;
}

int valid_repeat (struct seq_state* st, const int i) {
  struct divisor* d;

  for (d = st->divisors; d < st->divisors + st->n_divisors; d++) {
    if (all_true(d, i)) {
      return 0;
    }
  }
  return 1;
}

void update_match (struct ring* last_bases, struct divisor* div, const int i) {
  int c0;
  int c1;

  c0 = read_ring(last_bases, i);
  c1 = read_ring(last_bases, (i - div->d));

  write_ring(div->matches, i, c0 == c1);
}

void update_matches (struct seq_state* st, const int i) {
  for (int j = 0; j < st->div_index; j++) {
    update_match(st->last_bases, &st->divisors[j], i);
  }
}

/* int find_next_n (struct seq_state* st, int i, int r, char c) { */
/*   write_ring(st->last_bases, i, c); */
/*   update_matches(st, i); */
/*   return r - !valid_repeat(st, i); */
/* } */

void scan_seqN(FILE* fp, char* chr, const int r, const int l) {
  struct seq_state* st;

  int c;
  int c0;
  int i = 0;
  int n = 0;
  int j;
  int k;
  int shift;

  st = malloc(sizeof(*st));

  st->last_bases = init_ring(r);
  st->div_index = 0;
  st->length = l;
  st->chr = chr;
  st->unit_buffer = malloc((r + 1) * sizeof(char));
  st->unit_buffer[r] = '\0';

  find_divisors(st, r);

  while (1) {
    c = fgetc(fp);

    /* fprintf(stderr, "---\n"); */
    /* fprintf(stderr, "i = %i\n", i); */
    /* fprintf(stderr, "n = %i\n", n); */
    /* fprintf(stderr, "c = %c\n", c); */
    /* fprintf(stderr, "valid = %i\n", valid_repeat(st, i)); */
    /* for (int j = 0; j < r; j++) { */
    /*   fprintf(stderr, "%c", read_ring(st->last_bases, i - j)); */
    /* } */
    /* fprintf(stderr, "\n"); */
    /* for (int j = 0; j < st->n_divisors; j++) { */
    /*   fprintf(stderr, "d = %i\n", st->divisors[j].d); */
    /*   /\* for (int k = 0; k < st->divisors[j].n_required; k++) { *\/ */
    /*   for (int k = 0; k < st->divisors[j].matches->n; k++) { */
    /*     fprintf(stderr, "%i", read_ring(st->divisors[j].matches, i - k)); */
    /*   } */
    /*   fprintf(stderr, "\n"); */
    /* } */

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
        if (st->div_index <= st->n_divisors && st->divisors[st->div_index].d == n + 1) {
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
            for (j = 0; j < st->n_divisors; j++) {
              if (st->divisors[j].n_empty > 0) {
                for (k = max(1, shift - st->divisors[j].n_empty); k < shift; k++) {
                  update_match(st->last_bases, &st->divisors[j], i - k);
                }
              }
            }
          }

          n =  r - !valid_repeat(st, i);
        }
      }
      i++;
    }
  }

  /* TODO also free divisors? */
  free(st->unit_buffer);
  free(st);
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
