#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEADER_PREFIX '>'

void print_entry(char* chr, int start, int end, char* unit) {
  printf("%s\t%i\t%i\tunit=%s\n", chr, start, end, unit);
}

struct ring {
  int n;
  int* elements;
};

struct ring* init_ring (int n) {
  struct ring* r = malloc(sizeof(*r));
  r->n = n;
  r->elements = malloc(n * sizeof(int));
  return r;
}

void write_ring (struct ring* r, int i, int c) {
  r->elements[i % r->n] = c;
}

int read_ring (struct ring* r, int i) {
  return r->elements[i % r->n];
}

struct divisor {
  int d;
  struct ring* matches;
};

struct seq_state {
  struct ring* last_bases;
  // TODO make this an array of pointers and not literal divisors?
  struct divisor* divisors;
  int div_index;
  int n_divisors;
};

int compare (const void* a, const void* b) {
  const struct divisor* a0 = (const struct divisor*)a;
  const struct divisor* b0 = (const struct divisor*)b;
  return (a0->d - b0->d);
}

struct divisor init_divisor(int d, int r) {
  struct divisor div;

  div.d = d;
  div.matches = init_ring(d * ((r / d) - 1));
  return div;
}

void find_divisors(struct seq_state* st, int r) {
  double u = floor(sqrt(r));
  int i = 2;
  int n = 1;
  int one_div;
  struct divisor* divisors = malloc(n * sizeof(*divisors));

  divisors[0] = init_divisor(1, r);

  while (i <= u) {
    if (r % i == 0) {
      one_div = r / i == i;
      n = n + 1 + !one_div;
      // lets just pretend realloc won't return NULL here
      divisors = realloc(divisors, n * sizeof(*divisors));
      divisors[n - 1] = init_divisor(i, r);
      if (!one_div) {
        divisors[n - 2] = init_divisor(r / i, r);
      }
    }
    i++;
  }

  qsort(divisors, n, sizeof(*divisors), compare);

  /* for (i = 0; i < n; i++) { */
  /*   fprintf( stderr, "%i\n", divisors[i].d); */
  /* } */

  st->divisors = divisors;
  st->n_divisors = n;
}

void print_entryN (struct seq_state* st, char* chr, int i, int n, int l) {
  if (n >= l) {
    int r = st->last_bases->n;
    int o = (i - n) % r;
    char* s = malloc((r + 1) * sizeof(char));
    int j;

    for (j = 0; j < r; j++) {
      s[j] = read_ring(st->last_bases, (j + o));
    }

    s[r] = '\0';

    print_entry(chr, i - n, i, s);

    free(s);
  }
}

int all_true (struct ring* ring) {
  int i;

  /* fprintf(stderr, "all-true: %i\n", ring->n); */
  for (i = 0; i < ring->n; i++) {
    /* fprintf(stderr, "%i", ring->elements[i]); */
    if (ring->elements[i] == 0) {
      return 0;
    }
  }

  return 1;
}

int valid_repeat (struct seq_state* st) {
  int i;

  for (i = 0; i < st->n_divisors; i++) {
    if (all_true(st->divisors[i].matches) == 1) {
      /* fprintf(stderr, "false\n"); */
      return 0;
    }
  }
  /* fprintf(stderr, "true\n"); */

  return 1;
}

void update_match (struct ring* last_bases, struct divisor* div, int i) {
  int c0 = read_ring(last_bases, i);
  int c1 = read_ring(last_bases, (i - div->d));
  /* fprintf(stderr, "%i %c %c %i %i\n", i % div->matches->n, c0, c1, div->d, c0 == c1); */

  write_ring(div->matches, i, c0 == c1);
}

void update_matches (struct seq_state* st, int i) {
  int j;

  for (j = 0; j < st->div_index; j++) {
    update_match(st->last_bases, &st->divisors[j], i);
  }
}

int find_next_n (struct seq_state* st, int i, int r, char c) {
  write_ring(st->last_bases, i, c);
  update_matches(st, i);
  /* fprintf(stderr, "valid: %i\n", valid_repeat(st)); */
  return r - !valid_repeat(st);
}

void scan_seqN(FILE* fp, char* chr, int r, int l) {
  struct seq_state* st = malloc(sizeof(*st));
  int c;
  int c0;
  int i = 0;
  int n = 0;

  st->last_bases = init_ring(r);
  st->div_index = 0;
  find_divisors(st, r);

  /* for (i = 0; i < st->n_divisors; i++) { */
  /*   fprintf(stderr, "%i\n", st->divisors[i].matches->n); */
  /* } */

  /* i = 0; */
  /* int j; */
  /* int k; */

  while (1) {
    c = fgetc(fp);

    // ignore newlines
    if (c != '\n') {

      /* fprintf(stderr, "base = %c\n", c); */
      /* /\* fprintf(stderr, "divisor_index = %i\n", st->div_index); *\/ */

      /* for (j = 0; j < r; j++) { */
      /*   fprintf(stderr, "%c", st->last_bases->elements[j]); */
      /* } */
      /* fprintf(stderr, "\n"); */

      /* for (j = 0; j < st->n_divisors; j++) { */
      /*   for (k = 0; k < st->divisors[j].matches->n; k++) { */
      /*     fprintf(stderr, "%i", st->divisors[j].matches->elements[k]); */
      /*   } */
      /*   fprintf(stderr, "\n"); */
      /* } */

      /* fprintf(stderr, "n = %i\n", n); */

      if (c == EOF || c == HEADER_PREFIX) {
        print_entryN(st, chr, i, n, l);
        break;

      } else if (c == 'N') {
        print_entryN(st, chr, i, n, l);
        n = 0;
        st->div_index = 0;

      } else if (n < r - 1) {
        write_ring(st->last_bases, i, c);
        update_matches(st, i);
        /* fprintf(stderr, "%i\n", st->divisors[st->div_index].d); */
        if (st->div_index <= st->n_divisors && st->divisors[st->div_index].d == n + 1) {
          st->div_index++;
        }
        n++;

      } else if (n == r - 1) {
        n = find_next_n(st, i, r, c);

      } else {
        c0 = read_ring(st->last_bases, i);
        if (c0 == c) {
          update_matches(st, i);
          n++;
        } else {
          print_entryN(st, chr, i, n, l);
          n = find_next_n(st, i, r, c);
          /* fprintf(stderr, "%i\n", n); */
        }
      }
      i++;
      /* fprintf(stderr, "---\n"); */
    }
  }

  /* TODO also free divisors? */
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
