#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEADER_PREFIX '>'

/******************************************************************************* 
  Various helper functions
*/

void print_entry(char* chr, int p, int n, char* unit) {
  printf("%s\t%i\t%i\tunit=%s\n", chr, p - n, p, unit);
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

/******************************************************************************* 
  Divisor functions
*/

int int_compare (const void* a, const void* b) {
  const int* a0 = (const int*)a;
  const int* b0 = (const int*)b;

  return (a0 - b0);
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

/*
  'exclusive divisor' = 'divisor that doesn't divide any other divisors'
*/
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

/******************************************************************************* 
  Simple ring buffer
*/

typedef struct {
  int size;
  int* elements;
} Ring;

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

/******************************************************************************* 
  Submatch buffer (for testing if a given repeat has smaller repeats in it)
*/

typedef struct {
  int d;
  Ring* matches;
} SubMatch;

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

int all_match (SubMatch* sm) {
  for (int i = 0; i < sm->matches->size; i++) {
    if (!sm->matches->elements[i]) {
      return 0;
    }
  }

  return 1;
}

/******************************************************************************* 
  The state of the polynucleotide repeat scanning loop
*/

typedef struct {
  Ring* last_bases;
  SubMatch** submatches;
  int submatch_index;
  int n_submatches;
  int length;
  char* chr;
  char* unit_buffer;
} SeqState;

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

  fprintf(stderr, "Will ignore subrepeats of the following lengths:");
  for (int i = 0; i < st->n_submatches; i++) {
    fprintf(stderr, " %i", st->submatches[i]->d);
  }
  fprintf(stderr, "\n");
   

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

/******************************************************************************* 
  Scan chromosome for polynucleotide repeats.

  Algorithm overview:

  Repeat tracking: Allocate an array of size 'r' to keep track of the last bases
  encountered. Also keep track of the length of a current repeat ('n'). If n <=
  'r', write the current base to the buffer. If 'n' > 'r' and the candidate
  repeat sequence is valid (see below), compare the current base to its relative
  modular position in the buffer. Keep scanning and incrementing n until a
  non-match is encountered. If the non-match is an N, start over completely,
  else write the new base to the buffer to make a new repeat target and reset
  'n' to 'r' or 'r' - 1 depending on if the new repeat is valid (see next).

  Subpattern filtering: We need to not print entries that are composed of
  repeats with smaller units. To identify what smaller unit sizes to find,
  determine the exclusive divisors of 'r' (eg divisors which do not divide other
  divisors of 'r'). For each divisor 'd', make a boolean array to represent all
  matches between all pairs in the current candidate repeat that are 'd' bp away
  from each other (without wraparound). When scanning through the sequence,
  track the match for each divisor between the current base and the base 'd' bp
  behind. A candidate repeat sequence will be 'valid' if none of these arrays
  are entirely true.

  NOTE: the arrays for the repeat buffer and subpattern matches are all
  represented as circular arrays that are indexed modulo the current position in
  the chromosome.
 */

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

void update_submatch_index (SeqState* st, int n) {
  SubMatch** sms = st->submatches;
  int x = st->submatch_index;

  if (x < st->n_submatches && sms[x]->d == n + 1) {
    st->submatch_index++;
  }
}

int next_n (SeqState* st, int p, int c) {
  write_ring(st->last_bases, p, c);
  update_matches(st, p);
  return st->last_bases->size - invalid_repeat(st);
}

void scan_seqN(FILE* fp, SeqState* st) {
  int c, c0;
  int p = 0;
  int n = 0;
  int r = st->last_bases->size;

  /* this should be the only thing we need to reset for each chromosome */
  st->submatch_index = 0;

  while (1) {
    c = fgetc(fp);

    /* ignore newlines */
    if (c != '\n') {

      if (c == EOF || c == HEADER_PREFIX) {
        /* ensure last repeat is printed if long enough */
        print_entryN(st, p, n);
        break;

      } else if (c == 'N') {
        print_entryN(st, p, n);
        n = 0;
        st->submatch_index = 0;

      } else if (n < r - 1) {
        write_ring(st->last_bases, p, c);
        update_matches(st, p);
        update_submatch_index(st, n);
        n++;

      } else if (n == r - 1) {
        n = next_n(st, p, c);

      } else {
        c0 = read_ring(st->last_bases, p);
        if (c0 == c) {
          n++;
          /*
            NOTE: this may seem wasteful since we are updating the submatch
            buffer(s) even when we know the repeat sequence is valid. However,
            the buffers for the subpatterns do not evenly divide the repeat
            sequence and thus will become out of phase as n increases. While it
            is possible to correct the phase in the alternative to this branch,
            this requires lots of complex math, which will lead to a larger
            slowdown given that this branch runs much less frequently (since
            repeats are relatively rare).
          */
          update_matches(st, p);
        } else {
          print_entryN(st, p, n);
          n = next_n(st, p, c);
        }
      }
      p++;
    }
  }
}

/*******************************************************************************
  Scan chromosome for homopolymers.

  Basic algorithm: Scan through each character and count the number of times
  it doesn't change. If it changes to a non-N and the length of the last
  repeated sequences is longer than our minimum, print it.
 */

void scan_seq1(FILE* fp, char* chr, int len) {
  char lastBase[2] = {'N', '\0'};
  int p = 0;
  int n = 1;
  int c;

  while (1) {
    c = fgetc(fp);

    /* ignore newlines */
    if (c != '\n') {

      if (c == EOF || c == HEADER_PREFIX) {
        /* ensure we print the last repeat if it is valid */
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

/*******************************************************************************
  FASTA parser
*/

int seek_char(FILE* fp, char t) {
  int c;

  do {
    c = fgetc(fp);
  } while (EOF != c && t != c);

  return c;
}

int parse_header(FILE* fp, char* chr) {
  /*
    ASSUME we are on a header (eg after a '>')
  */
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
  SeqState* st;
  int is_homopoly;

  if (len <= rep) {
    fprintf(stderr, "Repeat length must be less than total length\n");
    exit(-1);
  }

  if (fp == NULL) {
    perror("Error in opening file");
    return (-1);
  }

  /*
    Find the first header (which should be the first character, but genomics
    is weird...
  */
  seek_char(fp, HEADER_PREFIX);


  /* Initialize polynuc struct once if needed */
  is_homopoly = rep == 1;

  if (!is_homopoly) {
    fprintf(stderr, "Finding polynuc repeats >=%ibp with unit size %ibp\n", rep, len);
    st = init_seq_state(chr, rep, len);
  } else {
    fprintf(stderr, "Finding homopolymers >=%ibp\n", len);
  }

  printf("#repeat_length: %i\n", rep);
  printf("#total_length: %i\n", len);

  /* Parse each chromosome */
  while (parse_header(fp, chr) != EOF) {
    fprintf(stderr, "Parsing chromosome %s\n", chr);

    if (is_homopoly) {
      scan_seq1(fp, chr, len);
    } else {
      scan_seqN(fp, st);
    }
  }

  if (!is_homopoly) {
    free_seq_state(st);
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
