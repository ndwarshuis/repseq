#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEADER_PREFIX '>'

/******************************************************************************* 
  Various helper functions
*/

void print_entry(char* chr, int p, char n, char* unit) {
  printf("%s\t%i\t%i\tunit=%s\n", chr, p - n, p, unit);
}

/******************************************************************************* 
  The state of the polynucleotide repeat scanning loop
*/

typedef struct SeqState {
  char* last_bases;
  int length;
  int repsize;
  int (*invalid_repeat)(struct SeqState *st);
  char* chr;
  char* unit_buffer;
} SeqState;

int invalid_repeat2 (SeqState* st) {
  char *lb;

  lb = st->last_bases;

  /* Repeat is invalid if it is 2-mer homopolymer */
  return lb[0] == lb[1];
}

int invalid_repeat3 (SeqState* st) {
  char *lb;

  lb = st->last_bases;

  /* Repeat is invalid if it is 3-mer homopolymer */
  return lb[0] == lb[1] && lb[1] == lb[2];
}

int invalid_repeat4 (SeqState* st) {
  char *lb;

  lb = st->last_bases;

  /*
    Repeat is invalid if it has two identical dinuc repeats (which will also
    exclude 4-mer homopolymers)
  */
  return lb[0] == lb[2] && lb[1] == lb[3];
}

SeqState* init_seq_state (char* chr, int rep, int len) {
  SeqState* st;

  st = malloc(sizeof(*st));

  st->last_bases = malloc(rep * sizeof(char));
  st->length = len;
  st->repsize = rep;
  st->chr = chr;
  st->unit_buffer = malloc((rep + 1) * sizeof(char));
  st->unit_buffer[rep] = '\0';

  switch (rep) {
    case 2:
      st->invalid_repeat = &invalid_repeat2;
      break;
    case 3: 
      st->invalid_repeat = &invalid_repeat3;
      break;
    case 4: 
      st->invalid_repeat = &invalid_repeat4;
      break;
  default: 
    fprintf(stderr, "invalid r (this should never happen)\n");
    exit(-1);
  }

  return st;
}

void free_seq_state (SeqState* st) {
  free(st->last_bases);
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

  Subpattern filtering: Just finding a repeat is not enough; we need to
  determine if the repeat has any "smaller" repeats inside it to consider it
  valid. For now, repeats of length 2, 3, and 4 are hardcoded for speed. 2 and
  3-mer repeats are checked to see if they are homopolymers (eg all bases the
  same) and 4-mers are checked for 2 identical 2-mer (which also excludes
  homopolymers).
 */

void print_entryN (SeqState* st, int p, int n) {
  int o, j;

  if (n >= st->length) {
    o = p - n;

    for (j = 0; j < st->repsize; j++) {
      st->unit_buffer[j] = st->last_bases[(j + o) % st->repsize];
    }

    print_entry(st->chr, p, n, st->unit_buffer);
  }
}

int next_n (SeqState *st, int i, int c) {
  st->last_bases[i] = c;
  return st->repsize - st->invalid_repeat(st);
}

void scan_seqN (FILE* fp, SeqState* st) {
  int c;
  int p = 0;
  int n = 0;
  int i = 0;
  int q = st->repsize - 1;

  while (1) {
    c = fgetc(fp);

    /* ignore newlines */
    if (c != '\n') {

      if (c == 'N') {
        print_entryN(st, p, n);
        n = 0;

      } else if (c == EOF || c == HEADER_PREFIX) {
        /* ensure last repeat is printed if long enough */
        print_entryN(st, p, n);
        break;

      } else {
        if (n < q) {
          st->last_bases[i] = c;
          n++;

        } else if (n > q) {
          if (st->last_bases[i] == c) {
            n++;
          } else {
            print_entryN(st, p, n);
            n = next_n(st, i, c);
          }

        } else {
          st->last_bases[i] = c;
          n = next_n(st, i, c);
        }
      }
      p++;
      /* super fast 'modulo operator' */
      if (i == q) {
        i = 0;
      } else {
        i++;
      }
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
  char last_base[2] = {'N', '\0'};
  int p = 0;
  int n = 1;
  int c;

  while (1) {
    c = fgetc(fp);

    /* ignore newlines */
    if (c != '\n') {
      if (c == last_base[0]) {
        n++;

      } else {
        if (n >= len && last_base[0] != 'N') {
          print_entry(chr, p, n, last_base);
        }

        if (c == EOF || c == HEADER_PREFIX) {
          break;
        } else {
          n = 1;
          last_base[0] = c;
        }
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

  if (4 < rep) {
    fprintf(stderr, "Repeat length must be in [1,4]\n");
    exit(-1);
  }

  if (fp == NULL) {
    fprintf(stderr, "Error in opening file\n");
    exit(-1);
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
