#include <stdio.h>
#include <stdlib.h>

#define HEADER_PREFIX '>'
#define MAXREP 4
#define BED_FMT(f) "%s\t%i\t%i\tunit=%" #f "\n"

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

int invalid_repeat2 (int lb[MAXREP]) {
  /* Repeat is invalid if it is 2-mer homopolymer */
  return lb[0] == lb[1];
}

int invalid_repeat3 (int lb[MAXREP]) {
  /* Repeat is invalid if it is 3-mer homopolymer */
  return lb[0] == lb[1] && lb[1] == lb[2];
}

int invalid_repeat4 (int lb[MAXREP]) {
  /*
    Repeat is invalid if it has two identical 2-mers (which also excludes 4-mer
    homopolymers)
  */
  return lb[0] == lb[2] && lb[1] == lb[3];
}

void print_entryN (char* chr, int rep, int len, int p, int n, int lb[MAXREP]) {
  int o, j;

  if (n >= len) {
    char unit_buffer[MAXREP + 1];
    o = p - n;

    for (j = 0; j < rep; j++) {
      unit_buffer[j] = lb[(j + o) % rep];
    }
    unit_buffer[rep] = '\0';

    printf(BED_FMT(s), chr, p - n, p, unit_buffer);
  }
}

void scan_seqN (FILE* fp, char* chr, int rep, int len) {
  int c;
  int p = 0;
  int n = 0;
  int i = 0;
  int q = rep - 1;
  int last_bases[MAXREP];
  int (*invalid_repeat)(int lb[MAXREP]);

  switch (rep) {
  case 2:
    invalid_repeat = &invalid_repeat2;
    break;
  case 3: 
    invalid_repeat = &invalid_repeat3;
    break;
  case 4: 
    invalid_repeat = &invalid_repeat4;
    break;
  default: 
    fprintf(stderr, "invalid r (this should never happen)\n");
    exit(-1);
  }

  while (1) {
    c = fgetc(fp);

    /* ignore newlines */
    if (c != '\n') {

      if (c == 'N') {
        print_entryN(chr, rep, len, p, n, last_bases);
        n = 0;

      } else if (c == EOF || c == HEADER_PREFIX) {
        /* ensure last repeat is printed if long enough */
        print_entryN(chr, rep, len, i, n, last_bases);
        break;

      } else {
        if (n < q) {
          last_bases[i] = c;
          n++;

        } else if (n == q) {
          last_bases[i] = c;
          n = rep - invalid_repeat(last_bases);

        } else {
          if (last_bases[i] == c) {
            n++;
          } else {
            print_entryN(chr, rep, len, p, n, last_bases);
            last_bases[i] = c;
            n = rep - invalid_repeat(last_bases);
          }
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

void scan_seq1 (FILE* fp, char* chr, int len) {
  char last_base = 'N';
  int p = 0;
  int n = 1;
  int c;

  while (1) {
    c = fgetc(fp);

    /* ignore newlines */
    if (c != '\n') {
      if (c == last_base) {
        n++;

      } else {
        if (n >= len && last_base != 'N') {
          printf(BED_FMT(c), chr, p - n, p, last_base);
        }
        if (c == EOF || c == HEADER_PREFIX) {
          break;
        }
        n = 1;
        last_base = c;

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
  int is_homopoly;

  if (!(1 <= rep && rep <= MAXREP)) {
    fprintf(stderr, "Repeat length must be in [1,4]\n");
    exit(-1);
  }

  if (!(len >= rep)) {
    fprintf(stderr, "Total length must be >= repeat length\n");
    exit(-1);
  }

  if (fp == NULL) {
    fprintf(stderr, "Error opening file\n");
    exit(-1);
  }

  /*
    Find the first header (which should be the first character, but genomics
    is weird...
  */
  seek_char(fp, HEADER_PREFIX);

  is_homopoly = rep == 1;

  if (!is_homopoly) {
    fprintf(stderr, "Finding %i-mer repeats >=%ibp\n", rep, len);
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
