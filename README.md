# Simple Repeat Finder

Find repeated nucleotide sequences in a FASTA file.

## Building

Ensure gcc is installed:

```
make
```

## Usage

To find repeats of at least `L` length with unit size `R` in FASTA file
`INPUTFILE`:

```
./repseq R L INPUTFILE
```

Detected regions will be printed to stdout in .bed format with an additional
fourth column specifying the unit sequence of the repeat.

NOTE: for a given `R` this will not find any repeats that could be detected by a
divisor of `R`. For example, if `R` is 4 it will not detect and print a
repeated sequence like `ATAT` since this is composed of two smaller repeats
(`AT`) which correspond to `R` of 2.
