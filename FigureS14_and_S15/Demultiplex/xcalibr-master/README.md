# XCALIBR

## Introduction

**xcalibr** - eXtracting Counting and LInking to Barcode References

**xcalibr** processes short reads and extracts a table of counts matched to index
sequences and tag sequences


## Method

Before the short reads are counted they are first processed into a binary
structure (hash), this ensures only unique reads processed which speeds up the
process.

**xcalibr** uses a template sequence to locate and extract the barcodes. Only exact matches to the reference sequence are currently supported. The reference sequences should be provided as fasta files. If no reference is present all barcodes at the requested location are reported which will result in large output files. You can use `xcalibr reduce` to reduce remove barcodes with low counts.
**xcalibr** supports three 

## Dependencies
The program is implemented in perl and currently uses only core modules.

## Command line tools

A wrapper bash script can be used to call the different tools in the package.
```
xcalibr hash
xcalibr analyze
xcalibr extract
xcalibr reduce
```

### hash
Run hash on a FastQ file. This file may be gzipped. Hash also reads (uncompressed FastQ) from STDIN if
no filename is provided.

### analyze
Use analyzehash to debug the input data. It will print a summary table given de
hash file and a range. You can use is to check constant parts or get a simple
overview when relatively few barcodes have been sequenced.

### extract
The workhorse program. Provide a template that xcalibr uses to locate and
extract the barcodes. `CGAT` are used to encode constant sequences and `X`, `Y`
and `Z` be be used as placeholder for the barcodes to be extracted. If only 1
or 2 placeholders are used a single table is generated. I three placeholders
are used one is used to generate multiple output files.

### reduce
Use 'xcalibr reduce' to filter the output table. Use the `--mincount` argument to specify the minimum number of reads that must be cocounted for a barcode to pass.

## Examples

If the xcalibr path is in your PATH you can use the following example:

```
cat testdata.fq | xcalibr hash testdata.bin
xcalibr analyze testdata.bin 7-22
xcalibr extract --template X6CAGGCGCTTAGGATCCY15  testdata.bin > res-all.txt
xcalibr extract --template X6CAGGCGCTTAGGATCCY15  testdata.bin | xcalibr reduce --mincount 5 > res-min5.txt
```
You can use zcat if your FastQ file is gzipped.
