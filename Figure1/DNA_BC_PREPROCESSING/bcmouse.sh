# get the fastq files or symlink to target data
find /path/to/fastq  -name "*.fastq.gz" -exec ln -s {} \;

# we assume the barcode is present in the filename. Normal demultiplexing allows 1bp distance,
# but we collect only the perfect matches with the index sequences.
# these reads are piped through 'xcalibr hash' to create the analysis data
ls -1 *_R1_001.fastq.gz |sed -e 'p; s/_S[0-9].*_R1_001.fastq.gz//p; s/.*_//  ' | parallel -N3 -j8 'zcat {1} | perl ./perfect.pl {3} | xcalibr hash {2}-perfect.bin'

# we test the hash for the expected constant part. The illumina sequencer tends to make mistakes here. We
# collect the top hits to create a constant template sequence that discards as little data as possible
ls -1 *perfect.bin |parallel -j8 'xcalibr analyze {} 13-38 > {.}-constant.txt'
cat *constant.txt |grep "^C.CG" |cut -f1 | sort | uniq -c | sort -n | tail -n20
#67		CTCGAGGCCATCGAAGTATCAAGTCC
#112	CTCGAGGTCATCGAAGAATCAAGTCC
#131	CNCGAGGNCATCGAAGTATCAAGTNC
#131	CNCGAGGNCATCGAAGTATCAAGTCC
#132	CTCGAGGTCATCGAAGTATCAAGTCC
#132	CTCGAGGNCATCGAAGTATCAAGTCC

# use constant: CTCGAGGNCATCGAAGTATCAAG
# note the N's in the variable positions

# extract a count table from the hash using the template. We use the sparse output to avoid creating a huge matrix with al lot of 0's
ls -1 *perfect.bin |parallel -j12 'xcalibr extract --template X12CTCGAGGNCATCGAAGTATCAAGYn --rows Y --cols X --requireseq --sparse --output {.}-long.txt {} &> {.}.log'


# apply post-processing to the sparse count table.
# we expect a constant sequence (template2.fa) to start somewhere in every sequence
# we use ncbi blast2 to find these hits
# we hope to find only one hit per sequence, but this is not always the case.
# a perl script uses the blast results to process the long output
# The long tables and the *clipped-template2-e6.txt are sent to Joost Beltman/Jos Urbanus

for F in *perfect-long.txt
do
BASE=${F%.txt}
#Blast against template clip tags
cut -f1,2 --output-delimiter=Y $F |perl count2fasta.pl >$BASE.fa
LD_LIBRARY_PATH=./blast2/lib/ ./blast2/bin/blast2 -p blastn -i $BASE.fa -j template2.fa  -W7 -K1  -e 1e-6 -m8 >$BASE-blast-template2-e6.txt
#check counts
echo "Tailing $BASE, should be 1 count:"
cut -f1 $BASE-blast-template2-e6.txt | sort | uniq -c | sort -n | tail -5
#cut tags
perl cuttagslong.pl $F $BASE-blast-template2-e6.txt >$BASE-clipped-template2-e6.txt
done




