# Input include: subset of a chromosome sequence, gff3 of gene annotation for within subset region, SNPs (chrom, pos, ref, alt. NO header), output file name, the start region the subset.

python ../mysnpeffv3.py chr4B-563625250..564907087.fasta HighConfidenceJuly2017-chr4B-563625250..564907087.gff3 SNPs.txt out.txt 563625250