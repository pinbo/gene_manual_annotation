# prepare variant file from the downloaded gff3 file from Jbrowse
 ../parse_variant_gff.py 1.\ All\ Accessions\ with\ SNPEff-chr1A-468604981..468608315.gff3 468607974 - > variants2.txt

# remember to change the gene names
../mysnpeffv2.py gene.fa annotation.txt variants.txt > SNP_effects.txt

