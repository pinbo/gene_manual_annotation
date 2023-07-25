#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  simple SNP effect prediction
#
#  Copyright 2020 Junli Zhang <zhjl86@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#

# function: SNP/indel effect prediction with input downloaded from JBrowse (gff3 gene annotation file), SNP file, and subset of the reference genome.


### Imported
import sys

AAcodon = {
	"I" : ["ATT", "ATC", "ATA"],
	"L" : ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
	"V" : ["GTT", "GTC", "GTA", "GTG"],
	"F" : ["TTT", "TTC"],
	"M" : ["ATG"],
	"C" : ["TGT", "TGC"],
	"A" : ["GCT", "GCC", "GCA", "GCG"],
	"G" : ["GGT", "GGC", "GGA", "GGG"],
	"P" : ["CCT", "CCC", "CCA", "CCG"],
	"T" : ["ACT", "ACC", "ACA", "ACG"],
	"S" : ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
	"Y" : ["TAT", "TAC"],
	"W" : ["TGG"],
	"Q" : ["CAA", "CAG"],
	"N" : ["AAT", "AAC"],
	"H" : ["CAT", "CAC"],
	"E" : ["GAA", "GAG"],
	"D" : ["GAT", "GAC"],
	"K" : ["AAA", "AAG"],
	"R" : ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
	"*" : ["TAA", "TAG", "TGA"]
}

AA3letter = {
	"I" : "Ile",
	"L" : "Leu",
	"V" : "Val",
	"F" : "Phe",
	"M" : "Met",
	"C" : "Cys",
	"A" : "Ala",
	"G" : "Gly",
	"P" : "Pro",
	"T" : "Thr",
	"S" : "Ser",
	"Y" : "Tyr",
	"W" : "Trp",
	"Q" : "Gln",
	"N" : "Asn",
	"H" : "His",
	"E" : "Glu",
	"D" : "Asp",
	"K" : "Lys",
	"R" : "Arg",
	"*" : "Stop"
}

# make a new codon table with 3-letter as keys
AA2 = {}
for k, v in AAcodon.items():
	for i in v:
		AA2[i] = k


# BLOSUM62
B62header = ["C", "S", "T", "A", "G", "P", "D", "E", "Q", "N", "H", "R", "K", "M", "I", "L", "V", "W", "Y", "F", "*"]
B62table = [[9, -1, -1, 0, -3, -3, -3, -4, -3, -3, -3, -3, -3, -1, -1, -1, -1, -2, -2, -2, -4],
[-1, 4, 1, 1, 0, -1, 0, 0, 0, 1, -1, -1, 0, -1, -2, -2, -2, -3, -2, -2, -4],
[-1, 1, 5, 0, -2, -1, -1, -1, -1, 0, -2, -1, -1, -1, -1, -1, 0, -2, -2, -2, -4],
[0, 1, 0, 4, 0, -1, -2, -1, -1, -2, -2, -1, -1, -1, -1, -1, 0, -3, -2, -2, -4],
[-3, 0, -2, 0, 6, -2, -1, -2, -2, 0, -2, -2, -2, -3, -4, -4, -3, -2, -3, -3, -4],
[-3, -1, -1, -1, -2, 7, -1, -1, -1, -2, -2, -2, -1, -2, -3, -3, -2, -4, -3, -4, -4],
[-3, 0, -1, -2, -1, -1, 6, 2, 0, 1, -1, -2, -1, -3, -3, -4, -3, -4, -3, -3, -4],
[-4, 0, -1, -1, -2, -1, 2, 5, 2, 0, 0, 0, 1, -2, -3, -3, -2, -3, -2, -3, -4],
[-3, 0, -1, -1, -2, -1, 0, 2, 5, 0, 0, 1, 1, 0, -3, -2, -2, -2, -1, -3, -4],
[-3, 1, 0, -2, 0, -2, 1, 0, 0, 6, 1, 0, 0, -2, -3, -3, -3, -4, -2, -3, -4],
[-3, -1, -2, -2, -2, -2, -1, 0, 0, 1, 8, 0, -1, -2, -3, -3, -3, -2, 2, -1, -4],
[-3, -1, -1, -1, -2, -2, -2, 0, 1, 0, 0, 5, 2, -1, -3, -2, -3, -3, -2, -3, -4],
[-3, 0, -1, -1, -2, -1, -1, 1, 1, 0, -1, 2, 5, -1, -3, -2, -2, -3, -2, -3, -4],
[-1, -1, -1, -1, -3, -2, -3, -2, 0, -2, -2, -1, -1, 5, 1, 2, 1, -1, -1, 0, -4],
[-1, -2, -1, -1, -4, -3, -3, -3, -3, -3, -3, -3, -3, 1, 4, 2, 3, -3, -1, 0, -4],
[-1, -2, -1, -1, -4, -3, -4, -3, -2, -3, -3, -2, -2, 2, 2, 4, 1, -2, -1, 0, -4],
[-1, -2, 0, 0, -3, -2, -3, -2, -2, -3, -3, -3, -2, 1, 3, 1, 4, -3, -1, -1, -4],
[-2, -3, -2, -3, -2, -4, -4, -3, -2, -4, -2, -3, -3, -1, -3, -2, -3, 11, 2, 1, -4],
[-2, -2, -2, -2, -3, -3, -3, -2, -1, -2, 2, -2, -2, -1, -1, -1, -1, 2, 7, 3, -4],
[-2, -2, -2, -2, -3, -4, -3, -3, -3, -3, -1, -3, -3, 0, 0, 0, -1, 1, 3, 6, -4],
[-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]]


# process
# function to extract sequences from a fasta file 
def get_fasta(infile):
	fasta = {} # dictionary for alignment
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line: # skip blank lines
				if line.startswith(">"):
					sequence_name = line.lstrip("> ").split()[0] # left strip > or space, so " > abc edf" will be "abc edf", then split by space to get "abc"
					fasta[sequence_name] = ""
				else:
					fasta[sequence_name] += line.replace(" ", "") # remove spaces in case
	return fasta

# get reverse complement sequence
def RC(seq):
	s1 = "BDHKMNRSVWYATGCbdhkmnrsvwyatgc-"
	s2 = "VHDMKNYSBWRTACGvhdmknysbwrtacg-"
	seq_dict = {s1[i]:s2[i] for i in range(len(s1))}
	return "".join([seq_dict[base] for base in reversed(seq)])

# variant object
class exon(object):
	"""Information of exons"""
	def __init__(self):
		self.gene = ""
		self.start = 0 # 0 based
		self.end =  0 # 0 based

# variant object
class mRNA(object):
	"""Information of exons"""
	def __init__(self):
		self.gene = ""
		self.strand = "+"
		self.start = 0
		self.end =  0
		self.seq = ""
		self.intron = []
		self.exome = []
		self.splice = []

# process gff3 file downloaded from JBrowse
def parse_gff3(infile):
	ann = {} # dict for CDS
	genelist = [] # ordered gene lists
	genes = {} # mRNA dict
	mRNA_name = ""
	with open(infile) as file_one:
		for line in file_one:
			if line. startswith("#"):
				continue
			line = line.strip()
			if line:
				ll = line.split()
				if ll[2] == "mRNA":
					mm = mRNA()
					mm.strand = ll[6]
					mm.start = int(ll[3]) # in RefSeq coordinate
					mm.end = int(ll[4]) # 0 based, but in RefSeq coordinate
					tt = ll[8].split(";")
					tempdict = {}
					for ii in tt:
						key, val = ii.split("=")
						tempdict[key] = val
					mRNA_name = tempdict["ID"]
					genelist.append(mRNA_name)
					mm.gene = mRNA_name
					ann[mRNA_name] = []
					genes[mRNA_name] = mm
				elif ll[2] == "CDS":
					anntemp = exon()
					anntemp.gene = mRNA_name
					anntemp.start = int(ll[3])
					anntemp.end = int(ll[4])
					ann[mRNA_name].append(anntemp)
	return ann, genes, genelist # CDS and mRNAs

# variant object: snps and indels
class var(object):
	"""Information of variants"""
	def __init__(self):
		self.gene = ["intergene"] # a list, the mRNAs it is located
		self.pos = 0
		self.ref = "" 
		self.alt = ""
		self.strand = ["+"]
		self.eff = ["intergene"]
		self.B62 = [""] # blosum62 score

# get the variation file
# return a list of variant object
def get_var(infile):
	vart = [] # list for variants
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line:
				vartemp = var()
				ll = line.split() # white space
				vartemp.pos = int(ll[1])
				vartemp.ref, vartemp.alt = ll[2:4]
				vart.append(vartemp)
	return vart

### read files
seq_file = sys.argv[1] # subset of one chromosome, for example chr4B:563784341..563787208
gff3_file = sys.argv[2] # gene annnotation file downloaded from JBrowse from the same region, such as chr4B:563784341..563787208
var_file = sys.argv[3] # a SNP file with 4 columns (no header): chromosome, position, ref, alt
outfile = sys.argv[4] # output file name
seq_start = int(sys.argv[5]) # an integar of the star position of sequences in seq_file, for example 563784341

allseq = get_fasta(seq_file)
print("sequence was parsed")
refseq = allseq[list(allseq.keys())[0]]
print(refseq[0:10])

#ee_all = get_annotation(ann_file)
ee_all, genes, genelist = parse_gff3(gff3_file) # CDSs and mRNAs
snplist = get_var(var_file) # SNPs and indels
# get mRNA sequences
for v in genes.values():
	v.seq = refseq[(v.start - seq_start):(v.end - seq_start + 1)]
	if v.strand == "-":
		v.seq = RC(v.seq)

# find the gene locations of each snp and its strand
for snp in snplist:
	#print "snp is ", snp.gene
	for gene in genelist:
		#print "Gene is", gene, snp.pos
		#print "Gene positions", genes[gene].start, genes[gene].end
		if snp.pos >= genes[gene].start and snp.pos <= genes[gene].end:
			snp.gene.append(gene)
			snp.strand.append(genes[gene].strand)
			#break

# sort exons for each mRNA splice
for kk in ee_all:
	ee = ee_all[kk] # exons
	gene = genes[kk] # mRNA
	i = -1 # in case there is only one exon
	for i in range(len(ee) - 1):
		gene.splice += range(ee[i].end + 1, ee[i].end + 3)
		gene.splice += range(ee[i + 1].start - 2, ee[i + 1].start)
		gene.intron += range(ee[i].end + 3, ee[i+1].start - 2) # intron did not include the two splice points
		gene.exome += range(ee[i].start, ee[i].end + 1)
	#print "i is, ", i
	gene.exome += range(ee[i+1].start, ee[i+1].end + 1)
	if gene.strand == "-":
		gene.exome.sort(reverse = True) # so the biggest number wil be the first one
		gene.start, gene.end = gene.end, gene.start

# write output
out = open(outfile, "w")
out.write("Gene\tPosition\tStrand\tRef_allele\tAlt_allele\tEffect\tAA_change\tBLOSUM62_score\n")

for j in snplist:
	nn = len(j.gene) # number of mRNAs
	if nn == 1: # only default "intergene"
		j.eff = "InterGene"
		out.write("\t".join([j.gene[0], str(j.pos), j.strand[0], j.ref, j.alt, j.eff[0], str(j.B62[0])]) + "\n")
		continue
	for i in range(1, nn):
		ID = j.gene[i] # mRNA ID
		mRNA = genes[ID]
		seq = mRNA.seq
		ee = ee_all[ID] # annotation, exon positions
		intron = mRNA.intron
		splice = mRNA.splice
		exome = mRNA.exome
		alt = j.alt # alternative allele
		B62 = ""
		eff = ""
		if j.strand[i] == "-":
			#ref = RC(j.ref)
			alt = RC(j.alt)
		if len(j.ref) == 1 and len(j.alt) == 1: # snps
			if j.pos in intron:
				eff = "intron_variant"
			elif j.pos in splice:
				eff = "splice_site_mutant"
			elif j.pos in exome: # in exon
				ii = exome.index(j.pos) # position in the cDNA
				relative_pos = abs(j.pos - mRNA.start) # position in the mRNA
				if ii%3 == 0:
					ref_codon = seq[relative_pos:(relative_pos+3)]
					alt_codon = alt + ref_codon[1:]
				elif ii%3 == 1:
					ref_codon = seq[(relative_pos - 1):(relative_pos+3-1)]
					alt_codon = ref_codon[0] + alt + ref_codon[2]
				else:
					ref_codon = seq[(relative_pos-2):(relative_pos+3-2)]
					alt_codon = ref_codon[0:2] + alt
				# check AA change
				ref_AA = AA2[ref_codon] + "-" + AA3letter[AA2[ref_codon]]
				alt_AA = AA2[alt_codon] + "-" + AA3letter[AA2[alt_codon]]
				if ref_AA == alt_AA:
					eff = "synonymous_variant\t" + ref_AA + str(int(ii/3) + 1) + alt_AA
				elif alt_AA == "-":
					eff = "early_stop_condon\t" + ref_AA + str(int(ii/3) + 1) + alt_AA
				else:
					eff = "missense_variant\t" + ref_AA + str(int(ii/3) + 1) + alt_AA
					B62 = B62table[B62header.index(ref_AA[0])][B62header.index(alt_AA[0])]
			else:
				eff = "UTR_variant"
		## if it is an insertion. I suppose the first letter of the ref and alt are the same
		elif len(j.ref) == 1 and len(j.alt) > 1:
			if j.pos in intron:
				eff = "intron_variant"
			elif j.pos in splice:
				if j.pos - 1 in exome: # first letter in the left splice site
					if alt[0:2] == "GT":
						eff = "intron_variant_near_splice_site"
					else:
						eff = "splice_donor_mutant"
				elif j.pos + 1 in intron: # 2nd letter in the left splice site
					eff = "intron_variant_near_splice_site"
				elif j.pos - 1 in intron: # 1st letter in the right splice site
					if alt[0:2] == "AG":
						if (len(j.alt) - 1)%3 == 0:
							eff = "inframe_insertion_near_splice"
						else:
							eff = "frame_shift_near_splice"
					else:
						eff = "splice_acceptor_mutant"
				elif j.pos + 1 in exome: # 2nd letter in the right splice site
					if (len(j.alt) - 1)%3 == 0:
						eff = "inframe_insertion_near_splice"
					else:
						eff = "frame_shift_near_splice"
			elif j.pos in exome:# in exon
				if (len(j.alt) - 1)%3 == 0: # in frame
					eff = "amino_acid_insertion"
				else:
					eff = "frame_shift"
			else:
				eff = "UTR_variant"
		## if it is a deletion
		else:
			if set(range(j.pos, j.pos + len(j.ref))) <= set(intron): # if only in intron
				eff = "intron_variant"
			elif set(range(j.pos, j.pos + len(j.ref))) <= set(exome): # if only in exome
				if (len(j.ref) - 1)%3 == 0: # in frame
					eff = "amino_acid_deletion"
				else:
					eff = "frame_shift"
			elif set(range(j.pos, j.pos + len(j.ref))) & set(splice): # affecting splice
				eff = "splice_site_mutant"
			else:
				eff = "UTR_variant"
		# assign eff to j.eff
		j.eff.append(eff)
		# write out the effect
		out.write("\t".join([j.gene[i], str(j.pos), j.strand[i], j.ref, j.alt, eff, str(B62)]) + "\n")

# get the cDNA
""" gene = vv[0].gene
cdna = "".join([ff[gene][i] for i in exome])
out.write("\n>" + gene + "_cDNA\n" + cdna + "\n\n")

n = len(cdna)/3 # number of amino acid
protein = ""
for i in range(n):
	codon = cdna[(i*3):(i*3+3)]
	aa = AA2[codon]
	protein += aa
out.write(">" + gene + "_protein\n" + protein + "\n\n") """

