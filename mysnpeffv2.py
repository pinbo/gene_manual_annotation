#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  simple SNP effect prediction
#
#  Copyright 2018 Junli Zhang <zhjl86@gmail.com>
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
	s1 = "BDHKMNRSVWYATGCbdhkmnrsvwyatgc"
	s2 = "VHDMKNYSBWRTACGvhdmknysbwrtacg"
	seq_dict = {s1[i]:s2[i] for i in range(len(s1))}
	return "".join([seq_dict[base] for base in reversed(seq)])

# variant object
class exon(object):
	"""Information of exons"""
	def __init__(self):
		self.gene = ""
		#self.feature = ""
		self.start = 0 # 0 based
		self.end =  0 # 0 based

# get the exon location file
def get_annotation(infile):
	ann = {} # dict for variants
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line:
				anntemp = exon()
				ll = line.split()
				if ll[0] not in ann:
					ann[ll[0]] = []
				anntemp.gene = ll[0]
				anntemp.start = int(ll[2]) - 1 # 0 based
				anntemp.end = int(ll[3]) - 1 # 0 based
				ann[ll[0]].append(anntemp)
	return ann

# variant object
class var(object):
	"""Information of variants"""
	def __init__(self):
		self.gene = ""
		self.pos = 0 # 0 based
		self.ref = "" # 0 based
		self.alt = ""
		self.strand = ""
		self.eff = ""
		self.B62 = "" # blosum62 score

# get the variation file
# return a list of variant object
def get_var(infile):
	vart = {} # dict for variants
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line:
				vartemp = var()
				ll = line.split() # white space
				if ll[0] not in vart:
					vart[ll[0]] = []
				vartemp.gene = ll[0]
				vartemp.pos = int(ll[1]) - 1 # 0 based
				vartemp.strand = ll[4]
				if ll[4] == "+":
					vartemp.ref, vartemp.alt = ll[2:4]
				else:
					vartemp.ref = RC(ll[2])
					vartemp.alt = RC(ll[3])
				vart[ll[0]].append(vartemp)
	return vart

### read files
seq_file = sys.argv[1]
ann_file = sys.argv[2]
var_file = sys.argv[3]
#strand = sys.argv[4] # "+" for forward strand, "-" for reverse strand
outfile = sys.argv[4]

ff = get_fasta(seq_file)
ee_all = get_annotation(ann_file)
vv_all = get_var(var_file)

# sort the annotation values by exon position
for kk in ee_all:
	ll = ee_all[kk]
	#print kk, [x.start for x in ll]
	ll.sort(key=lambda x: x.start)

out = open(outfile, "w")

for kk in ff:
	ee = ee_all[kk] # annotation, exon positions
	vv = vv_all[kk] # variants
	## get the intron and splice range
	intron = [] # intron range. Here my intron did not include the two splice points
	splice = [] # splice point range
	exome = []
	i = -1 # in case there is only one exon
	for i in range(len(ee) - 1):
		splice += range(ee[i].end + 1, ee[i].end + 3)
		splice += range(ee[i + 1].start - 2, ee[i + 1].start)
		intron += range(ee[i].end + 3, ee[i+1].start - 2) # intron did not include the two splice points
		exome += range(ee[i].start, ee[i].end + 1)
	#print "i is, ", i

	exome += range(ee[i+1].start, ee[i+1].end + 1)

	#print "intron is ", intron
	#print "splice is ", splice
	#print "exon is ", exome
	## loop through the snps
	if vv[0].strand == "+":
		for j in vv: # loop through variants
			if len(j.ref) == 1 and len(j.alt) == 1: # snps
				if j.pos in intron:
					j.eff = "intron_variant"
				elif j.pos in splice:
					j.eff = "splice_site_mutant"
				elif j.pos in exome: # in exon
					ii = exome.index(j.pos) # position in the cDNA
					if ii%3 == 0:
						ref_codon = ff[j.gene][j.pos:(j.pos+3)]
						alt_codon = j.alt + ref_codon[1:]
					elif ii%3 == 1:
						ref_codon = ff[j.gene][(j.pos - 1):(j.pos+3-1)]
						alt_codon = ref_codon[0] + j.alt + ref_codon[2]
					else:
						ref_codon = ff[j.gene][(j.pos-2):(j.pos+3-2)]
						alt_codon = ref_codon[0:2] + j.alt
					# check AA change
					ref_AA = AA2[ref_codon] + "-" + AA3letter[AA2[ref_codon]]
					alt_AA = AA2[alt_codon] + "-" + AA3letter[AA2[alt_codon]]
					if ref_AA == alt_AA:
						j.eff = "synonymous_variant\t" + ref_AA + str(ii/3 + 1) + alt_AA
					elif alt_AA == "-":
						j.eff = "early_stop_condon\t" + ref_AA + str(ii/3 + 1) + alt_AA
					else:
						j.eff = "missense_variant\t" + ref_AA + str(ii/3 + 1) + alt_AA
						j.B62 = B62table[B62header.index(ref_AA[0])][B62header.index(alt_AA[0])]
				else:
					j.eff = "UTR_variant"
			## if it is an insertion. I suppose the first letter of the ref and alt are the same
			elif len(j.ref) == 1 and len(j.alt) > 1:
				if j.pos in intron:
					j.eff = "intron_variant"
				elif j.pos in splice:
					if j.pos - 1 in exome: # first letter in the left splice site
						if j.alt[0:2] == "GT":
							j.eff = "intron_variant_near_splice_site"
						else:
							j.eff = "splice_donor_mutant"
					elif j.pos + 1 in intron: # 2nd letter in the left splice site
						j.eff = "intron_variant_near_splice_site"
					elif j.pos - 1 in intron: # 1st letter in the right splice site
						if j.alt[0:2] == "AG":
							if (len(j.alt) - 1)%3 == 0:
								j.eff = "inframe_insertion_near_splice"
							else:
								j.eff = "frame_shift_near_splice"
						else:
							j.eff = "splice_acceptor_mutant"
					elif j.pos + 1 in exome: # 2nd letter in the right splice site
						if (len(j.alt) - 1)%3 == 0:
							j.eff = "inframe_insertion_near_splice"
						else:
							j.eff = "frame_shift_near_splice"
				elif j.pos in exome:# in exon
					if (len(j.alt) - 1)%3 == 0: # in frame
						j.eff = "amino_acid_insertion"
					else:
						j.eff = "frame_shift"
				else:
					j.eff = "UTR_variant"
			## if it is a deletion
			else:
				if set(range(j.pos, j.pos + len(j.ref))) <= set(intron): # if only in intron
					j.eff = "intron_variant"
				elif set(range(j.pos, j.pos + len(j.ref))) <= set(exome): # if only in exome
					if (len(j.ref) - 1)%3 == 0: # in frame
						j.eff = "amino_acid_deletion"
					else:
						j.eff = "frame_shift"
				elif set(range(j.pos, j.pos + len(j.ref))) & set(splice): # affecting splice
					j.eff = "splice_site_mutant"
				else:
					j.eff = "UTR_variant"
	else: # reverse strand
		for j in vv:
			if len(j.ref) == 1 and len(j.alt) == 1: # snps
				if j.pos in intron:
					j.eff = "intron_variant"
				elif j.pos in splice:
					j.eff = "splice_site_mutant"
				elif j.pos in exome: # in exon
					ii = exome.index(j.pos) # position in the cDNA
					if ii%3 == 0:
						ref_codon = ff[j.gene][j.pos:(j.pos+3)]
						alt_codon = j.alt + ref_codon[1:]
					elif ii%3 == 1:
						ref_codon = ff[j.gene][(j.pos - 1):(j.pos+3-1)]
						alt_codon = ref_codon[0] + j.alt + ref_codon[2]
					else:
						ref_codon = ff[j.gene][(j.pos-2):(j.pos+3-2)]
						alt_codon = ref_codon[0:2] + j.alt
					# check AA change
					ref_AA = AA2[ref_codon] + "-" + AA3letter[AA2[ref_codon]]
					alt_AA = AA2[alt_codon] + "-" + AA3letter[AA2[alt_codon]]
					if ref_AA == alt_AA:
						j.eff = "synonymous_variant\t" + ref_AA + str(ii/3 + 1) + alt_AA
					elif alt_AA == "-":
						j.eff = "early_stop_condon\t" + ref_AA + str(ii/3 + 1) + "*"
					else:
						j.eff = "missense_variant\t" + ref_AA + str(ii/3 + 1) + alt_AA
						j.B62 = B62table[B62header.index(ref_AA[0])][B62header.index(alt_AA[0])]
				else:
					j.eff = "UTR_variant"
			## if it is an insertion. I suppose the first letter of the ref and alt are the same
			elif len(j.ref) == 1 and len(j.alt) > 1:
				if j.pos in intron:
					j.eff = "intron_variant"
				elif j.pos in splice:
					if j.pos - 1 in exome: # first letter in the left splice site
						if (len(j.alt) - 1)%3 == 0:
							j.eff = "inframe_insertion_near_splice"
						else:
							j.eff = "frame_shift_near_splice"
					elif j.pos + 1 in intron: # 2nd letter in the left splice site
						if j.alt[-2:] == "GT":
							if (len(j.alt) - 1)%3 == 0:
								j.eff = "inframe_insertion_near_splice"
							else:
								j.eff = "frame_shift_near_splice"
						else:
							j.eff = "splice_donor_mutant"
					elif j.pos - 1 in intron: # 1st letter in the right splice site
						j.eff = "intron_variant_near_splice_site"
					elif j.pos + 1 in exome: # 2nd letter in the right splice site
						if j.alt[-2:] == "AG":
							j.eff = "intron_variant_near_splice_site"
						else:
							j.eff = "splice_acceptor_mutant"
				elif j.pos in exome:# in exon
					if (len(j.alt) - 1)%3 == 0: # in frame
						j.eff = "inframe_insertion"
					else:
						j.eff = "frame_shift_insertion"
				else:
					j.eff = "UTR_variant"
			## if it is a deletion
			else:
				if set(range(j.pos - len(j.ref) + 1, j.pos+1)) <= set(intron): # if only in intron
					j.eff = "intron_variant"
				elif set(range(j.pos - len(j.ref) + 1, j.pos+1)) <= set(exome): # if only in exome
					if (len(j.ref) - 1)%3 == 0: # in frame
						j.eff = "inframe_deletion"
					else:
						j.eff = "frame_shift_deletion"
				elif set(range(j.pos - len(j.ref) + 1, j.pos+1)) & set(splice): # affecting splice
					j.eff = "splice_acceptor_mutant"
				else:
					j.eff = "UTR_variant"

	# write out the effect
	out.write("Gene\tPosition\tRef_allele\tAlt_allele\tEffect\tAA_change\tBLOSUM62_score\n")
	for k in vv:
		out.write("\t".join([k.gene, str(k.pos + 1), k.ref, k.alt, k.eff, str(k.B62)]) + "\n")

	# get the cDNA
	gene = vv[0].gene
	cdna = "".join([ff[gene][i] for i in exome])
	out.write("\n>" + gene + "_cDNA\n" + cdna + "\n\n")

	n = len(cdna)/3 # number of amino acid
	protein = ""
	for i in range(int(n)):
		codon = cdna[(i*3):(i*3+3)]
		aa = AA2[codon]
		protein += aa
	out.write(">" + gene + "_protein\n" + protein + "\n\n")

