#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parse_variant_gff.py
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
#  


### Imported
import sys

### read files
gff_file = sys.argv[1]
start_pos = int(sys.argv[2]) # first exon starting position
strand = sys.argv[3] # + or - to indicate gene direction

#print "chromosome\tposition\tref\talt"
# get reverse complement sequence
def RC(seq):
	s1 = "BDHKMNRSVWYATGCbdhkmnrsvwyatgc"
	s2 = "VHDMKNYSBWRTACGvhdmknysbwrtacg"
	seq_dict = {s1[i]:s2[i] for i in range(len(s1))}
	return "".join([seq_dict[base] for base in reversed(seq)])

# get the exon location file
def parse_gff(infile, start_pos):
	with open(infile) as file_one:
		for line in file_one:
			if line.startswith("#"):
				continue
			line = line.strip()
			if line:
				ll = line.split("\t")
				chrom = ll[0]
				pos = abs(int(ll[3]) - start_pos) + 1
				details = ll[8].split(";")
				infor = dict(item.split("=") for item in details)
				ref = infor["reference_allele"]
				alt = infor["alternative_alleles"]
				print "\t".join([chrom, str(pos), ref, alt, strand])

parse_gff(gff_file, start_pos)









