#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
#  

# extract sequences and gff3 files from the whole set and run mysnpeffv3.py
# example: run_getkasp.py for_polymarker.csv 3 200 1 1 1 65

#########################

from subprocess import call
import sys

# variant object
class region(): # for each chromosome
	"""Information of exons"""
	def __init__(self):
		self.chr = ""
		self.min = 0 # start 
		self.max =  0 # end

# get regions
# return a list of region object
def get_region(infile):
	regionlist = [] # list for regions
	posdict = {} # chromsome with a list of positions
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line:
				ll = line.split() # white space
				chr, pos = ll[0:2]
				if chr in posdict:
					posdict[chr].append(pos)
				else:
					posdict[chr] = [pos]
	# put chrom, start, end into a list.
	for k, v in posdict.items():
		rr = region()
		rr.chr = k
		rr.min = min(v)
		rr.max = max(v)
		regionlist.append(rr)
	return regionlist

# get min and max from the gff3 subset
def getgff3range(gff3file): # suppose only for one chromsome
	col4 = []
	col5 = []
	with open(gff3file) as file_one:
		for line in file_one:
			if line. startswith("#"):
				continue
			line = line.strip()
			if line:
				ll = line.split()
				col4.append(ll[3])
				col5.append(ll[4])
	return min(col4), max(col5)

def main(args):
	script_path = os.path.dirname(os.path.realpath(__file__))
	reference, gff3file, snpfile = args[1:4]
	# step 1: get the flanking sequences
	#cmd = "blastdbcmd -entry_batch " + infile + " -db " + reference + " > " + flanking_file
	#call(cmd, shell=True)
	# step 1: get chromsome and position information from the snpfile
	regionlist = get_region(snpfile)
	# extract subset of gff3
	for rr in regionlist:
		gff3subset = rr.chr + ".subset.gff3"
		chrsubset = rr.chr + ".subset.fasta"
		snpsubset = snpfile + ".subset.txt"
		output = "out." + rr.chr + ".txt"
		# awk '$1=="4B"{$5 >= 1 && $4 <= 10}' input > output
		cmd1 = "gawk 'BEGIN{mm=" + str(rr.max) + "} $1==\"" + rr.chr + "\" && $5>=" + str(rr.min) + " && $4<=mm{print; if ($5>mm) mm=$5}' " + gff3file + " > " + gff3subset
		print(cmd1)
		call(cmd1, shell=True)
		# update rr.min and rr.max from the gff3 subset file
		newmin, newmax = getgff3range(gff3subset)
		if newmin < rr.min:
			rr.min = newmin
		if newmax > rr.max:
			rr.max = newmax
		cmd2 = "blastdbcmd" + " -db '" + reference + "' -entry " + rr.chr + " -range " + str(rr.min) + "-" + str(rr.max) + " > " + chrsubset
		print(cmd2)
		call(cmd2, shell=True)
		# subset SNP input file
		cmd3 = "gawk '$1==\"" + rr.chr + "\"' " + snpfile + " > " + snpsubset
		print(cmd3)
		call(cmd3, shell=True)
		# run mysnpeffv3
		cmd4 = script_path + "/mysnpeffv3.py " + chrsubset + " " + gff3subset + " " + snpsubset + " " + output + " " + str(rr.min)
		print(cmd4)
		call(cmd4, shell=True)
	
	# merge all the results
	cmd4 = "cat out.*.txt > all.out.txt"
	call(cmd4, shell=True)
	return 0

if __name__ == '__main__':
	from subprocess import call
	import sys, os
	sys.exit(main(sys.argv))
