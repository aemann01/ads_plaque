#!/usr/bin/env python3

'''
Reads in a gzipped fastq file and discard id list, writes cleaned fastq file in gzip format -- this WILL NOT WORK WITH THE NOVASEQ DATA
Full file is read into memory (>64GB) while trying to write to file
'''

import gzip
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fastq', help="Input fastq file", required=True)
parser.add_argument('-i', '--ids', help="List of ids to remove", required=True)
parser.add_argument('-o', '--output', help="Output file name", required=True)
args = parser.parse_args()

from Bio import SeqIO, bgzf

discard = set(line.rstrip("\n").split(None,1)[0] for line in open(args.ids))
print("Found %i identifiers in %s to remove" % (len(discard), args.ids))

count = []

with gzip.open(args.fastq, "rt") as handle:
	records = list((r for r in SeqIO.parse(handle, "fastq") if r.id not in discard))
	for i in records:
		count.append(i.id)
with gzip.open(args.output, "wb") as outgz:
	SeqIO.write(records, outgz, "fastq")

print("Saved %i records from %s to %s" % (len(count), args.fastq, args.output))
if len(count) < len(discard): print("Warning %i IDs not found in %s" % (len(discard)-len(count), args.fastq))




