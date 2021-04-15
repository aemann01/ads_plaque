#!/usr/bin/env python3

'''
Reads in fastq file and discard id list, writes cleaned fastq file
'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fastq', help="Input fastq file", required=True)
parser.add_argument('-i', '--ids', help="List of ids to remove", required=True)
parser.add_argument('-o', '--output', help="Output file name", required=True)
args = parser.parse_args()

from Bio import SeqIO

discard = set(line.rstrip("\n").split(None,1)[0] for line in open(args.ids))
print("Found %i unique identifiers in %s to remove" % (len(discard), args.ids))
records = (r fo r in SeqIO.parse(args.fastq, "fastq") if r.id not in discard)
count = SeqIO.write(records, args.output, "fastq")
print("Saved %i records from %s to %s" % (count, args.fastq, args.output))
if count < len(discard): print("Warning %i IDs not found in %s" % (len(discard)-count, args.fastq))




