''' Adds missing fields to list of taxonomy strings. Use: python fix_taxonomy.py my_taxa.txt > new_taxa.txt '''

import sys

with open(sys.argv[1], "r") as f:
	for line in f:
		x = line.rstrip("\n")
		y = x.count(";")
		if y < 10:
			add = 10 - y
			new = (";" + x.split(";")[-1] + "_unknown") * add
			print(x + new)
		else:
			print(x)
