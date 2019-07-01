#!/usr/bin/env python

""" Select only those lines from file 2, if value x in column y
matches value a in column b of file 1. """

import argparse

def filter_file(f1, f2, s1="\t", s2="\t", X=0, Y=0):
	
	# Get filter from file 1
	fltr = set([ line.rstrip().split(s1)[X] for line in open(f1) ])

	# Output those lines in file 2 for which the value in column Y
	# matches with the filter from file 1. 
	for line in open(f2):
		if line.rstrip().split('\t')[Y] in fltr:
			print(line.rstrip())

if __name__ == '__main__':

	def zero_or_more(v):
		if int(v) < 0:
			raise argparse.ArgumentTypeError("--X and --Y can't be lower than 0.")
		return int(v)

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f1', '--file1', required=True, help="File 1, the file with the final selection in column X.")
	parser.add_argument('-f2', '--file2', required=True, help="File 2, the file to select from based on match with file 1.")
	parser.add_argument('--s1', default="\t", help="File separator for file 1.")
	parser.add_argument('--s2', default="\t", help="File separator for file 2.")
	parser.add_argument('--X', type=zero_or_more, default=0, help="Column X in file 1, containing the values to match in file 2.")
	parser.add_argument('--Y', type=zero_or_more, default=0, help="Column Y in file 2, that should match with values in column X in file 1.")
	args = parser.parse_args()

	filter_file(args.file1, args.file2, args.s1, args.s2, args.X, args.Y)