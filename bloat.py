#!/usr/bin/env python

""" bloat.py => Supplement file1 with information from file2, based on
column X in file1 and column Y in file2. The data in column Y in file2
should be unique. Line can be commented out by adding a hash-symbol (#) 
at the start of the line. The bloated output will be directed to stdout. """

import argparse

def bloat(file1, file2, s1="\t", s2="\t", X=0, Y=0, bloats='*'):

	# Get the bloating content from file2
	bloating = {}
	for line in open(file2):
		if not line.startswith('#'):
			cols = line.rstrip().split(s2)
			colY = cols[Y]

			# Only add the first occurence of bloated stuff
			# for the identifier at column Y
			if colY not in bloating:

				# Add all columns, except column Y
				if bloats == '*':
					del cols[Y]
					bloating[colY] = s2.join(cols)

				# Add specific columns only
				else:
					bloating[colY] = s2.join([ cols[int(x)] for x in bloats ])


	# Output the bloated file1 to stdout
	for line in open(file1):
		if not line.startswith('#'):
			cols = line.rstrip().split(s1)
			colX = cols[X]
			if cols[X] in bloating:
				print "{}{}{}".format(line.rstrip(), s1, bloating[colX])
			else: 
				print line.rstrip()

if __name__ == '__main__':

	def zero_or_more(val):
		if int(val) < 0:
			raise argparse.ArgumentTypeError("--column can't be lower than 0.")
		return int(val)

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--file1', required=True, help="File 1, the file to be bloated.")
	parser.add_argument('--file2', required=True, help="File 2, containing the bloat information.")
	parser.add_argument('--s1', default="\t", help="File separator for file1.")
	parser.add_argument('--s2', default="\t", help="File separator for file2.")
	parser.add_argument('--X', type=zero_or_more, default=0, help="ColumnX in file1 to match bloat information from file2.")
	parser.add_argument('--Y', type=zero_or_more, default=0, help="ColumnY in file2 that contains the identifier to match the bloat content to ColumnX in file1.")
	parser.add_argument('--bloats', default='*', nargs='+', help="Specify the column indices to bloat by, with the first column being 0.")
	args = parser.parse_args()

	bloat(args.file1, args.file2, args.s1, args.s2, args.X, args.Y, args.bloats)