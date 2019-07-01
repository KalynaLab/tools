#!/usr/bin/env python

""" Sequence info: Pipe a fasta file and get sequence information,
like GC content and sequence length. """

from __future__ import division
import argparse

def yield_fasta(f):
	''' Simple fasta parser that yield's the identifier and sequence of each record '''
	
	class SeqRecord:
		def __init__(self, seq_id, seq):
			self.id = seq_id
			self.seq = seq

	seq = ''
	for line in f:
		if line.startswith('>'):
			if len(seq):
				yield(SeqRecord(identifier, seq))
				seq = ''
			identifier = line[1:].rstrip()
		else: seq += line.rstrip()
	yield(SeqRecord(identifier, seq))

def seqInfo(fasta_file, get_length=False, get_dinu=False, get_GC=False):

	if any([ get_length, get_dinu, get_GC ]):
		
		for record in yield_fasta(fasta_file):
			output = [record.id]
			
			# Get sequence length
			if get_length: output.append(str(len(record.seq)))

			# Get dinucleotides (assumption is that the 
			# fasta sequence is the intron sequence)
			if get_dinu: output.append('{}-{}'.format(record.seq[:2], record.seq[-2:]))

			# Get GC content
			if get_GC: output.append('{:.3f}'.format(sum([ record.seq.count(char) for char in ['G', 'C'] ]) / len(record.seq) * 100))

			# Output
			print '\t'.join(output)

if __name__ == '__main__':

	import sys

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--fasta', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('--length', action='store_true', help="Output the sequence length.")	
	parser.add_argument('--dinu', action='store_true', help="Assumes the fasta sequence is an intron sequence. Output the dinucleotides.")
	parser.add_argument('--gc', action='store_true', help="Output the sequence GC content.")
	args = parser.parse_args()

	seqInfo(args.fasta, args.length, args.dinu, args.gc)