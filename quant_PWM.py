#!/usr/bin/env python

""" Calculate the splice site strength for a set of sequences,
scaled based on the provided background sequences. 
"""

from __future__ import division
from natsort import natsorted
from math import log
import argparse

def yield_fasta(f):
	''' Simple fasta parser that yield's the identifier and sequence of each record '''
	
	class SeqRecord:
		def __init__(self, seq_id, seq):
			self.id = seq_id
			self.seq = seq

	seq = ''
	for line in open(f):
		if line.startswith('>'):
			if len(seq):
				yield(SeqRecord(identifier, seq))
				seq = ''
			identifier = line[1:].rstrip()
		else: seq += line.rstrip()
	yield(SeqRecord(identifier, seq))

def read_PWM_file(f):

	PWM, parse, bases = [], True, []
	for line in open(f):

		if not line.startswith('#'):

			if line in ['\n', '\r\n']: parse = False
			elif parse:

				cols = line.rstrip().split('\t')
				try: 
					total = sum([ float(x) for x in cols[1:] ])
					PWM.append( { bases[i-1]: float(cols[i])/total for i in xrange(1, len(cols)) })
				except ValueError: 
					bases = cols[1:]

	return PWM

def score_splice_site(seq, PWM):

	return sum([ log((PWM[i][base]/0.25)+0.0001, 2) for i, base in enumerate(seq) if base in PWM[i].keys() ])

def update_min_max(seq, PWM):

	score = score_splice_site(seq, PWM['pwm'])
	PWM['min'] = score if score < PWM['min'] else PWM['min']
	PWM['max'] = score if score > PWM['max'] else PWM['max']
	
	return PWM

def rescale_score(PWM_score, ss_min, ss_max):
    ''' Scale the PWM LOD score to 0-100 
    
        ((b-a)*(PWM_score-min)/(max-min))+a
        
        If the PWM score is positive -> a = 50, b = 100, min = 0
        If the PWM score is negative -> a = 0, b = 50, max = 0
    '''

    if PWM_score > 0: norm_score = ((50*PWM_score)/ss_max)+50
    else: norm_score = -(50*(PWM_score-ss_min)/ss_min)
        
    return norm_score

def quant_PWM(PWM_file, bkgd_fa, test_fa):

	# Parse the PWM
	PWM = { 'pwm': read_PWM_file(PWM_file), 'min': 0, 'max': 0 }
	
	# Quantify the background sequences
	for record in yield_fasta(bkgd_fa):
		PWM = update_min_max(record.seq.upper(), PWM)

	# Quantify the sequences of interest and output
	# the scaled score based on everything quantified
	test = {}
	for record in yield_fasta(test_fa):
		test[record.id] = { 'seq': record.seq.upper(), 'score': score_splice_site(record.seq.upper(), PWM['pwm']) }
		PWM = update_min_max(record.seq.upper(), PWM)

	print "#id\tsequence\tscore"
	for t in natsorted(test):
		print "{}\t{}\t{}".format(t, test[t]['seq'], test[t]['score'])

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--PWM', required=True, help="PWM file for either a 5' or 3' splice site.")
	parser.add_argument('--bkgd', required=True, help="Background 5' or 3' sequences (in fasta format).")
	parser.add_argument('--test', required=True, help="5' or 3' sequences of interest (in fasta format).")
	args = parser.parse_args()

	quant_PWM(args.PWM, args.bkgd, args.test)