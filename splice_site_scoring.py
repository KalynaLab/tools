#!/usr/bin/env python

"""
Splice site analysis based on position weight matrices (PWMs)
"""

import os
import argparse
from math import log

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

# Get splice site sequences
def get_fasta(bed_file, genome_fasta, output_file, upstream_nt=-14, downstream_nt=3):
	""" Extract the intronic fasta sequence, with going -x nt and +y nt in the neighboring exons."""

	import uuid
	import subprocess

	# Create a temporary bed file w/ adjustment for
	# the upstream_nt and downstream_nt parameters
	tmpBed = uuid.uuid4().hex + ".bed"
	while os.path.isfile(tmpBed):
		tmpBed = uuid.uuid4().hex + ".bed"

	with open(tmpBed, 'w') as fout:
		for line in open(bed_file):
			c, s, e, bed_id, score, strand = line.rstrip().split('\t')
			fout.write( "{}\t{}\t{}\t{}\t{}\t{}\n".format(c, int(s)-upstream_nt, int(e)+downstream_nt, bed_id, score, strand) )

	# Extract the fasta sequences
	subprocess.call("bedtools getfasta -s -name -fi {0} -bed {1} | sed 's/::.*//g' > {2}".format(genome_fasta, tmpBed, output_file), shell=True)

	# Clean-up temporary bed file
	os.remove(tmpBed)

def read_PWM_from_file(PWM_file, region=None):

	# Parse only a specific region of the PWM
	if region: 
		region = [ str(x+1) if x >= 0 else str(x) for x in [ i for i in xrange(*map(int, region.split(','))) ] ]

	PWM, parse, bases = [], True, []
	for line in open(PWM_file):
		if not line.startswith('#'):

			if line in ['\n', '\r\n']: parse = False
			elif parse:
				cols = line.rstrip().split('\t')
				try:
					total = sum([ float(x) for x in cols[1:] ])

					# Only parse position specified in region
					if region:
						if cols[0] in region:
							PWM.append( { bases[i-1]: float(cols[i])/total for i in xrange(1, len(cols)) } )	
					
					# No region specified, so parse everything
					else:
						PWM.append( { bases[i-1]: float(cols[i])/total for i in xrange(1, len(cols)) } )

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

def score_sites(bkgd_fa, test_fa, PWM5_file, PWM3_file, five_prime_region, three_prime_region, output_file):
	
	# Read PWM, take the splice site regions
	# into account
	PWM5 = { 'pwm': read_PWM_from_file(PWM5_file, five_prime_region), 'min': 0, 'max': 0 }
	PWM3 = { 'pwm': read_PWM_from_file(PWM3_file, three_prime_region), 'min': 0, 'max': 0 }

	five_size = sum([ abs(x) for x in map(int, five_prime_region.split(',')) ])
	three_size = sum([ abs(x) for x in map(int, three_prime_region.split(',')) ])

	# Score the background splice sites
	for record in yield_fasta(bkgd_fa):
		PWM5 = update_min_max(record.seq[:five_size].upper(), PWM5)
		PWM3 = update_min_max(record.seq[-three_size:].upper(), PWM3)

	# Score the test splice sites
	test, order = {}, []
	for record in yield_fasta(test_fa):
		
		five_seq = record.seq[:five_size].upper()
		three_seq = record.seq[-three_size:].upper()

		test[record.id] = { 
			'5seq': five_seq, 
			'3seq': three_seq,
			'5score': score_splice_site(five_seq, PWM5['pwm']),
			'3score': score_splice_site(three_seq, PWM3['pwm'])
		}
		PWM5 = update_min_max(five_seq, PWM5)
		PWM3 = update_min_max(three_seq, PWM3)

		order.append(record.id)

	# Output 
	with open(output_file, 'w') as fout:
		for x in order:
			fout.write( "{}\t{}\t{}\t{:.3f}\t{:.3f}\n".format(x, 
				test[x]['5seq'], 
				test[x]['3seq'], 
				rescale_score(test[x]['5score'], PWM5['min'], PWM5['max']),
				rescale_score(test[x]['3score'], PWM3['min'], PWM3['max']))
			)

if __name__ == '__main__':

	version = "0.0.1"
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-v', '--version', action='version', version=version, default=version)
	parser.add_argument('-w', '--work-dir', default="./", help="Output working directory.")

	subparsers = parser.add_subparsers(dest='command', help="Sub-command help.")

	# Extract fasta sequences
	parser_a = subparsers.add_parser('get-fasta', help="Extract the intronic fasta sequence, with going -x nt and +y nt in the neighboring exons (requires bedtools getfasta).")
	parser_a.add_argument('-b', '--bed', required=True, help="Bed annotation of intron coordinates.")
	parser_a.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta file.")
	parser_a.add_argument('-o', '--output-file', required=True, help="Output file.")
	parser_a.add_argument('--upstream-nt', type=int, default=-14, help="Number of nucleotides to read into the upstream (5'SS) exon.")
	parser_a.add_argument('--downstream-nt', type=int, default=3, help="Number of nucleotides to read into the downstream (3'SS) exon.")

	# Create PWM based on fasta file
	parser_b = subparsers.add_parser('make-PWM', help="Create a PWM from an input fasta file (in the format from the get-fasta command).")


	# Score splice sites
	parser_c = subparsers.add_parser('score-sites', help="Score splice sites based on supplied PWMs.")
	parser_c.add_argument('-b', '--bkgd', required=True, help="Background fasta sequences for score scaling (from get-fasta).")
	parser_c.add_argument('-t', '--test', required=True, help="Fasta sequences to score (from get-fasta).")
	parser_c.add_argument('--PWM5', required=True, help="5' splice site PWM file.")
	parser_c.add_argument('--PWM3', required=True, help="3' splice site PWM file.")
	parser_c.add_argument('-o', '--output-file', required=True, help="Output file.")
	parser_c.add_argument('--five-prime-region', default='-3,+10', help="Specify the number of exonic bases (-x) and intronic bases (+y) for the 5' splice site.")
	parser_c.add_argument('--three-prime-region', default='-14,+3', help="Specify the number of intronic base (-x) and exonic bases (+y) for the 3' splice site.")
	group = parser_c.add_mutually_exclusive_group()
	group.add_argument('--ath', action='store_true', help="Default preset for Arabidopsis thaliana: 5'ss => -3,+10, 3'ss => -14,+3.")
	group.add_argument('--hsa', action='store_true', help="Default preset for Homo sapiens: 5'ss => -3,+6, 3'ss => -6,+3.")

	args = parser.parse_args()

	work_dir = args.work_dir if args.work_dir.endswith('/') else args.work_dir+'/'
	if not os.path.exists(work_dir): os.makedirs(work_dir)

	if args.command == "get-fasta":
		get_fasta(args.bed, args.genome_fasta, args.output_file, args.upstream_nt, args.downstream_nt)

	elif args.command == "score-sites":

		# Set presets
		if args.ath:
			args.five_prime_region = "-3,+10"
			args.three_prime_region = "-14,+3"
		elif args.hsa:
			args.five_prime_region = "-3,+6"
			args.three_prime_region = "-6,+3"

		score_sites(args.bkgd, args.test, args.PWM5, args.PWM3, args.five_prime_region, args.three_prime_region, args.output_file)