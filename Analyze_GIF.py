"""Import Modules"""
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
import os, sys
from collections import Counter
import numpy as np
import operator
import string
import itertools
import fuzzysearch
import xlsxwriter
from PIL import Image
import png

"""Line Arguments"""
run_number = sys.argv[1] #MiSeq run number for Data_Path
condition = sys.argv[2] #file name (minus .fastq) shoud be the only argument following SPCR_blast.py when running script
re_blast = sys.argv[3] #either yes or no, (to re-do the blast)

"""Globals"""
user_profile = os.environ ['USERPROFILE']
Data_Path = '%s/Data Analysis/MS%s/%s_Results' % (user_profile,run_number,condition) #this should be a folder with the fastq in it
Blast_Data_Path = '%s/Blast_Databases' % user_profile  #this should be a folder with NCBI blast databases
Proto_Data_Path = '%s/Image_SPCR_Databases' % user_profile #this should be a folder with fasta files containing the intended SPCRs for checking errors

bin_Seq = {'C': '00', 'T': '01', 'A': '10', 'G': '11'}
Trip_bin = {'0000': 'TGA', '1000': 'AGA', '0100': 'TAA', '1100': 'CAC', '0010': 'GAC', '1010': 'AGC', '0110': 'GGA', '1110': 'TGG'}
bin_Trip = {y:x for x,y in Trip_bin.iteritems()}
#This is the nucleotide triplet to pixel value conversion, numbers depend on the palette used, in this case extracted from the image in ImageJ
trip_dict = {20: ['TTT','CCC','AAA'],
			 224: ['TCT','CAC','AGA'],
			 179: ['TAT','CGC','ATG'],
			 240: ['TGT','CTA','ACG'],
			 100: ['TTC','CCA','AGG'],
			 232: ['TCC','CAA','GTT'],
			 216: ['TAC','CGA','GCT'],
			 255: ['TGC','CTG','GAT'],
			 59: ['TTA','CCG','GGT'],
			 208: ['TCA','CAG','GTC'],
			 148:['TAA','CGG','GCC'],
			 43:['TGA','ATT','GAC'],
			 200:['TTG','ACT','GGC'],
			 127:['TCG','AAT','GTA'],
			 83:['TAG','AGT','GCA'],
			 168:['TGG','ATC','GAA'],
			 32:['CTT','ACC','GGA'],
			 192:['CCT','AAC','GTG'],
			 112:['CAT','AGC','GCG'],
			 72:['CGT','ATA','GAG'],
			 160:['CTC','ACA','GGG'], 
			 0:['AAG']} 

trip_to_num = dict( (v,k) for k in trip_dict for v in trip_dict[k] )
PixEt_limited_codon_set = ['TGA','AGA','TAA','CAC','GAC','AGC','GGA','TGG']

unaligned_IDs = []
unaligned_SPCRs = []
unaligned_SPCR_nuc = []
Pixet_dict = {}
Pixet_avg_count_dict = {}
Pixet_least_count_dict = {}
Too_few_Pixets_list = []
Pixets_by_count_list = []
Pixets_by_least_count_list = []
Score_dict = {}
double_list = []
triple_list = []
Pixets_5count_noPerfect_test1 = []
Pixets_5count_Perfect_test1 = []
Set_Frame1_SPCRs = []
Set_Frame2_SPCRs = []
Set_Frame3_SPCRs = []
Set_Frame4_SPCRs = []
Set_Frame5_SPCRs = []

dist_repeat = 4
dist_SPCRs = 5
Repeat = 'GTGTTCCCCGCGCCAGCGGGGATAAACC'

Frame_SPCR_dict = {	'Frame1': [],
					'Frame2': [],
					'Frame3': [],
					'Frame4': [],
					'Frame5': []}

"""Defs"""
def get_context(record):
	"""input is a blastn record
		returns info and protospacer context"""
	direction = []
	if not blast_record.alignments:
		unaligned_IDs.append (blast_record.query)

def seq_to_bin_num(seq):
	seq_code = [bin_Seq[seq[0]]]
	seq_code.append(bin_Seq[seq[1]])
	seq_code.append(bin_Trip[seq[2:5]])
	PixEt_bin = ''.join(seq_code)
	PixEt_bin_rev = str(PixEt_bin)[::-1]
	PixEt_int = int(PixEt_bin_rev, 2)
	return PixEt_int

# @profile
def spcr_order_doubles1(doubles, Fed_dict):
    First = 'N'
    Second = 'N'
    if Fed_dict['A'] == doubles[0]:
    	First = 'A'
    elif Fed_dict['B'] == doubles[0]:
    	First = 'B'
    elif Fed_dict['C'] == doubles[0]:
    	First = 'C'
    elif Fed_dict['D'] == doubles[0]:
    	First = 'D'
    elif Fed_dict['E'] == doubles[0]:
    	First = 'E'
    
    if Fed_dict['A'] == doubles[1]:
    	Second = 'A'
    elif Fed_dict['B'] == doubles[1]:
    	Second = 'B'
    elif Fed_dict['C'] == doubles[1]:
    	Second = 'C'
    elif Fed_dict['D'] == doubles[1]:
    	Second = 'D'
    elif Fed_dict['E'] == doubles[1]:
    	Second = 'E'

    order = First+Second
    return order

# @profile
def spcr_order_doubles_allbut1(doubles, Fed_dict):
	First = 'N'
	Second = 'N'
	if Fed_dict['A'] == doubles[0]:
		First = 'A'
	elif Fed_dict['B'] == doubles[0]:
		First = 'B'
	elif Fed_dict['C'] == doubles[0]:
		First = 'C'
	elif Fed_dict['D'] == doubles[0]:
		First = 'D'
	elif Fed_dict['E'] == doubles[0]:
		First = 'E'
	elif doubles[0] in Set_Frame1_SPCRs:
		First = '1'
	elif doubles[0] in Set_Frame2_SPCRs:
		First = '2'
	elif doubles[0] in Set_Frame3_SPCRs:
		First = '3'
	elif doubles[0] in Set_Frame4_SPCRs:
		First = '4'
	elif doubles[0] in Set_Frame5_SPCRs:
		First = '5'

	if Fed_dict['A'] == doubles[1]:
		Second = 'A'
	elif Fed_dict['B'] == doubles[1]:
		Second = 'B'
	elif Fed_dict['C'] == doubles[1]:
		Second = 'C'
	elif Fed_dict['D'] == doubles[1]:
		Second = 'D'
	elif Fed_dict['E'] == doubles[1]:
		Second = 'E'
	elif doubles[1] in Set_Frame1_SPCRs:
		Second = '1'
	elif doubles[1] in Set_Frame2_SPCRs:
		Second = '2'
	elif doubles[1] in Set_Frame3_SPCRs:
		Second = '3'
	elif doubles[1] in Set_Frame4_SPCRs:
		Second = '4'
	elif doubles[1] in Set_Frame5_SPCRs:
		Second = '5'

	order = First+Second
	return order

# @profile
def spcr_order_triples1(triples, Fed_dict):
    First = 'N'
    Second = 'N'
    Third = 'N'
    if Fed_dict['A'] == triples[0]:
    	First = 'A'
    elif Fed_dict['B'] == triples[0]:
    	First = 'B'
    elif Fed_dict['C'] == triples[0]:
    	First = 'C'
    elif Fed_dict['D'] == triples[0]:
    	First = 'D'
    elif Fed_dict['E'] == triples[0]:
    	First = 'E'
    
    if Fed_dict['A'] == triples[1]:
    	Second = 'A'
    elif Fed_dict['B'] == triples[1]:
    	Second = 'B'
    elif Fed_dict['C'] == triples[1]:
    	Second = 'C'
    elif Fed_dict['D'] == triples[1]:
    	Second = 'D'
    elif Fed_dict['E'] == triples[1]:
    	Second = 'E'
    
    if Fed_dict['A'] == triples[2]:
    	Third = 'A'
    elif Fed_dict['B'] == triples[2]:
    	Third = 'B'
    elif Fed_dict['C'] == triples[2]:
    	Third = 'C'
    elif Fed_dict['D'] == triples[2]:
    	Third = 'D'
    elif Fed_dict['E'] == triples[2]:
    	Third = 'E'

    order = [First+Second]
    order.append(Second+Third)
    return order

# @profile
def spcr_order_triples_allbut1(triples, Fed_dict):
	First = 'N'
	Second = 'N'
	Third = 'N'
	if Fed_dict['A'] == triples[0]:
		First = 'A'
	elif Fed_dict['B'] == triples[0]:
		First = 'B'
	elif Fed_dict['C'] == triples[0]:
		First = 'C'
	elif Fed_dict['D'] == triples[0]:
		First = 'D'
	elif Fed_dict['E'] == triples[0]:
		First = 'E'
	elif triples[0] in Set_Frame1_SPCRs:
		First = '1'
	elif triples[0] in Set_Frame2_SPCRs:
		First = '2'
	elif triples[0] in Set_Frame3_SPCRs:
		First = '3'
	elif triples[0] in Set_Frame3_SPCRs:
		First = '4'
	elif triples[0] in Set_Frame5_SPCRs:
		First = '5'
    
	if Fed_dict['A'] == triples[1]:
		Second = 'A'
	elif Fed_dict['B'] == triples[1]:
		Second = 'B'
	elif Fed_dict['C'] == triples[1]:
		Second = 'C'
	elif Fed_dict['D'] == triples[1]:
		Second = 'D'
	elif Fed_dict['E'] == triples[1]:
		Second = 'E'
	elif triples[1] in Set_Frame1_SPCRs:
		Second = '1'
	elif triples[1] in Set_Frame2_SPCRs:
		Second = '2'
	elif triples[1] in Set_Frame3_SPCRs:
		Second = '3'
	elif triples[1] in Set_Frame4_SPCRs:
		Second = '4'
	elif triples[1] in Set_Frame5_SPCRs:
		Second = '5'
    
	if Fed_dict['A'] == triples[2]:
		Third = 'A'
	elif Fed_dict['B'] == triples[2]:
		Third = 'B'
	elif Fed_dict['C'] == triples[2]:
		Third = 'C'
	elif Fed_dict['D'] == triples[2]:
		Third = 'D'
	elif Fed_dict['E'] == triples[2]:
		Third = 'E'
	elif triples[2] in Set_Frame1_SPCRs:
		Third = '1'
	elif triples[2] in Set_Frame2_SPCRs:
		Third = '2'
	elif triples[2] in Set_Frame3_SPCRs:
		Third = '3'
	elif triples[2] in Set_Frame4_SPCRs:
		Third = '4'
	elif triples[2] in Set_Frame5_SPCRs:
		Third = '5'

	order = [First+Second]
	order.append(Second+Third)
	return order

# @profile
def test_order(permutation): #tests against other spacers in pixet
	test_dict = {}
	ratio_list = []
	if order_dict['%s%s' % (permutation[1], permutation[0])] > order_dict['%s%s' % (permutation[0], permutation[1])]:
	    test_dict['Test_1'] = 1
	else:
	    test_dict['Test_1'] = 0
	if order_dict['%s%s' % (permutation[2], permutation[1])] > order_dict['%s%s' % (permutation[1], permutation[2])]:
	    test_dict['Test_2'] = 1
	else:
	    test_dict['Test_2'] = 0
	if order_dict['%s%s' % (permutation[3], permutation[2])] > order_dict['%s%s' % (permutation[2], permutation[3])]:
	    test_dict['Test_3'] = 1
	else:
	    test_dict['Test_3'] = 0
	if order_dict['%s%s' % (permutation[4], permutation[3])] > order_dict['%s%s' % (permutation[3], permutation[4])]:
	    test_dict['Test_4'] = 1
	else:
	    test_dict['Test_4'] = 0
	if order_dict['%s%s' % (permutation[2], permutation[0])] > order_dict['%s%s' % (permutation[0], permutation[2])]:
	    test_dict['Test_5'] = 1
	else:
	    test_dict['Test_5'] = 0
	if order_dict['%s%s' % (permutation[3], permutation[1])] > order_dict['%s%s' % (permutation[1], permutation[3])]:
	    test_dict['Test_6'] = 1
	else:
	    test_dict['Test_6'] = 0
	if order_dict['%s%s' % (permutation[4], permutation[2])] > order_dict['%s%s' % (permutation[2], permutation[4])]:
	    test_dict['Test_7'] = 1
	else:
	    test_dict['Test_7'] = 0
	if order_dict['%s%s' % (permutation[3], permutation[0])] > order_dict['%s%s' % (permutation[0], permutation[3])]:
	    test_dict['Test_8'] = 1
	else:
	    test_dict['Test_8'] = 0
	if order_dict['%s%s' % (permutation[4], permutation[1])] > order_dict['%s%s' % (permutation[1], permutation[4])]:
	    test_dict['Test_9'] = 1
	else:
	    test_dict['Test_9'] = 0
	if order_dict['%s%s' % (permutation[4], permutation[0])] > order_dict['%s%s' % (permutation[0], permutation[4])]:
	    test_dict['Test_10'] = 1
	else:
	    test_dict['Test_10'] = 0
	if order_dict['N%s' % (permutation[0])] > 0 and order_dict['%sN' % (permutation[0])] > 0:
		ratio_list.append(float(order_dict['N%s' % (permutation[0])])/order_dict['%sN' % (permutation[0])])
	if order_dict['N%s' % (permutation[1])] > 0 and order_dict['%sN' % (permutation[1])] > 0:
		ratio_list.append(float(order_dict['N%s' % (permutation[1])])/order_dict['%sN' % (permutation[1])])
	if order_dict['N%s' % (permutation[2])] > 0 and order_dict['%sN' % (permutation[2])] > 0:
		ratio_list.append(float(order_dict['N%s' % (permutation[2])])/order_dict['%sN' % (permutation[2])])
	if order_dict['N%s' % (permutation[3])] > 0 and order_dict['%sN' % (permutation[3])] > 0:
		ratio_list.append(float(order_dict['N%s' % (permutation[3])])/order_dict['%sN' % (permutation[3])])
	if order_dict['N%s' % (permutation[4])] > 0 and order_dict['%sN' % (permutation[4])] > 0:
		ratio_list.append(float(order_dict['N%s' % (permutation[4])])/order_dict['%sN' % (permutation[4])])
	if ratio_list == sorted(ratio_list, reverse=True): #checks for a descending ratio of NtoSpcrs
		test_dict['Test_11'] = len(ratio_list)
	else:
		test_dict['Test_11'] = 0
	return test_dict

# @profile
def test_order2(permutation): #tests against spacers that have already been asigned to a frame
	test_dict = {}
	if 0.7*(order_dict['1%s' % permutation[0]]+order_dict['%s1' % permutation[0]]) > order_dict['1%s' % permutation[0]] > 0.3*(order_dict['1%s' % permutation[0]]+order_dict['%s1' % permutation[0]]):
	    test_dict['Test_1'] = 1
	else:
	    test_dict['Test_1'] = 0
	if order_dict['2%s' % permutation[0]] > order_dict['%s2' % permutation[0]]:
	    test_dict['Test_2'] = 1
	else:
	    test_dict['Test_2'] = 0
	if order_dict['3%s' % permutation[0]] > order_dict['%s3' % permutation[0]]:
	    test_dict['Test_3'] = 1
	else:
	    test_dict['Test_3'] = 0
	if order_dict['4%s' % permutation[0]] > order_dict['%s4' % permutation[0]]:
	    test_dict['Test_4'] = 1
	else:
	    test_dict['Test_4'] = 0
	if order_dict['5%s' % permutation[0]] > order_dict['%s5' % permutation[0]]:
	    test_dict['Test_5'] = 1
	else:
	    test_dict['Test_5'] = 0

	if order_dict['1%s' % permutation[1]] < order_dict['%s1' % permutation[1]]:
	    test_dict['Test_6'] = 1
	else:
	    test_dict['Test_6'] = 0
	if 0.7*(order_dict['2%s' % permutation[1]]+order_dict['%s2' % permutation[1]]) > order_dict['2%s' % permutation[1]] > 0.3*(order_dict['2%s' % permutation[1]]+order_dict['%s2' % permutation[1]]):
	    test_dict['Test_7'] = 1
	else:
	    test_dict['Test_7'] = 0
	if order_dict['3%s' % permutation[1]] > order_dict['%s3' % permutation[1]]:
	    test_dict['Test_8'] = 1
	else:
	    test_dict['Test_8'] = 0
	if order_dict['4%s' % permutation[1]] > order_dict['%s4' % permutation[1]]:
	    test_dict['Test_9'] = 1
	else:
	    test_dict['Test_9'] = 0
	if order_dict['5%s' % permutation[1]] > order_dict['%s5' % permutation[1]]:
	    test_dict['Test_10'] = 1
	else:
	    test_dict['Test_10'] = 0

	if order_dict['1%s' % permutation[2]] < order_dict['%s1' % permutation[2]]:
	    test_dict['Test_11'] = 1
	else:
	    test_dict['Test_11'] = 0
	if order_dict['2%s' % permutation[2]] < order_dict['%s2' % permutation[2]]:
	    test_dict['Test_12'] = 1
	else:
	    test_dict['Test_12'] = 0
	if 0.7*(order_dict['3%s' % permutation[2]]+order_dict['%s3' % permutation[2]]) > order_dict['3%s' % permutation[2]] > 0.3*(order_dict['3%s' % permutation[2]]+order_dict['%s3' % permutation[2]]):
	    test_dict['Test_13'] = 1
	else:
	    test_dict['Test_13'] = 0
	if order_dict['4%s' % permutation[2]] > order_dict['%s4' % permutation[2]]:
	    test_dict['Test_14'] = 1
	else:
	    test_dict['Test_14'] = 0
	if order_dict['5%s' % permutation[2]] > order_dict['%s5' % permutation[2]]:
	    test_dict['Test_15'] = 1
	else:
	    test_dict['Test_15'] = 0

	if order_dict['1%s' % permutation[3]] < order_dict['%s1' % permutation[3]]:
	    test_dict['Test_16'] = 1
	else:
	    test_dict['Test_16'] = 0
	if order_dict['2%s' % permutation[3]] < order_dict['%s2' % permutation[3]]:
	    test_dict['Test_17'] = 1
	else:
	    test_dict['Test_17'] = 0
	if order_dict['3%s' % permutation[3]] < order_dict['%s3' % permutation[3]]:
	    test_dict['Test_18'] = 1
	else:
	    test_dict['Test_18'] = 0
	if 0.7*(order_dict['4%s' % permutation[3]]+order_dict['%s4' % permutation[3]]) > order_dict['4%s' % permutation[3]] > 0.3*(order_dict['4%s' % permutation[3]]+order_dict['%s4' % permutation[3]]):
	    test_dict['Test_19'] = 1
	else:
	    test_dict['Test_19'] = 0
	if order_dict['5%s' % permutation[3]] > order_dict['%s5' % permutation[3]]:
	    test_dict['Test_20'] = 1
	else:
	    test_dict['Test_20'] = 0

	if order_dict['1%s' % permutation[4]] < order_dict['%s1' % permutation[4]]:
	    test_dict['Test_21'] = 1
	else:
	    test_dict['Test_21'] = 0
	if order_dict['2%s' % permutation[4]] < order_dict['%s2' % permutation[4]]:
	    test_dict['Test_22'] = 1
	else:
	    test_dict['Test_22'] = 0
	if order_dict['3%s' % permutation[4]] < order_dict['%s3' % permutation[4]]:
	    test_dict['Test_23'] = 1
	else:
	    test_dict['Test_23'] = 0
	if order_dict['4%s' % permutation[4]] < order_dict['%s4' % permutation[4]]:
	    test_dict['Test_24'] = 1
	else:
	    test_dict['Test_24'] = 0
	if 0.7*(order_dict['5%s' % permutation[4]]+order_dict['%s5' % permutation[4]]) > order_dict['5%s' % permutation[4]] > 0.3*(order_dict['5%s' % permutation[4]]+order_dict['%s5' % permutation[4]]):
	    test_dict['Test_25'] = 1
	else:
	    test_dict['Test_25'] = 0

	return test_dict

"""Run"""
#First, blast (if not already done)
if re_blast == 'yes':
	blastn_cline = NcbiblastnCommandline(task='blastn-short', query="%s/new_SPCRs_seqs.fasta" % Data_Path, db="%s/Specific_Blast_Database" % Blast_Data_Path, out="%s/blastn_output.xml" % Data_Path, evalue=0.0001, outfmt=5)
	stdout, stdeff = blastn_cline()
#open blastn file
SPCR_blast = open("%s/blastn_output.xml" % Data_Path)
from Bio.Blast import NCBIXML
blast_records = NCBIXML.parse(SPCR_blast) #iterator
#analyze spcrs
for blast_record in blast_records:
	spcr_context = get_context(blast_record)
record_dict = SeqIO.index("%s/new_SPCRs_seqs.fasta" % Data_Path, "fasta")
# record_dict = SeqIO.to_dict(SeqIO.parse("%s/new_SPCRs_seqs.fasta" % Data_Path, "fasta"))
for ID in unaligned_IDs:
	unaligned_SPCRs.append (record_dict[ID.split(' ')[0]])
	unaligned_SPCR_nuc.append (record_dict[ID.split(' ')[0]].seq)
c = Counter(unaligned_SPCR_nuc)
freq_list = c.most_common()
# remove spacers that have inserted in reversed orientation (if spacer is found in other orientations, take only the one with the higher count) 
freq_list_revcomp = []
freq_list_less_revcomps = []
for seq in freq_list:
	if not any(seq[0][2:31] in s for s in freq_list_revcomp):
		freq_list_less_revcomps.append(seq)
	freq_list_revcomp.append(seq[0].reverse_complement())
#make a dictionary with all potential SPCRs for each pixet along with their freq_list
for seq in freq_list_less_revcomps:
	if len(seq[0]) == 33 and seq[0][0] == 'G':
		PixEt_ID = seq[0][28:33]
		if PixEt_ID[2:5] in PixEt_limited_codon_set:
			PixEt_ID_num = seq_to_bin_num(PixEt_ID)
			if PixEt_ID_num not in [8,40,120]:
				if PixEt_ID_num == 126:		#this fixes disallowed Pixets because of CTT
					PixEt = 8
				elif PixEt_ID_num == 127:
					PixEt = 40
				elif PixEt_ID_num == 128:
					PixEt = 120
				else:
					PixEt = PixEt_ID_num
				if PixEt in Pixet_dict and len(Pixet_dict[PixEt]) < 5:
					Pixet_dict[PixEt].append(seq)
				elif PixEt not in Pixet_dict:
					Pixet_dict[PixEt] = [seq]
#rank all Pixets by the average count aross frames (considering only those that have 5 or more "Frames")
for PixEt in range(1,105):
	counts = []
	if PixEt in Pixet_dict and len(Pixet_dict[PixEt]) == 5:
		for i in range(0,5):
			counts.append(Pixet_dict[PixEt][i][1])
		Pixet_least_count_dict[PixEt] = counts[4]
	elif PixEt in Pixet_dict and len(Pixet_dict[PixEt]) < 5:
		Too_few_Pixets_list.append(PixEt)

print "Pixets with 5 or more sequences: %s" % len(Pixet_least_count_dict)
print "Pixets with fewer than 5 sequences: %s" % len(Too_few_Pixets_list)

Pixets_ordered_by_least_count = sorted(Pixet_least_count_dict.items(), key=operator.itemgetter(1))
Pixets_ordered_by_least_count.reverse() #makes into a descending tuple
for P in Pixets_ordered_by_least_count:
	Pixets_by_least_count_list.append(P[0])

for doubles in SeqIO.parse("%s/double_expansion_sequences_two_read_seqs.fastq" % Data_Path, "fastq"):
    repeats = fuzzysearch.find_near_matches(Repeat, doubles.seq, max_l_dist=dist_repeat)
    spacer_list = [doubles.seq[repeats[0].end:repeats[1].start]]
    spacer_list.append (doubles.seq[repeats[1].end:repeats[2].start])
    double_list.append(spacer_list)
for doubles in SeqIO.parse("%s/double_expansion_sequences_three_read_seqs.fastq" % Data_Path, "fastq"):
    repeats = fuzzysearch.find_near_matches(Repeat, doubles.seq, max_l_dist=dist_repeat)
    spacer_list = [doubles.seq[repeats[0].end:repeats[1].start]]
    spacer_list.append (doubles.seq[repeats[1].end:repeats[2].start])
    double_list.append(spacer_list)
for triples in SeqIO.parse("%s/triple_expansion_sequences_three_read_seqs.fastq" % Data_Path, "fastq"):
    repeats = fuzzysearch.find_near_matches(Repeat, triples.seq, max_l_dist=dist_repeat)
    spacer_list = [triples.seq[repeats[0].end:repeats[1].start]]
    spacer_list.append (triples.seq[repeats[1].end:repeats[2].start])
    spacer_list.append (triples.seq[repeats[2].end:repeats[3].start])
    triple_list.append(spacer_list)

#get the first score for all pixets with 5 or more pixets, put spcrs into frames if perfect, else save pixet
for ranked_pixet in Pixets_by_least_count_list:
	Current_set_dict = {}
	order_dict = {}
	internal_score_dict = {}
	internal_score_dict_summed = {}
	# order_list = []
	for tup in itertools.product('ABCDEN12345',repeat=2):
		string_tup = ''.join(str(i) for i in tup)
		order_dict[string_tup] = 0
		# order_list.append(string)

	for i in range(0,5):
		Current_set_dict[string.ascii_uppercase[i]] = Pixet_dict[ranked_pixet][i][0]

	for doubles in double_list:
	    order = spcr_order_doubles1(doubles,Current_set_dict)
	    order_dict[order] +=1
	for triples in triple_list:
	    order = spcr_order_triples1(triples,Current_set_dict)
	    order_dict[order[0]] +=1
	    order_dict[order[1]] +=1

	for permutation in itertools.permutations(string.ascii_uppercase[:len(Current_set_dict)]):
		test_dict_current = test_order(permutation)
		score = sum(test_dict_current.values())
		internal_score_dict[permutation] = score

	if max(internal_score_dict.iteritems(), key=operator.itemgetter(1))[0] == 15:
		predicted_order = max(internal_score_dict.iteritems(), key=operator.itemgetter(1))[0]
		Frame_SPCR_dict['Frame1'].append(Current_set_dict[predicted_order[0]])
		Frame_SPCR_dict['Frame2'].append(Current_set_dict[predicted_order[1]])
		Frame_SPCR_dict['Frame3'].append(Current_set_dict[predicted_order[2]])
		Frame_SPCR_dict['Frame4'].append(Current_set_dict[predicted_order[3]])
		Frame_SPCR_dict['Frame5'].append(Current_set_dict[predicted_order[4]])
		Score_dict[ranked_pixet] = [internal_score_dict[predicted_order],'n/a']
		Pixets_5count_Perfect_test1.append(ranked_pixet)
	else:
		Pixets_5count_noPerfect_test1.append(ranked_pixet)

print 'Pixets with 5 or more and a perfect score on test1 (analyzed): %s' % len(Pixets_5count_Perfect_test1)
print 'Pixets with 5 or more without a perfect test1: %s' % len(Pixets_5count_noPerfect_test1)

#get order ranked pixets without a perfect test1
for ranked_pixet in Pixets_5count_noPerfect_test1:
	Current_set_dict = {}
	order_dict = {}
	internal_score_dict = {}
	internal_score_dict_rev = {}
	internal_score_dict_summed = {}
	Set_Frame1_SPCRs = set(Frame_SPCR_dict['Frame1'])
	Set_Frame2_SPCRs = set(Frame_SPCR_dict['Frame2'])
	Set_Frame3_SPCRs = set(Frame_SPCR_dict['Frame3'])
	Set_Frame4_SPCRs = set(Frame_SPCR_dict['Frame4'])
	Set_Frame5_SPCRs = set(Frame_SPCR_dict['Frame5'])
	# order_list = []
	for tup in itertools.product('ABCDEN12345',repeat=2):
		string_tup = ''.join(str(i) for i in tup)
		order_dict[string_tup] = 0
		# order_list.append(string)

	for i in range(0,5):
		Current_set_dict[string.ascii_uppercase[i]] = Pixet_dict[ranked_pixet][i][0]

	for doubles in double_list:
	    order = spcr_order_doubles_allbut1(doubles,Current_set_dict)
	    order_dict[order] +=1
	for triples in triple_list:
	    order = spcr_order_triples_allbut1(triples,Current_set_dict)
	    order_dict[order[0]] +=1
	    order_dict[order[1]] +=1

	for permutation in itertools.permutations(string.ascii_uppercase[:len(Current_set_dict)]):
		test_dict1 = test_order(permutation)
		test_dict2 = test_order2(permutation)
		score1 = sum(test_dict1.values())
		score2 = sum(test_dict2.values())
		internal_score_dict[permutation] = [score1,score2] 
		internal_score_dict_rev[permutation] = [score2,score1] #this orderr will give the higest second score with the highest first score

	predicted_order = max(internal_score_dict_rev.iteritems(), key=operator.itemgetter(1))[0]
	Frame_SPCR_dict['Frame1'].append(Current_set_dict[predicted_order[0]])
	Frame_SPCR_dict['Frame2'].append(Current_set_dict[predicted_order[1]])
	Frame_SPCR_dict['Frame3'].append(Current_set_dict[predicted_order[2]])
	Frame_SPCR_dict['Frame4'].append(Current_set_dict[predicted_order[3]])
	Frame_SPCR_dict['Frame5'].append(Current_set_dict[predicted_order[4]])
	Score_dict[ranked_pixet] = internal_score_dict[predicted_order]

"""next, need to deal with pixets that have fewer than 5 spcrs"""
#maybe this won't happen, fix by adding 'dummy' sequences
for pixet in Too_few_Pixets_list:
	Current_set_dict = {}
	order_dict = {}
	internal_score_dict = {}
	internal_score_dict_summed = {}
	Set_Frame1_SPCRs = set(Frame_SPCR_dict['Frame1'])
	Set_Frame2_SPCRs = set(Frame_SPCR_dict['Frame2'])
	Set_Frame3_SPCRs = set(Frame_SPCR_dict['Frame3'])
	Set_Frame4_SPCRs = set(Frame_SPCR_dict['Frame4'])
	Set_Frame5_SPCRs = set(Frame_SPCR_dict['Frame5'])
	for tup in itertools.product('ABCDEN12345',repeat=2):
		string_tup = ''.join(str(i) for i in tup)
		order_dict[string_tup] = 0

	for i in range(0,len(Pixet_dict[pixet])):
		Current_set_dict[string.ascii_uppercase[i]] = Pixet_dict[pixet][i][0]

	for i in range(len(Pixet_dict[pixet]),5):
		Current_set_dict[string.ascii_uppercase[i]] = 'dummy'

	for doubles in double_list:
	    order = spcr_order_doubles_allbut1(doubles,Current_set_dict)
	    order_dict[order] +=1
	for triples in triple_list:
	    order = spcr_order_triples_allbut1(triples,Current_set_dict)
	    order_dict[order[0]] +=1
	    order_dict[order[1]] +=1

	#test these based on only test set 2    
	for permutation in itertools.permutations(string.ascii_uppercase[:len(Current_set_dict)]):
		# test_dict1 = test_order(permutation)
		test_dict2 = test_order2(permutation)
		# score1 = sum(test_dict1.values())
		score2 = sum(test_dict2.values())
		internal_score_dict[permutation] = score2

	predicted_order = max(internal_score_dict.iteritems(), key=operator.itemgetter(1))[0]
	if Current_set_dict[predicted_order[0]] is not 'dummy':
		Frame_SPCR_dict['Frame1'].append(Current_set_dict[predicted_order[0]])
	if Current_set_dict[predicted_order[1]] is not 'dummy':
		Frame_SPCR_dict['Frame2'].append(Current_set_dict[predicted_order[1]])
	if Current_set_dict[predicted_order[2]] is not 'dummy':
		Frame_SPCR_dict['Frame3'].append(Current_set_dict[predicted_order[2]])
	if Current_set_dict[predicted_order[3]] is not 'dummy':
		Frame_SPCR_dict['Frame4'].append(Current_set_dict[predicted_order[3]])
	if Current_set_dict[predicted_order[4]] is not 'dummy':
		Frame_SPCR_dict['Frame5'].append(Current_set_dict[predicted_order[4]])
	Score_dict[pixet] = ['n/a',internal_score_dict[predicted_order]]

print "Next %s pixets analyzed" % len(Too_few_Pixets_list)
Pixets_analyzed = Pixets_5count_Perfect_test1+Pixets_5count_noPerfect_test1+Too_few_Pixets_list

#compare against correct pixets
proto_frame1 = SeqIO.index("%s/GIF_protos_frame1.fasta" % Proto_Data_Path, "fasta")
proto_frame2 = SeqIO.index("%s/GIF_protos_frame2.fasta" % Proto_Data_Path, "fasta")
proto_frame3 = SeqIO.index("%s/GIF_protos_frame3.fasta" % Proto_Data_Path, "fasta")
proto_frame4 = SeqIO.index("%s/GIF_protos_frame4.fasta" % Proto_Data_Path, "fasta")
proto_frame5 = SeqIO.index("%s/GIF_protos_frame5.fasta" % Proto_Data_Path, "fasta")
Frame1_score = 0
Frame2_score = 0
Frame3_score = 0
Frame4_score = 0
Frame5_score = 0
for key in proto_frame1:
	if str(proto_frame1[key].seq[2:35]) in Frame_SPCR_dict['Frame1']:
		Frame1_score += 1
Frame1_score /= float(len(proto_frame1))
Frame1_score *= 100
for key in proto_frame2:
	if str(proto_frame2[key].seq[2:35]) in Frame_SPCR_dict['Frame2']:
		Frame2_score += 1
Frame2_score /= float(len(proto_frame2))
Frame2_score *= 100
for key in proto_frame3:
	if str(proto_frame3[key].seq[2:35]) in Frame_SPCR_dict['Frame3']:
		Frame3_score += 1
Frame3_score /= float(len(proto_frame3))
Frame3_score *= 100
for key in proto_frame4:
	if str(proto_frame4[key].seq[2:35]) in Frame_SPCR_dict['Frame4']:
		Frame4_score += 1
Frame4_score /= float(len(proto_frame4))
Frame4_score *= 100
for key in proto_frame5:
	if str(proto_frame5[key].seq[2:35]) in Frame_SPCR_dict['Frame5']:
		Frame5_score += 1
Frame5_score /= float(len(proto_frame5))
Frame5_score *= 100

"""Output"""
#Score_dict
workbook = xlsxwriter.Workbook('%s/Scores.xlsx' % (Data_Path))
worksheet = workbook.add_worksheet('Scores_by_Pixet')
bold = workbook.add_format({'bold': True})
#headings
worksheet.write(0,0,'Pixet', bold)
worksheet.write(0,1,'Test Set 1', bold)
worksheet.write(0,2,'Test Set 2', bold)
worksheet.write(0,3,'Frame1', bold)
worksheet.write(0,4,'Frame2', bold)
worksheet.write(0,5,'Frame3', bold)
worksheet.write(0,6,'Frame4', bold)
worksheet.write(0,7,'Frame5', bold)
#Data
row = 1
col = 0
for pixet in Pixets_analyzed:
	worksheet.write(row,col,pixet)
	worksheet.write(row,col+1,Score_dict[pixet][0])
	worksheet.write(row,col+2,Score_dict[pixet][1])
	if proto_frame1['%s' % pixet][2:35].seq in set(Frame_SPCR_dict['Frame1']):
		worksheet.write(row,col+3,1)
	else:
		worksheet.write(row,col+3,0)
	if proto_frame2['%s' % pixet][2:35].seq in set(Frame_SPCR_dict['Frame2']):
		worksheet.write(row,col+4,1)
	else:
		worksheet.write(row,col+4,0)
	if proto_frame3['%s' % pixet][2:35].seq in set(Frame_SPCR_dict['Frame3']):
		worksheet.write(row,col+5,1)
	else:
		worksheet.write(row,col+5,0)
	if proto_frame4['%s' % pixet][2:35].seq in set(Frame_SPCR_dict['Frame4']):
		worksheet.write(row,col+6,1)
	else:
		worksheet.write(row,col+6,0)
	if proto_frame5['%s' % pixet][2:35].seq in set(Frame_SPCR_dict['Frame5']):
		worksheet.write(row,col+7,1)
	else:
		worksheet.write(row,col+7,0)
	row += 1

worksheet2 = workbook.add_worksheet('Scores_by_Frame')
#headings
worksheet2.write(0,1,'Pixets Correctly Identified (Percentage)', bold)
row = 1
col = 0
for i in range(1,6):
	worksheet2.write(row,col,'Frame %s' % i)
	row +=1
#data
worksheet2.write(1,1,Frame1_score)
worksheet2.write(2,1,Frame2_score)
worksheet2.write(3,1,Frame3_score)
worksheet2.write(4,1,Frame4_score)
worksheet2.write(5,1,Frame5_score)
workbook.close()


#make the images
for frame in range(1,6): #1,6 after debugging
	Seq_dict = {}
	frame_spcrs = Frame_SPCR_dict['Frame%s' % frame]

	im = Image.new('L', (36,26), 20)
	pix = im.load()

	for record in frame_spcrs:
		PixEt_ID = record[28:33]
		PixEt_ID_num = seq_to_bin_num(PixEt_ID)
		if PixEt_ID_num == 126:		#this fixes disallowed Pixets because of CTT
			PixEt = 8
		elif PixEt_ID_num == 127:
			PixEt = 40
		elif PixEt_ID_num == 128:
			PixEt = 120
		else:
			PixEt = PixEt_ID_num
		Seq_dict[PixEt] = record[1:28]

	PixEt = 1
	for row in range(0,13):
		for i in range(0,8):
			if PixEt in Seq_dict:
				current = Seq_dict[PixEt]
				pix[i,row] = trip_to_num[current[0:3]]
				pix[i+12,row+13] = trip_to_num[current[3:6]]
				pix[i+24,row] = trip_to_num[current[6:9]]
				pix[i+4,row+13] = trip_to_num[current[9:12]]
				pix[i+16,row] = trip_to_num[current[12:15]]
				pix[i+28,row+13] = trip_to_num[current[15:18]]
				pix[i+8,row] = trip_to_num[current[18:21]]
				pix[i+20,row+13] = trip_to_num[current[21:24]]
				if PixEt in [1,2,3,4,9,10,11,12,17,18,19,20,25,26,27,28,33,34,35,36,41,42,43,44,49,50,51,52,57,58,59,60,65,66,67,68,73,74,75,76,81,82,83,84,89,90,91,92,97,98,99,100]:
					pix[i+32,row] = trip_to_num[current[24:27]]
				else:
					pix[i-4,row+13] = trip_to_num[current[24:27]]

			PixEt +=1	

	png.fromarray(np.asarray(im),'L;8').save("%s/reconstructed_image_frame%s.png" % (Data_Path,frame))


#Called_vs_Right_analysis
workbook = xlsxwriter.Workbook('%s/Called_vs_Right.xlsx' % Data_Path)
bold = workbook.add_format({'bold': True})

worksheet1 = workbook.add_worksheet('Frame1')
worksheet1.write(0,0,'Pixet',bold)
worksheet1.write(0,1,'Intended',bold)
worksheet1.write(0,2,'Called',bold)
worksheet1.write(0,3,'Correct?',bold)
row = 1
col = 0
Seq_dict = {}
for record in Frame_SPCR_dict['Frame1']:
	PixEt_ID = record[28:33]
	PixEt_ID_num = seq_to_bin_num(PixEt_ID)
	if PixEt_ID_num == 126:		#this fixes disallowed Pixets because of CTT
		PixEt = 8
	elif PixEt_ID_num == 127:
		PixEt = 40
	elif PixEt_ID_num == 128:
		PixEt = 120
	else:
		PixEt = PixEt_ID_num
	Seq_dict[PixEt] = record[1:28]
for key in proto_frame1:
	worksheet1.write(row,col,key)
	worksheet1.write(row,col+1,str(proto_frame1[key][3:30].seq))
	if int(key) in Seq_dict:
		worksheet1.write(row,col+2,str(Seq_dict[int(key)]))
		if proto_frame1[key][3:30].seq == Seq_dict[int(key)]:
			worksheet1.write(row,col+3,'1')
		else:
			worksheet1.write(row,col+3,'0')
	else:
		worksheet1.write(row,col+2,'none')
		worksheet1.write(row,col+3,'0')
	row +=1

worksheet2 = workbook.add_worksheet('Frame2')
worksheet2.write(0,0,'Pixet',bold)
worksheet2.write(0,1,'Intended',bold)
worksheet2.write(0,2,'Called',bold)
worksheet2.write(0,3,'Correct?',bold)
row = 1
col = 0
Seq_dict = {}
for record in Frame_SPCR_dict['Frame2']:
	PixEt_ID = record[28:33]
	PixEt_ID_num = seq_to_bin_num(PixEt_ID)
	if PixEt_ID_num == 126:		#this fixes disallowed Pixets because of CTT
		PixEt = 8
	elif PixEt_ID_num == 127:
		PixEt = 40
	elif PixEt_ID_num == 128:
		PixEt = 120
	else:
		PixEt = PixEt_ID_num
	Seq_dict[PixEt] = record[1:28]
for key in proto_frame2:
	worksheet2.write(row,col,key)
	worksheet2.write(row,col+1,str(proto_frame2[key][3:30].seq))
	if int(key) in Seq_dict:
		worksheet2.write(row,col+2,str(Seq_dict[int(key)]))
		if proto_frame2[key][3:30].seq == Seq_dict[int(key)]:
			worksheet2.write(row,col+3,'1')
		else:
			worksheet2.write(row,col+3,'0')
	else:
		worksheet2.write(row,col+2,'none')
		worksheet2.write(row,col+3,'0')
	row +=1

worksheet3 = workbook.add_worksheet('Frame3')
worksheet3.write(0,0,'Pixet',bold)
worksheet3.write(0,1,'Intended',bold)
worksheet3.write(0,2,'Called',bold)
worksheet3.write(0,3,'Correct?',bold)
row = 1
col = 0
Seq_dict = {}
for record in Frame_SPCR_dict['Frame3']:
	PixEt_ID = record[28:33]
	PixEt_ID_num = seq_to_bin_num(PixEt_ID)
	if PixEt_ID_num == 126:		#this fixes disallowed Pixets because of CTT
		PixEt = 8
	elif PixEt_ID_num == 127:
		PixEt = 40
	elif PixEt_ID_num == 128:
		PixEt = 120
	else:
		PixEt = PixEt_ID_num
	Seq_dict[PixEt] = record[1:28]
for key in proto_frame3:
	worksheet3.write(row,col,key)
	worksheet3.write(row,col+1,str(proto_frame3[key][3:30].seq))
	if int(key) in Seq_dict:
		worksheet3.write(row,col+2,str(Seq_dict[int(key)]))
		if proto_frame3[key][3:30].seq == Seq_dict[int(key)]:
			worksheet3.write(row,col+3,'1')
		else:
			worksheet3.write(row,col+3,'0')
	else:
		worksheet3.write(row,col+2,'none')
		worksheet3.write(row,col+3,'0')
	row +=1

worksheet4 = workbook.add_worksheet('Frame4')
worksheet4.write(0,0,'Pixet',bold)
worksheet4.write(0,1,'Intended',bold)
worksheet4.write(0,2,'Called',bold)
worksheet4.write(0,3,'Correct?',bold)
row = 1
col = 0
Seq_dict = {}
for record in Frame_SPCR_dict['Frame4']:
	PixEt_ID = record[28:33]
	PixEt_ID_num = seq_to_bin_num(PixEt_ID)
	if PixEt_ID_num == 126:		#this fixes disallowed Pixets because of CTT
		PixEt = 8
	elif PixEt_ID_num == 127:
		PixEt = 40
	elif PixEt_ID_num == 128:
		PixEt = 120
	else:
		PixEt = PixEt_ID_num
	Seq_dict[PixEt] = record[1:28]
for key in proto_frame4:
	worksheet4.write(row,col,key)
	worksheet4.write(row,col+1,str(proto_frame4[key][3:30].seq))
	if int(key) in Seq_dict:
		worksheet4.write(row,col+2,str(Seq_dict[int(key)]))
		if proto_frame4[key][3:30].seq == Seq_dict[int(key)]:
			worksheet4.write(row,col+3,'1')
		else:
			worksheet4.write(row,col+3,'0')
	else:
		worksheet4.write(row,col+2,'none')
		worksheet4.write(row,col+3,'0')
	row +=1

worksheet5 = workbook.add_worksheet('Frame5')
worksheet5.write(0,0,'Pixet',bold)
worksheet5.write(0,1,'Intended',bold)
worksheet5.write(0,2,'Called',bold)
worksheet5.write(0,3,'Correct?',bold)
row = 1
col = 0
Seq_dict = {}
for record in Frame_SPCR_dict['Frame5']:
	PixEt_ID = record[28:33]
	PixEt_ID_num = seq_to_bin_num(PixEt_ID)
	if PixEt_ID_num == 126:		#this fixes disallowed Pixets because of CTT
		PixEt = 8
	elif PixEt_ID_num == 127:
		PixEt = 40
	elif PixEt_ID_num == 128:
		PixEt = 120
	else:
		PixEt = PixEt_ID_num
	Seq_dict[PixEt] = record[1:28]
for key in proto_frame5:
	worksheet5.write(row,col,key)
	worksheet5.write(row,col+1,str(proto_frame5[key][3:30].seq))
	if int(key) in Seq_dict:
		worksheet5.write(row,col+2,str(Seq_dict[int(key)]))
		if proto_frame5[key][3:30].seq == Seq_dict[int(key)]:
			worksheet5.write(row,col+3,'1')
		else:
			worksheet5.write(row,col+3,'0')
	else:
		worksheet5.write(row,col+2,'none')
		worksheet5.write(row,col+3,'0')
	row +=1

workbook.close()