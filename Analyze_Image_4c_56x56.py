"""Import Modules"""
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
import os, sys
import csv
import xlsxwriter
from collections import Counter
from PIL import Image
import png

"""Line Arguments"""
run_number = sys.argv[1] #MiSeq run number for Data_Path
condition = sys.argv[2] #file name of sequencing data (minus .fastq)

"""Globals"""
user_profile = os.environ ['USERPROFILE']
Data_Path = '%s/Data Analysis/MS%s/%s_Results' % (user_profile,run_number,condition) #this should be a folder with the fastq in it
Blast_Data_Path = '%s/Blast_Databases' % user_profile  #this should be a folder with NCBI blast databases

unaligned_IDs = []
unaligned_SPCRs = []
unaligned_SPCR_nuc = []
count_alignments = []
Seq_dict = {}

handle = open("%s/Specific_Blast_Database.fasta" % Blast_Data_Path, "rU") #This should contain the bacterial genome, lacI, and the plasmid used to express Cas1+2, in that order
records = list(SeqIO.parse(handle,"fasta",generic_dna))
handle.close()
genome_seq = records[0]
lacI_seq = records[1]
plasmid_seq = records[2]

Seq_bin = {'00': 'C', '01': 'T', '10': 'A', '11': 'G'}
bin_Seq = {'C': '00', 'T': '01', 'A': '10', 'G': '11'}

Seq_sing = {0: 'C', 1: 'T', 2: 'A', 3: 'G'}
# Sing_seq = {'C': 0, 'T': 7, 'A': 3, 'G': 11}		#higher contrast
Sing_seq = {'C': 0, 'T': 6, 'A': 4, 'G': 8}

PixEt_2_ID_dict = {}	#this makes PixEts next to each other have more distinct sequences in the identifier
for i in range(1,113,2):
	PixEt_2_ID_dict[i] = i
for i in range(2,57,2):
	PixEt_2_ID_dict[i] = i+56
for i in range(58,113,2):
	PixEt_2_ID_dict[i] = i-56
ID_2_PixEt_dict = {v:k for k, v in PixEt_2_ID_dict.items()}

Seq_dict = {}

im = Image.new('L', (56,56), 15)
pix = im.load()

"""Defs"""
def get_context(record):
	"""input is a blastn record
		returns info and protospacer context"""
	direction = []
	if not blast_record.alignments:
		unaligned_IDs.append (blast_record.query)

def bin_num_to_seq(num): 
	num_bin = '{0:08b}'.format(num)  #gives number in 8 digit binary
	nuc = [Seq_bin[num_bin[0:2]]]
	nuc.append(Seq_bin[num_bin[2:4]])
	nuc.append(Seq_bin[num_bin[4:6]])
	nuc.append(Seq_bin[num_bin[6:8]])
	PixEt_nuc = ''.join(nuc)
	return PixEt_nuc

def seq_to_bin_num(seq):
	seq_code = [bin_Seq[seq[0]]]
	seq_code.append(bin_Seq[seq[1]])
	seq_code.append(bin_Seq[seq[2]])
	seq_code.append(bin_Seq[seq[3]])
	PixEt_bin = ''.join(seq_code)
	PixEt_int = int(PixEt_bin, 2)
	return PixEt_int

"""Run"""
#First, blast
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
for seq in freq_list:
	if len(seq[0]) == 33:
		PixEt_ID = seq[0][1:5]
		PixEt_ID_num = seq_to_bin_num(PixEt_ID)
		if 0 < PixEt_ID_num < 113:
			PixEt = ID_2_PixEt_dict[PixEt_ID_num]
			if PixEt not in Seq_dict:
				Seq_dict[PixEt] = seq[0][5:34]

PixEt = 1
for row in range(0,14):
	for i in range(0,8):
		if PixEt in Seq_dict:
			current = Seq_dict[PixEt]
			pix[i,row] = Sing_seq[current[0]]
			pix[i,row+28] = Sing_seq[current[1]]
			pix[i+16,row] = Sing_seq[current[2]]
			pix[i+16,row+28] = Sing_seq[current[3]]
			pix[i+32,row] = Sing_seq[current[4]]
			pix[i+32,row+28] = Sing_seq[current[5]]
			pix[i+48,row] = Sing_seq[current[6]]
			pix[i+48,row+28] = Sing_seq[current[7]]
			pix[i+8,row] = Sing_seq[current[8]]
			pix[i+8,row+28] = Sing_seq[current[9]]
			pix[i+24,row] = Sing_seq[current[10]]
			pix[i+24,row+28] = Sing_seq[current[11]]
			pix[i+40,row] = Sing_seq[current[12]]
			pix[i+40,row+28] = Sing_seq[current[13]]

			pix[i,row+14] = Sing_seq[current[14]]
			pix[i,row+42] = Sing_seq[current[15]]
			pix[i+16,row+14] = Sing_seq[current[16]]
			pix[i+16,row+42] = Sing_seq[current[17]]
			pix[i+32,row+14] = Sing_seq[current[18]]
			pix[i+32,row+42] = Sing_seq[current[19]]
			pix[i+48,row+14] = Sing_seq[current[20]]
			pix[i+48,row+42] = Sing_seq[current[21]]
			pix[i+8,row+14] = Sing_seq[current[22]]
			pix[i+8,row+42] = Sing_seq[current[23]]
			pix[i+24,row+14] = Sing_seq[current[24]]
			pix[i+24,row+42] = Sing_seq[current[25]]
			pix[i+40,row+14] = Sing_seq[current[26]]
			pix[i+40,row+42] = Sing_seq[current[27]]

		PixEt +=1	


"""Output"""
png.fromarray(np.asarray(im, np.uint8),'L;4').save("%s/reconstructed_image.png" % Data_Path)

