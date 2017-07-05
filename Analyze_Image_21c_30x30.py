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

PixEt_2_ID_dict = {}	#this makes PixEts next to each other have more distinct sequences in the identifier
for i in range(1,101,2):
	PixEt_2_ID_dict[i] = i
for i in range(2,51,2):
	PixEt_2_ID_dict[i] = i+50
for i in range(52,101,2):
	PixEt_2_ID_dict[i] = i-50
ID_2_PixEt_dict = {v:k for k, v in PixEt_2_ID_dict.items()}

#This is the nucleotide triplet to pixel value conversion, numbers depend on the palette used, in this case extracted from the image in ImageJ
trip_dict = {96: ['TTT','CCC','AAA'],			
			 192: ['TCT','CAC','AGA'],
			 136: ['TAT','CGC','ATG'],
			 208: ['TGT','CTA','ACG'],
			 128: ['TTC','CCA','AGG'],
			 200: ['TCC','CAA','GTT'],
			 184: ['TAC','CGA','GCT'],
			 120: ['TGC','CTG','GAT'],
			 176: ['TTA','CCG','GGT'],
			 112: ['TCA','CAG','GTC'],
			 168:['TAA','CGG','GCC'],
			 104:['TGA','ATT','GAC'],
			 160:['TTG','ACT','GGC'],
			 152:['TCG','AAT','GTA'],
			 144:['TAG','AGT','GCA'],
			 15:['TGG','ATC','GAA'],
			 16:['CTT','ACC','GGA'],
			 17:['CCT','AAC','GTG'],
			 18:['CAT','AGC','GCG'],
			 19:['CGT','ATA','GAG'],
			 20:['CTC','ACA','GGG'], 
			 0:['AAG']} 

trip_to_num = dict( (v,k) for k in trip_dict for v in trip_dict[k] )

Seq_dict = {}

im = Image.new('L', (30,30), 20)
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
		if 0 < PixEt_ID_num < 101:
			PixEt = ID_2_PixEt_dict[PixEt_ID_num]
			if PixEt not in Seq_dict:
				Seq_dict[PixEt] = seq[0][5:34]

PixEt = 1
for row in range(0,10):
	for i in range(0,10):
		if PixEt in Seq_dict:
			current = Seq_dict[PixEt]
			pix[i,row] = trip_to_num[current[0:3]]
			pix[i+10,row+10] = trip_to_num[current[3:6]]
			pix[i+20,row+20] = trip_to_num[current[6:9]]
			pix[i+10,row] = trip_to_num[current[9:12]]
			pix[i+20,row+10] = trip_to_num[current[12:15]]
			pix[i,row+20] = trip_to_num[current[15:18]]
			pix[i+20,row] = trip_to_num[current[18:21]]
			pix[i,row+10] = trip_to_num[current[21:24]]
			pix[i+10,row+20] = trip_to_num[current[24:27]]
		PixEt +=1	

"""Output"""
png.fromarray(np.asarray(im),'L;8').save("%s/reconstructed_image.png" % Data_Path)


