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
import fuzzysearch

"""Line Arguments"""
run_number = sys.argv[1] #MiSeq run number for Data_Path
condition = sys.argv[2] #file name (minus .fastq) shoud be the only argument following SPCR_blast.py when running script
blast_db = sys.argv[3] #which combo of genome + plasmids to define blast database

"""Globals"""
user_profile = os.environ ['USERPROFILE']
Data_Path = '%s/Data Analysis/MS%s/%s_Results' % (user_profile,run_number,condition) #this should be a folder with the fastq in it
Blast_Data_Path = '%s/Blast_Databases' % user_profile  #this should be a folder with NCBI blast databases
protospacer_BL21_buff_15 = []
protospacer_plasmid_buff_15 = []
All_protospacer_buff_15 = []
unaligned_IDs = []
unaligned_SPCRs = []
unaligned_SPCR_nuc = []
count_alignments = []
E_VALUE_THRESH = 0.00000001
Top_fed_SPCRs = []
Top_fed_SPCRs_counts = []
Fed_pos_dict = {}
Repeat = 'GTGTTCCCCGCGCCAGCGGGGATAAACC'
dist_repeat = 4


handle = open("%s/%s.fasta" % (Blast_Data_Path,blast_db), "rU")
records = list(SeqIO.parse(handle,"fasta",generic_dna))
handle.close()
genome_seq = records[0]
lacI_seq = records[1]
Fplasmid_seq = records[2]
if len(records) > 3:
	Splasmid_seq = records[3]

BL21_genome_line_fow = [0]*len(genome_seq)
BL21_genome_line_rev = [0]*len(genome_seq)
BL21_genome_line_single_fow = [0]*len(genome_seq)
BL21_genome_line_single_rev = [0]*len(genome_seq)
Fplasmid_line_fow = [0]*len(Fplasmid_seq)
Fplasmid_line_rev = [0]*len(Fplasmid_seq)
if len(records) > 3:
	Splasmid_line_fow = [0]*len(Splasmid_seq)
	Splasmid_line_rev = [0]*len(Splasmid_seq)

"""Defs"""
def get_context(record):
	"""input is a blastn record
		returns info and protospacer context"""
	direction = []
	if not blast_record.alignments:
		unaligned_IDs.append (blast_record.query)
	elif blast_record.alignments[0].hsps[0].expect < E_VALUE_THRESH and blast_record.alignments[0].title[16:19] =='gen':
		if blast_record.alignments[0].hsps[0].sbjct_start < blast_record.alignments[0].hsps[0].sbjct_end:
			direction = "forward"
		else:
			direction = "reverse"
		if direction == 'forward':
			protospacer_BL21_buff_15.append (genome_seq[blast_record.alignments[0].hsps[0].sbjct_start-16:blast_record.alignments[0].hsps[0].sbjct_end+15])					
			BL21_genome_line_single_fow[blast_record.alignments[0].hsps[0].sbjct_start-1] += 1
			for i in range(blast_record.alignments[0].hsps[0].sbjct_start-1,blast_record.alignments[0].hsps[0].sbjct_end):
				BL21_genome_line_fow[i] += 1
		elif direction == 'reverse':
			protospacer_BL21_buff_15.append (genome_seq[blast_record.alignments[0].hsps[0].sbjct_end-16:blast_record.alignments[0].hsps[0].sbjct_start+15].reverse_complement())
			BL21_genome_line_single_rev[blast_record.alignments[0].hsps[0].sbjct_end-1] -= 1
			for i in range(blast_record.alignments[0].hsps[0].sbjct_end-1,blast_record.alignments[0].hsps[0].sbjct_start):
				BL21_genome_line_rev[i] -=1
	elif blast_record.alignments[0].hsps[0].expect < E_VALUE_THRESH and blast_record.alignments[0].title[16:19] =='Fpl':
		if blast_record.alignments[0].hsps[0].sbjct_start < blast_record.alignments[0].hsps[0].sbjct_end:
			direction = "forward"
		else:
			direction = "reverse"
		if direction == 'forward':
			protospacer_plasmid_buff_15.append (Fplasmid_seq[blast_record.alignments[0].hsps[0].sbjct_start-16:blast_record.alignments[0].hsps[0].sbjct_end+15])
			for i in range(blast_record.alignments[0].hsps[0].sbjct_start-1,blast_record.alignments[0].hsps[0].sbjct_end):
				Fplasmid_line_fow[i] += 1
		elif direction == 'reverse':
			protospacer_plasmid_buff_15.append (Fplasmid_seq[blast_record.alignments[0].hsps[0].sbjct_end-16:blast_record.alignments[0].hsps[0].sbjct_start+15].reverse_complement())
			for i in range(blast_record.alignments[0].hsps[0].sbjct_end-1,blast_record.alignments[0].hsps[0].sbjct_start):
				Fplasmid_line_rev[i] -=1
	elif blast_record.alignments[0].hsps[0].expect < E_VALUE_THRESH and blast_record.alignments[0].title[16:19] =='Spl':
		if blast_record.alignments[0].hsps[0].sbjct_start < blast_record.alignments[0].hsps[0].sbjct_end:
			direction = "forward"
		else:
			direction = "reverse"
		if direction == 'forward':
			protospacer_plasmid_buff_15.append (Splasmid_seq[blast_record.alignments[0].hsps[0].sbjct_start-16:blast_record.alignments[0].hsps[0].sbjct_end+15])
			for i in range(blast_record.alignments[0].hsps[0].sbjct_start-1,blast_record.alignments[0].hsps[0].sbjct_end):
				Splasmid_line_fow[i] += 1
		elif direction == 'reverse':
			protospacer_plasmid_buff_15.append (Splasmid_seq[blast_record.alignments[0].hsps[0].sbjct_end-16:blast_record.alignments[0].hsps[0].sbjct_start+15].reverse_complement())
			for i in range(blast_record.alignments[0].hsps[0].sbjct_end-1,blast_record.alignments[0].hsps[0].sbjct_start):
				Splasmid_line_rev[i] -=1


"""Run"""
#First, blast
blastn_cline = NcbiblastnCommandline(task='blastn-short', query="%s/new_SPCRs_seqs.fasta" % Data_Path, db="%s/%s" % (Blast_Data_Path,blast_db), out="%s/blastn_output.xml" % Data_Path, evalue=E_VALUE_THRESH, outfmt=5)
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
for seq in c.most_common(20):
	if 30 < len(seq[0]) < 35:
		Top_fed_SPCRs.append (seq[0])
		Top_fed_SPCRs_counts.append (seq[1])
#Figure out where the top fed spacers are, not super helpful
for seq in Top_fed_SPCRs:
	Fed_pos_dict[seq] = [0,0,0] #will get counts of time in the first second and third position

"""Output"""
print "Number of Mapped Genomic SPCRs:", len(protospacer_BL21_buff_15)
print "Number of Mapped Plasmid SPCRs:", len(protospacer_plasmid_buff_15)
print "Unaligned SPCRS:", len(unaligned_SPCRs)

#Write protoSPCRs to a fasta file (for pLOGO)
Genomic_protoSPCRs_q = open("%s/Genomic_protoSPCRs_seqs.fasta" % Data_Path, "w")
SeqIO.write(protospacer_BL21_buff_15, Genomic_protoSPCRs_q, "fasta")
Genomic_protoSPCRs_q.close()

Plasmid_protoSPCRs_q = open("%s/Plasmid_protoSPCRs_seqs.fasta" % Data_Path, "w")
SeqIO.write(protospacer_plasmid_buff_15, Plasmid_protoSPCRs_q, "fasta")
Plasmid_protoSPCRs_q.close()

All_protospacer_buff_15 = protospacer_plasmid_buff_15 + protospacer_BL21_buff_15
All_protoSPCRs_q = open("%s/All_protoSPCRs_seqs.fasta" % Data_Path, "w")
SeqIO.write(All_protospacer_buff_15, All_protoSPCRs_q, "fasta")
All_protoSPCRs_q.close()

unaligned_SPCR_reads = open("%s/Unaligned_SPCR_seqs.fasta" % Data_Path, "w")
SeqIO.write(unaligned_SPCRs, unaligned_SPCR_reads, "fasta")
unaligned_SPCR_reads.close()

#Write a excel file with all the relevant data
workbook = xlsxwriter.Workbook('%s/Blast_analysis.xlsx' % Data_Path)
worksheet = workbook.add_worksheet()
#Add titles (row, col: zero referenced)
worksheet.write(0,1,condition)
worksheet.write(0,2,'# Pos1')
worksheet.write(0,3,'# Pos2')
worksheet.write(0,4,'# Pos3')
worksheet.write(1,0,'Number of Mapped Genomic SPCRs:')
worksheet.write(2,0,'Number of Mapped Plasmid SPCRs:')
worksheet.write(3,0,'Number of Unaligned SPCRs (fed):')
row = 4
for i in range (0,len(Top_fed_SPCRs)):
	worksheet.write(row+i,0,'%s' % Top_fed_SPCRs[i])
#Add data
worksheet.write(1,1, len(protospacer_BL21_buff_15))
worksheet.write(2,1, len(protospacer_plasmid_buff_15))
worksheet.write(3,1, len(unaligned_SPCRs))
row = 4
for i in range (0,len(Top_fed_SPCRs_counts)):
	worksheet.write(row+i,1,'%s' % Top_fed_SPCRs_counts[i])
row = 4
col = 2
for i in range (0,len(Top_fed_SPCRs_counts)):
	worksheet.write(row+i,col,'%s' % Fed_pos_dict[Top_fed_SPCRs[i]][0])
	worksheet.write(row+i,col+1,'%s' % Fed_pos_dict[Top_fed_SPCRs[i]][1])
	worksheet.write(row+i,col+2,'%s' % Fed_pos_dict[Top_fed_SPCRs[i]][2])
workbook.close()

