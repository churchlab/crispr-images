"""Import Modules"""
import sys,os
from PIL import Image
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import xlsxwriter

"""Arguments"""
image_set = sys.argv[1] #image set, defines folder


"""Globals"""
user_profile = os.environ ['USERPROFILE']
Data_Path = '%s/Image/Image_set_%s' % (user_profile,image_set) #for running locally
Seq_bin = {'00': 'C', '01': 'T', '10': 'A', '11': 'G'}
Seq_sing = {0: 'C', 1: 'T', 2: 'A', 3: 'G'}
Seq_list = []
Seq_list_recs = []
hp_list = []
Seq_list_recs_scram = []
PixEt_2_ID_dict = {}	#this makes PixEts next to each other have more distinct sequences in the identifier
for i in range(1,113,2):
	PixEt_2_ID_dict[i] = i
for i in range(2,57,2):
	PixEt_2_ID_dict[i] = i+56
for i in range(58,113,2):
	PixEt_2_ID_dict[i] = i-56

"""Defs"""
def bin_num_to_seq(num): 
	num_bin = '{0:08b}'.format(num)  #gives number in 8 digit binary
	nuc = [Seq_bin[num_bin[0:2]]]
	nuc.append(Seq_bin[num_bin[2:4]])
	nuc.append(Seq_bin[num_bin[4:6]])
	nuc.append(Seq_bin[num_bin[6:8]])
	PixEt_nuc = ''.join(nuc)
	return PixEt_nuc

"""Run"""
im = Image.open("%s/source.tif" % Data_Path)
pix = im.load() #for pix[x,y], x starts left and goes right, y starts at top and goes down 
PixEt = 1
for row in range(0,14):
	for i in range(0,8):
		PixEt_ID = PixEt_2_ID_dict[PixEt]
		seq = ['AAG']
		seq.append(bin_num_to_seq(PixEt_ID))  #this adds the PixEt identifier, NOTE: change to PixEt ID to count by 2 to make ids more different
		seq.append(Seq_sing[pix[i,row]])	#this starts adding the pixel values
		seq.append(Seq_sing[pix[i,row+28]])
		seq.append(Seq_sing[pix[i+16,row]])
		seq.append(Seq_sing[pix[i+16,row+28]])
		seq.append(Seq_sing[pix[i+32,row]])
		seq.append(Seq_sing[pix[i+32,row+28]])
		seq.append(Seq_sing[pix[i+48,row]])
		seq.append(Seq_sing[pix[i+48,row+28]])

		seq.append(Seq_sing[pix[i+8,row]])
		seq.append(Seq_sing[pix[i+8,row+28]])
		seq.append(Seq_sing[pix[i+24,row]])
		seq.append(Seq_sing[pix[i+24,row+28]])
		seq.append(Seq_sing[pix[i+40,row]])
		seq.append(Seq_sing[pix[i+40,row+28]])

		seq.append(Seq_sing[pix[i,row+14]])	
		seq.append(Seq_sing[pix[i,row+42]])
		seq.append(Seq_sing[pix[i+16,row+14]])
		seq.append(Seq_sing[pix[i+16,row+42]])
		seq.append(Seq_sing[pix[i+32,row+14]])
		seq.append(Seq_sing[pix[i+32,row+42]])
		seq.append(Seq_sing[pix[i+48,row+14]])
		seq.append(Seq_sing[pix[i+48,row+42]])

		seq.append(Seq_sing[pix[i+8,row+14]])
		seq.append(Seq_sing[pix[i+8,row+42]])
		seq.append(Seq_sing[pix[i+24,row+14]])
		seq.append(Seq_sing[pix[i+24,row+42]])
		seq.append(Seq_sing[pix[i+40,row+14]])
		seq.append(Seq_sing[pix[i+40,row+42]])	

		full_seq = ''.join(seq)
		seq_scram = [full_seq[2:7]]
		seq_scram.append(''.join(random.sample(full_seq[7:35],len(full_seq[7:35]))))
		scramble = ''.join(seq_scram)
		seq_obj = Seq(full_seq)
		seq_rec = SeqRecord(seq_obj, id='%s' % PixEt, description='ID=%s' % PixEt_ID)
		seq_obj_scram = Seq(scramble)
		seq_rec_scram = SeqRecord(seq_obj_scram, id='%s' % PixEt, description='ID=%s' % PixEt_ID)

		print PixEt
		print full_seq

		Seq_list.append(full_seq)
		Seq_list_recs.append(seq_rec)
		Seq_list_recs_scram.append(seq_rec_scram)
		PixEt +=1

#Change to minimal hairpin format
for seq in Seq_list_recs:
	hp = [str(seq[7:35].seq)]
	hp.append(str(seq[0:30].seq.reverse_complement()))
	full_hp = ''.join(hp)
	hp_list.append(full_hp)

"""Addtl Output"""
Seqs_q = open("%s/Protospacer_seqs.fasta" % (Data_Path), "w")
SeqIO.write(Seq_list_recs, Seqs_q, "fasta")
Seqs_q.close()

Seqs_scram_q = open("%s/Protospacer_seqs_scrambled.fasta" % (Data_Path), "w")
SeqIO.write(Seq_list_recs_scram, Seqs_scram_q, "fasta")
Seqs_scram_q.close()

#For IDT ordering, write to excel, 2X 96 well w/ 56 each
#Output format is as a minimal hairpin
workbook1 = xlsxwriter.Workbook('%s/IDT_plate1.xlsx' % Data_Path)
worksheet = workbook1.add_worksheet()
bold = workbook1.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,0,'Well Position',bold)
worksheet.write(0,1,'Name',bold)
worksheet.write(0,2,'Sequence',bold)
Well_letters = ['A','B','C','D']
row = 1
col = 0
for letter in Well_letters:
	for i in range(1,13):
		worksheet.write(row,col,'%s%s' % (letter,i))
		row += 1
worksheet.write(row,col,'E1')
worksheet.write(row+1,col,'E2')
worksheet.write(row+2,col,'E3')
worksheet.write(row+3,col,'E4')
worksheet.write(row+4,col,'E5')
worksheet.write(row+5,col,'E6')
worksheet.write(row+6,col,'E7')
worksheet.write(row+7,col,'E8')
row = 0
col = 1
for i in range(1,57):
	worksheet.write(row+i,col,i)
row = 1
col = 2
for i in range(0,56):
	worksheet.write(row+i,col,'%s' % hp_list[i])
workbook1.close()

workbook2 = xlsxwriter.Workbook('%s/IDT_plate2.xlsx' % Data_Path)
worksheet2 = workbook2.add_worksheet()
bold = workbook2.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet2.write(0,0,'Well Position',bold)
worksheet2.write(0,1,'Name',bold)
worksheet2.write(0,2,'Sequence',bold)
Well_letters = ['A','B','C','D']
row = 1
col = 0
for letter in Well_letters:
	for i in range(1,13):
		worksheet2.write(row,col,'%s%s' % (letter,i))
		row += 1
worksheet2.write(row,col,'E1')
worksheet2.write(row+1,col,'E2')
worksheet2.write(row+2,col,'E3')
worksheet2.write(row+3,col,'E4')
worksheet2.write(row+4,col,'E5')
worksheet2.write(row+5,col,'E6')
worksheet2.write(row+6,col,'E7')
worksheet2.write(row+7,col,'E8')
row = 0
col = 1
for i in range(1,57):
	worksheet2.write(row+i,col,i+56)
row = 1
col = 2
for i in range(0,56):
	worksheet2.write(row+i,col,'%s' % hp_list[i+56])
workbook2.close()