"""Import Modules"""
import sys,os
from PIL import Image
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import re
from collections import Counter
import xlsxwriter

"""Arguments"""
image_set = sys.argv[1] #image set, defines folder

"""Globals"""
user_profile = os.environ ['USERPROFILE']
Data_Path = '%s/Image/Image_set_%s' % (user_profile,image_set) #for running locally
Seq_bin = {'00': 'C', '01': 'T', '10': 'A', '11': 'G'}
Seq_list = []
Seq_list_recs = []
hp_list = []
Seq_list_recs_scram = []
PixEt_2_ID_dict = {}	#this makes PixEts next to each other have more distinct sequences in the identifier
for i in range(1,101,2):
	PixEt_2_ID_dict[i] = i
for i in range(2,51,2):
	PixEt_2_ID_dict[i] = i+50
for i in range(52,101,2):
	PixEt_2_ID_dict[i] = i-50
exceptions = []

trip_dict = {0: ['TTT','CCC','AAA'],
			 1: ['TCT','CAC','AGA'],
			 2: ['TAT','CGC','ATG'],
			 3: ['TGT','CTA','ACG'],
			 4: ['TTC','CCA','AGG'],
			 5: ['TCC','CAA','GTT'],
			 6: ['TAC','CGA','GCT'],
			 7: ['TGC','CTG','GAT'],
			 8: ['TTA','CCG','GGT'],
			 9: ['TCA','CAG','GTC'],
			 10:['TAA','CGG','GCC'],
			 11:['TGA','ATT','GAC'],
			 12:['TTG','ACT','GGC'],
			 13:['TCG','AAT','GTA'],
			 14:['TAG','AGT','GCA'],
			 15:['TGG','ATC','GAA'],
			 16:['CTT','ACC','GGA'],
			 17:['CCT','AAC','GTG'],
			 18:['CAT','AGC','GCG'],
			 19:['CGT','ATA','GAG'],
			 20:['CTC','ACA','GGG']} 

"""Defs"""
def bin_num_to_seq(num): 
	num_bin = '{0:08b}'.format(num)  #gives number in 8 digit binary
	nuc = [Seq_bin[num_bin[0:2]]]
	nuc.append(Seq_bin[num_bin[2:4]])
	nuc.append(Seq_bin[num_bin[4:6]])
	nuc.append(Seq_bin[num_bin[6:8]])
	PixEt_nuc = ''.join(nuc)
	return PixEt_nuc

def select_triplet(spcr,color,PxEt):
	full_spcr = ''.join(spcr)
	# print full_spcr
	GC_spcr = ((full_spcr.count('G')+full_spcr.count('C'))/float(len(full_spcr)))*100
	# print GC_spcr
	trip1 = trip_dict[color][0]
	trip2 = trip_dict[color][1]
	trip3 = trip_dict[color][2]
	GC_trip1 = ((trip1.count('G')+trip1.count('C'))/float(len(trip1)))*100
	GC_trip2 = ((trip2.count('G')+trip2.count('C'))/float(len(trip2)))*100
	GC_trip3 = ((trip3.count('G')+trip3.count('C'))/float(len(trip3)))*100
	sorted_trips = [(trip1,GC_trip1),(trip2,GC_trip2),(trip3,GC_trip3)]
	random.shuffle(sorted_trips)
	sorted_trips.sort(key=lambda tup: tup[1])
	# print sorted_trips
	if GC_spcr <= 35:
		ordered_trips = [sorted_trips[2][0],sorted_trips[1][0],sorted_trips[0][0]]
	elif 35 < GC_spcr < 50:
		ordered_trips = [sorted_trips[1][0],sorted_trips[2][0],sorted_trips[0][0]]
	elif 50 <= GC_spcr < 65:
		ordered_trips = [sorted_trips[1][0],sorted_trips[0][0],sorted_trips[2][0]]
	elif GC_spcr >= 65:
		ordered_trips = [sorted_trips[0][0],sorted_trips[1][0],sorted_trips[2][0]]
	# print ordered_trips
	poss1 = full_spcr[len(full_spcr)-3:len(full_spcr)]+ordered_trips[0]
	poss2 = full_spcr[len(full_spcr)-3:len(full_spcr)]+ordered_trips[1]
	poss3 = full_spcr[len(full_spcr)-3:len(full_spcr)]+ordered_trips[2]
	# if consecutive(poss1) == True:
	# 	print 'Fail'
	if 'AAG' not in poss1 and 'CTT' not in poss1 and consecutive(poss1) == False:
		triplet = ordered_trips[0]
	elif 'AAG' not in poss2 and 'CTT' not in poss2 and consecutive(poss2) == False:
		triplet = ordered_trips[1]
	elif 'AAG' not in poss3 and 'CTT' not in poss2 and consecutive(poss3) == False:
		triplet = ordered_trips[2]
	else:
		triplet = ordered_trips[0]
		exceptions.append(PxEt)
	return triplet

def consecutive(string):
	return re.search(r'(.)\1\1\1', string) is not None

def pick_last(spcr):
	full_spcr = ''.join(spcr)
	c = Counter(full_spcr)
	poss1 = full_spcr[len(full_spcr)-3:len(full_spcr)]+c.most_common(4)[3][0]
	poss2 = full_spcr[len(full_spcr)-3:len(full_spcr)]+c.most_common(4)[2][0]
	poss3 = full_spcr[len(full_spcr)-3:len(full_spcr)]+c.most_common(4)[1][0]
	if 'AAG' not in poss1 and 'CTT' not in poss1 and consecutive(poss1) == False:
		final = c.most_common(4)[3][0]
	elif 'AAG' not in poss2 and 'CTT' not in poss2 and consecutive(poss2) == False:
		final = c.most_common(4)[2][0]
	elif 'AAG' not in poss3 and 'CTT' not in poss3 and consecutive(poss3) == False:
		final = c.most_common(4)[1][0]
	return final


"""Run"""
im = Image.open("%s/source.tif" % Data_Path)
pix = im.load() #for pix[x,y], x starts left and goes right, y starts at top and goes down 
PixEt = 1
for row in range(0,10):
	for i in range(0,10):
		PixEt_ID = PixEt_2_ID_dict[PixEt]
		seq = ['AAG']
		seq.append(bin_num_to_seq(PixEt_ID))  #this adds the PixEt identifier, NOTE: change to PixEt ID to make ids of nearby pixEts more different
		seq.append(select_triplet(seq,pix[i,row],PixEt))	#this starts adding the pixel values
		seq.append(select_triplet(seq,pix[i+10,row+10],PixEt))
		seq.append(select_triplet(seq,pix[i+20,row+20],PixEt))
		seq.append(select_triplet(seq,pix[i+10,row],PixEt))
		seq.append(select_triplet(seq,pix[i+20,row+10],PixEt))
		seq.append(select_triplet(seq,pix[i,row+20],PixEt))
		seq.append(select_triplet(seq,pix[i+20,row],PixEt))
		seq.append(select_triplet(seq,pix[i,row+10],PixEt))
		seq.append(select_triplet(seq,pix[i+10,row+20],PixEt))
		seq.append(pick_last(seq))

		full_seq = ''.join(seq)
		seq_scram = [full_seq[2:7]]
		seq_scram.append(''.join(random.sample(full_seq[7:35],len(full_seq[7:35]))))
		scramble = ''.join(seq_scram)
		seq_obj = Seq(full_seq)
		seq_rec = SeqRecord(seq_obj, id='%s' % PixEt, description='ID=%s' % PixEt_ID)
		seq_obj_scram = Seq(scramble)
		seq_rec_scram = SeqRecord(seq_obj_scram, id='%s' % PixEt, description='ID=%s' % PixEt_ID)

		# print PixEt
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

#For IDT ordering, write to excel, 2X 96 well w/ 50 each
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
row = 0
col = 1
for i in range(1,51):
	worksheet.write(row+i,col,i)
row = 1
col = 2
for i in range(0,50):
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
row = 0
col = 1
for i in range(1,51):
	worksheet2.write(row+i,col,i+50)
row = 1
col = 2
for i in range(0,50):
	worksheet2.write(row+i,col,'%s' % hp_list[i+50])
workbook2.close()