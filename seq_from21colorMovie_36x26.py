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
Trip_bin = {'0000': 'TGA', '1000': 'AGA', '0100': 'TAA', '1100': 'CAC', '0010': 'GAC', '1010': 'AGC', '0110': 'GGA', '1110': 'TGG'}

exceptions_dict = {}

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
	num_bin_rev = str(num_bin)[::-1]
	nuc = [Seq_bin[num_bin_rev[0:2]]]
	nuc.append(Seq_bin[num_bin_rev[2:4]])
	nuc.append(Trip_bin[num_bin_rev[4:8]])
	PixEt_nuc = ''.join(nuc)
	return PixEt_nuc

def select_triplet(spcr,color,PxEt):
	full_spcr = ''.join(spcr)
	GC_spcr = ((full_spcr.count('G')+full_spcr.count('C'))/float(len(full_spcr)))*100
	trip1 = trip_dict[color][0]
	trip2 = trip_dict[color][1]
	trip3 = trip_dict[color][2]
	GC_trip1 = ((trip1.count('G')+trip1.count('C'))/float(len(trip1)))*100
	GC_trip2 = ((trip2.count('G')+trip2.count('C'))/float(len(trip2)))*100
	GC_trip3 = ((trip3.count('G')+trip3.count('C'))/float(len(trip3)))*100
	sorted_trips = [(trip1,GC_trip1),(trip2,GC_trip2),(trip3,GC_trip3)]
	random.shuffle(sorted_trips)
	sorted_trips.sort(key=lambda tup: tup[1])
	if GC_spcr <= 35:
		ordered_trips = [sorted_trips[2][0],sorted_trips[1][0],sorted_trips[0][0]]
	elif 35 < GC_spcr < 50:
		ordered_trips = [sorted_trips[1][0],sorted_trips[2][0],sorted_trips[0][0]]
	elif 50 <= GC_spcr < 65:
		ordered_trips = [sorted_trips[1][0],sorted_trips[0][0],sorted_trips[2][0]]
	elif GC_spcr >= 65:
		ordered_trips = [sorted_trips[0][0],sorted_trips[1][0],sorted_trips[2][0]]
	poss1 = full_spcr[len(full_spcr)-3:len(full_spcr)]+ordered_trips[0]
	poss2 = full_spcr[len(full_spcr)-3:len(full_spcr)]+ordered_trips[1]
	poss3 = full_spcr[len(full_spcr)-3:len(full_spcr)]+ordered_trips[2]
	if len(full_spcr) == 3:
		if 'CTT' not in poss1 and consecutive(poss1) == False:
			triplet = ordered_trips[0]
		elif 'CTT' not in poss2 and consecutive(poss2) == False:
			triplet = ordered_trips[1]
		elif 'CTT' not in poss3 and consecutive(poss3) == False:
			triplet = ordered_trips[2]
		else:
			triplet = ordered_trips[0]
			exceptions_dict[frame].append(PxEt)
	elif len(full_spcr) > 3:
		if 'AAG' not in poss1 and 'CTT' not in poss1 and consecutive(poss1) == False:
			triplet = ordered_trips[0]
		elif 'AAG' not in poss2 and 'CTT' not in poss2 and consecutive(poss2) == False:
			triplet = ordered_trips[1]
		elif 'AAG' not in poss3 and 'CTT' not in poss3 and consecutive(poss3) == False:
			triplet = ordered_trips[2]
		else:
			triplet = ordered_trips[0]
			exceptions_dict[frame].append(PxEt)
	return triplet

def select_triplet_last(spcr,color,PxEt,PxEt_ID):
	full_spcr = ''.join(spcr)
	GC_spcr = ((full_spcr.count('G')+full_spcr.count('C'))/float(len(full_spcr)))*100
	trip1 = trip_dict[color][0]
	trip2 = trip_dict[color][1]
	trip3 = trip_dict[color][2]
	GC_trip1 = ((trip1.count('G')+trip1.count('C'))/float(len(trip1)))*100
	GC_trip2 = ((trip2.count('G')+trip2.count('C'))/float(len(trip2)))*100
	GC_trip3 = ((trip3.count('G')+trip3.count('C'))/float(len(trip3)))*100
	sorted_trips = [(trip1,GC_trip1),(trip2,GC_trip2),(trip3,GC_trip3)]
	random.shuffle(sorted_trips)
	sorted_trips.sort(key=lambda tup: tup[1])
	if GC_spcr <= 35:
		ordered_trips = [sorted_trips[2][0],sorted_trips[1][0],sorted_trips[0][0]]
	elif 35 < GC_spcr < 50:
		ordered_trips = [sorted_trips[1][0],sorted_trips[2][0],sorted_trips[0][0]]
	elif 50 <= GC_spcr < 65:
		ordered_trips = [sorted_trips[1][0],sorted_trips[0][0],sorted_trips[2][0]]
	elif GC_spcr >= 65:
		ordered_trips = [sorted_trips[0][0],sorted_trips[1][0],sorted_trips[2][0]]
	poss1 = full_spcr[len(full_spcr)-3:len(full_spcr)]+ordered_trips[0]+bin_num_to_seq(PxEt_ID)
	poss2 = full_spcr[len(full_spcr)-3:len(full_spcr)]+ordered_trips[1]+bin_num_to_seq(PxEt_ID)
	poss3 = full_spcr[len(full_spcr)-3:len(full_spcr)]+ordered_trips[2]+bin_num_to_seq(PxEt_ID)
	if 'AAG' not in poss1 and 'CTT' not in poss1 and consecutive(poss1) == False:
		triplet = ordered_trips[0]
	elif 'AAG' not in poss2 and 'CTT' not in poss2 and consecutive(poss2) == False:
		triplet = ordered_trips[1]
	elif 'AAG' not in poss3 and 'CTT' not in poss3 and consecutive(poss3) == False:
		triplet = ordered_trips[2]
	elif 'CTT' not in poss1 and consecutive(poss1) == False:
		triplet = ordered_trips[0]
	elif 'CTT' not in poss2 and consecutive(poss2) == False:
		triplet = ordered_trips[1]
	elif 'CTT' not in poss3 and consecutive(poss3) == False:
		triplet = ordered_trips[2]
	elif 'CTT' not in poss1:
		triplet = ordered_trips[0]
	elif 'CTT' not in poss2:
		triplet = ordered_trips[1]
	elif 'CTT' not in poss3:
		triplet = ordered_trips[2]
	else:
		triplet = ordered_trips[0]
		exceptions_dict[frame].append(PxEt)
	return triplet

def consecutive(string):
	return re.search(r'(.)\1\1\1', string) is not None

"""Run"""
for frame in range(1,6): 
	Seq_list = []
	Seq_list_recs = []
	hp_list = []
	exceptions_dict[frame] = []
	im = Image.open("%s/source%s.tif" % (Data_Path,frame))
	pix = im.load() #for pix[x,y], x starts left and goes right, y starts at top and goes down 
	PixEt = 1
	for row in range(0,13):
		for i in range(0,8):

			if PixEt == 8:		#this fixes disallowed Pixets because of CTT
				PixEt_ID = 126
			elif PixEt == 40:
				PixEt_ID = 127
			elif PixEt == 120:
				PixEt_ID = 128
			else:
				PixEt_ID = PixEt

			seq = ['AAG']
			seq.append(select_triplet(seq,pix[i,row],PixEt))	#this starts adding the pixel values
			seq.append(select_triplet(seq,pix[i+12,row+13],PixEt))
			seq.append(select_triplet(seq,pix[i+24,row],PixEt))
			seq.append(select_triplet(seq,pix[i+4,row+13],PixEt))
			seq.append(select_triplet(seq,pix[i+16,row],PixEt))
			seq.append(select_triplet(seq,pix[i+28,row+13],PixEt))
			seq.append(select_triplet(seq,pix[i+8,row],PixEt))
			seq.append(select_triplet(seq,pix[i+20,row+13],PixEt))
			if PixEt in [1,2,3,4,9,10,11,12,17,18,19,20,25,26,27,28,33,34,35,36,41,42,43,44,49,50,51,52,57,58,59,60,65,66,67,68,73,74,75,76,81,82,83,84,89,90,91,92,97,98,99,100]:
				seq.append(select_triplet_last(seq,pix[i+32,row],PixEt,PixEt_ID))
			else:
				seq.append(select_triplet_last(seq,pix[i-4,row+13],PixEt,PixEt_ID))				
			seq.append(bin_num_to_seq(PixEt_ID))  #this adds the PixEt identifier

			full_seq = ''.join(seq)
			seq_obj = Seq(full_seq)
			seq_rec = SeqRecord(seq_obj, id='%s' % PixEt, description='ID=%s' % PixEt_ID)

			Seq_list.append(full_seq)
			Seq_list_recs.append(seq_rec)
			PixEt +=1

	#Change to minimal hairpin format
	for seq in Seq_list_recs:
		hp = [str(seq[7:35].seq)]
		hp.append(str(seq[0:30].seq.reverse_complement()))
		full_hp = ''.join(hp)
		hp_list.append(full_hp)

	"""Addtl Output"""
	Seqs_q = open("%s/Protospacer_seqs_frame%s.fasta" % (Data_Path,frame), "w")
	SeqIO.write(Seq_list_recs, Seqs_q, "fasta")
	Seqs_q.close()

	#For IDT ordering, write to excel, 2X 96 well w/ 50 each
	#Output format is as a minimal hairpin
	workbook1 = xlsxwriter.Workbook('%s/IDT_plate1_frame%s.xlsx' % (Data_Path,frame))
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
	row = 0
	col = 1
	for i in range(1,53):
		worksheet.write(row+i,col,i)
	row = 1
	col = 2
	for i in range(0,52):
		worksheet.write(row+i,col,'%s' % hp_list[i])
	workbook1.close()

	workbook2 = xlsxwriter.Workbook('%s/IDT_plate2_frame%s.xlsx' % (Data_Path,frame))
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
	row = 0
	col = 1
	for i in range(1,53):
		worksheet2.write(row+i,col,i+52)
	row = 1
	col = 2
	for i in range(0,52):
		worksheet2.write(row+i,col,'%s' % hp_list[i+52])
	workbook2.close()


