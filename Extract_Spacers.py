"""Import Modules"""
import fuzzysearch
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import numpy
import xlsxwriter
import os, sys
from difflib import SequenceMatcher

"""File"""
run_number = sys.argv[1] #MiSeq run number for Data_Path
condition = sys.argv[2] #file name of sequencing data (minus .fastq)

"""Globals"""
write_binned_whole_reads_to_fastq = 'yes'  #can be helpful for debugging, necessary for recording analysis
user_profile = os.environ ['USERPROFILE']
Data_Path = '%s/Data Analysis/MS%s' % (user_profile,run_number)     #this should be a folder with the fastq in it
fastq_reads = '%s/%s.fastq' % (Data_Path,condition)
Repeat = 'GTGTTCCCCGCGCCAGCGGGGATAAACC'
Old_SPCR1 = 'GAGCACAAATATCATCGCTCAAACCACTTACGG'
Old_SPCR2 = 'GCCTCGCTGTAAATTCCAAAAACGATCTCTATA'
Old_SPCR3 = 'GACAGTACCGGAGTTTGACGGTGCCAACGGCGG'
Old_SPCR4 = 'GACAATCAGGGAACGATTGTTGACACTGTAAAA'
#Sets the overall fuzziness
dist_repeat = 4
dist_oldSPCRs = 5
wildtype_sequences_three_read = []
wildtype_sequences_two_read = []
wildtype_sequences_one_read = []
single_expansion_sequences_three_read_pos_one = []
single_expansion_sequences_three_read_pos_two = []
single_expansion_sequences_two_read_pos_one = []
single_expansion_sequences_two_read_pos_two = []
single_expansion_sequences_one_read = []
single_replacement_sequences_two_read_pos_one = []
double_expansion_sequences_two_read = []
double_expansion_sequences_three_read = []
triple_expansion_sequences_three_read = []
weird_sequences_three_read = []
weird_sequences_two_read = []
weird_sequences_one_read = []
weird_sequences_none_read = []
new_SPCRs = []
new_SPCRs_just_seqs = []
unique_new_SPCRs = []
Unaligned_SPCRs = []
new_SPCR_lengths = []
number_nonwierd_reads = []
nonfed_lengths = []
fed_lengths = []

"""Defs"""
def get_spcrs(sequence):
    last_rep = fuzzysearch.find_near_matches(Repeat[0:15], sequence.seq, max_l_dist=3)
    results = fuzzysearch.find_near_matches(Repeat, sequence.seq, max_l_dist=6)
    if len(results) == 3 and len(last_rep) >= 1 and last_rep[len(last_rep)-1].start > results[len(results)-1].start:
        spacer_list = [sequence[results[0].end:results[1].start]]
        spacer_list.append (sequence[results[1].end:results[2].start])
        spacer_list.append (sequence[results[2].end:last_rep[len(last_rep)-1].start])
    elif len(results) == 2 and len(last_rep) >= 1 and last_rep[len(last_rep)-1].start > results[len(results)-1].start:
        spacer_list = [sequence[results[0].end:results[1].start]]
        spacer_list.append (sequence[results[1].end:last_rep[len(last_rep)-1].start])
    elif len(results) == 1 and len(last_rep) >= 1 and last_rep[len(last_rep)-1].start > results[len(results)-1].start:
        spacer_list = [sequence[results[0].end:last_rep[len(last_rep)-1].start]]

    elif len(results) == 4:
        spacer_list = [sequence[results[0].end:results[1].start]]
        spacer_list.append (sequence[results[1].end:results[2].start])
        spacer_list.append (sequence[results[2].end:results[3].start])
    elif len(results) == 3:
        spacer_list = [sequence[results[0].end:results[1].start]]
        spacer_list.append (sequence[results[1].end:results[2].start])
    elif len(results) == 2:
        spacer_list = [sequence[results[0].end:results[1].start]]
    else:
        spacer_list = []
    return spacer_list

def not_existing(spacer):
    if SequenceMatcher(None,Old_SPCR1,spacer).ratio() < 0.83 and SequenceMatcher(None,Old_SPCR2,spacer).ratio() < 0.83 and SequenceMatcher(None,Old_SPCR3,spacer).ratio() < 0.83 and SequenceMatcher(None,Old_SPCR4,spacer).ratio() < 0.83 and SequenceMatcher(None,Repeat,spacer).ratio() < 0.83:
        return True
    else: return False

"""Run"""
#Create Results folder
newpath = ((r'%s/%s_Results') % (Data_Path,condition))
if not os.path.exists(newpath): os.makedirs(newpath)
#Read in fastq
for seq_record in SeqIO.parse(fastq_reads, "fastq"):
    spacer_list = get_spcrs(seq_record)
    #reads that have three clean spacers
    if len(spacer_list) == 3:
        if SequenceMatcher(None,Old_SPCR1,spacer_list[0]).ratio() > 0.83:
            if SequenceMatcher(None,Old_SPCR2,spacer_list[1]).ratio() > 0.83:
                if SequenceMatcher(None,Old_SPCR3,spacer_list[2]).ratio() > 0.83:
                    wildtype_sequences_three_read.append(seq_record)                             
        elif not_existing(spacer_list[0]) == True:
            if not_existing(spacer_list[1]) == True:
                if not_existing(spacer_list[2]) == True:
                    triple_expansion_sequences_three_read.append(seq_record)
                    spacer_list[1].id+='_2'
                    spacer_list[2].id+='_3'
                    if len(spacer_list[0]) < 61:
                        new_SPCRs.append (spacer_list[0])
                    if len(spacer_list[1]) < 61:
                        new_SPCRs.append (spacer_list[1])
                    if len(spacer_list[2]) < 61:
                        new_SPCRs.append (spacer_list[2])
                else: 
                    double_expansion_sequences_three_read.append(seq_record)
                    spacer_list[1].id+='_2'
                    if len(spacer_list[0]) < 61:
                        new_SPCRs.append (spacer_list[0])
                    if len(spacer_list[1]) < 61:
                        new_SPCRs.append (spacer_list[1])                   
            else:
                single_expansion_sequences_three_read_pos_one.append(seq_record)
                if len(spacer_list[0]) < 61:
                    new_SPCRs.append (spacer_list[0])

        elif not_existing(spacer_list[0]) == False and not_existing(spacer_list[1]) == True:
            single_expansion_sequences_three_read_pos_two.append(seq_record) #not taking spacers from this
        else:
            weird_sequences_three_read.append(seq_record)

    if len(spacer_list) == 2:
        if SequenceMatcher(None,Old_SPCR1,spacer_list[0]).ratio() > 0.83:
            if SequenceMatcher(None,Old_SPCR2,spacer_list[1]).ratio() > 0.83:
                wildtype_sequences_two_read.append(seq_record)
        elif not_existing(spacer_list[0]) == True:
            if not_existing(spacer_list[1]) == True:
                double_expansion_sequences_two_read.append(seq_record)
                spacer_list[1].id+='_2'
                if len(spacer_list[0]) < 61:
                    new_SPCRs.append (spacer_list[0])
                if len(spacer_list[1]) < 61:
                    new_SPCRs.append (spacer_list[1])                                  
            else:
                single_expansion_sequences_two_read_pos_one.append(seq_record)
                if len(spacer_list[0]) < 61:
                    new_SPCRs.append (spacer_list[0])                
        elif not_existing(spacer_list[0]) == False and not_existing(spacer_list[1]) == True:
            single_expansion_sequences_two_read_pos_two.append(seq_record) #not taking spacers from this
        else:
            weird_sequences_two_read.append(seq_record)

    if len(spacer_list) == 1:
        if SequenceMatcher(None,Old_SPCR1,spacer_list[0]).ratio() > 0.83:
            wildtype_sequences_one_read.append(seq_record)
        elif not_existing(spacer_list[0]) == True:
            single_expansion_sequences_one_read.append(seq_record)
            if len(spacer_list[0]) < 61:
                new_SPCRs.append (spacer_list[0])                        
        else:
            weird_sequences_one_read.append(seq_record)

    if len(spacer_list) == 0:
        weird_sequences_none_read.append(seq_record)


#Dump the data back into fastq files
if write_binned_whole_reads_to_fastq == 'yes':
    wildtype_sequences_three_read_q = open("%s/%s_Results/wildtype_sequences_three_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(wildtype_sequences_three_read, wildtype_sequences_three_read_q, "fastq")
    wildtype_sequences_three_read_q.close()

    wildtype_sequences_two_read_q = open("%s/%s_Results/wildtype_sequences_two_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(wildtype_sequences_two_read, wildtype_sequences_two_read_q, "fastq")
    wildtype_sequences_two_read_q.close()

    wildtype_sequences_one_read_q = open("%s/%s_Results/wildtype_sequences_one_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(wildtype_sequences_one_read, wildtype_sequences_one_read_q, "fastq")
    wildtype_sequences_one_read_q.close()

    single_expansion_sequences_two_read_pos_one_q = open("%s/%s_Results/single_expansion_sequences_two_read_pos_one_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(single_expansion_sequences_two_read_pos_one, single_expansion_sequences_two_read_pos_one_q, "fastq")
    single_expansion_sequences_two_read_pos_one_q.close()

    single_expansion_sequences_two_read_pos_two_q = open("%s/%s_Results/single_expansion_sequences_two_read_pos_two_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(single_expansion_sequences_two_read_pos_two, single_expansion_sequences_two_read_pos_two_q, "fastq")
    single_expansion_sequences_two_read_pos_two_q.close()

    single_expansion_sequences_three_read_pos_one_q = open("%s/%s_Results/single_expansion_sequences_three_read_pos_one_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(single_expansion_sequences_three_read_pos_one, single_expansion_sequences_three_read_pos_one_q, "fastq")
    single_expansion_sequences_three_read_pos_one_q.close()

    single_expansion_sequences_three_read_pos_two_q = open("%s/%s_Results/single_expansion_sequences_three_read_pos_two_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(single_expansion_sequences_three_read_pos_two, single_expansion_sequences_three_read_pos_two_q, "fastq")
    single_expansion_sequences_three_read_pos_two_q.close()

    single_expansion_sequences_one_read_q = open("%s/%s_Results/single_expansion_sequences_one_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(single_expansion_sequences_one_read, single_expansion_sequences_one_read_q, "fastq")
    single_expansion_sequences_one_read_q.close()

    double_expansion_sequences_two_read_q = open("%s/%s_Results/double_expansion_sequences_two_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(double_expansion_sequences_two_read, double_expansion_sequences_two_read_q, "fastq")
    double_expansion_sequences_two_read_q.close()

    double_expansion_sequences_three_read_q = open("%s/%s_Results/double_expansion_sequences_three_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(double_expansion_sequences_three_read, double_expansion_sequences_three_read_q, "fastq")
    double_expansion_sequences_three_read_q.close()

    triple_expansion_sequences_three_read_q = open("%s/%s_Results/triple_expansion_sequences_three_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(triple_expansion_sequences_three_read, triple_expansion_sequences_three_read_q, "fastq")
    triple_expansion_sequences_three_read_q.close()

    weird_sequences_three_read_q = open("%s/%s_Results/weird_sequences_three_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(weird_sequences_three_read, weird_sequences_three_read_q, "fastq")
    weird_sequences_three_read_q.close()

    weird_sequences_two_read_q = open("%s/%s_Results/weird_sequences_two_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(weird_sequences_two_read, weird_sequences_two_read_q, "fastq")
    weird_sequences_two_read_q.close()

    weird_sequences_one_read_q = open("%s/%s_Results/weird_sequences_one_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(weird_sequences_one_read, weird_sequences_one_read_q, "fastq")
    weird_sequences_one_read_q.close()

    weird_output_none_read_q = open("%s/%s_Results/weird_output_none_read_seqs.fastq" % (Data_Path,condition), "w")
    SeqIO.write(weird_sequences_none_read, weird_output_none_read_q, "fastq")
    weird_output_none_read_q.close()


#Write new SPCRs to a fasta file
new_SPCRs_q = open("%s/%s_Results/new_SPCRs_seqs.fasta" % (Data_Path,condition), "w")
SeqIO.write(new_SPCRs, new_SPCRs_q, "fasta")
new_SPCRs_q.close()

#Write just unique new SPCRs to a fasta file
for whole_seq in new_SPCRs:
    new_SPCRs_just_seqs.append (whole_seq.seq)
unique_new_SPCRs_set = set(new_SPCRs_just_seqs)
unique_spcr_number = 0
for spcr in unique_new_SPCRs_set:
    unique_new_SPCRs.append (SeqRecord(spcr, id='%s' % unique_spcr_number))
    unique_spcr_number +=1
unique_new_SPCRs_q = open("%s/%s_Results/unique_new_SPCRs_seqs.fasta" % (Data_Path,condition), "w")
SeqIO.write(unique_new_SPCRs, unique_new_SPCRs_q, "fasta")
unique_new_SPCRs_q.close()

# Find length of new SPCRs
new_SPCR_lengths = [len(seq_record) for seq_record in new_SPCRs]

number_nonwierd_reads = len(wildtype_sequences_three_read)+len(wildtype_sequences_two_read)+len(wildtype_sequences_one_read)+len(single_expansion_sequences_two_read_pos_one)+len(single_expansion_sequences_two_read_pos_two)+len(single_expansion_sequences_one_read)+len(single_replacement_sequences_two_read_pos_one)+len(double_expansion_sequences_two_read)+len(single_expansion_sequences_three_read_pos_one)+len(single_expansion_sequences_three_read_pos_two)+len(double_expansion_sequences_three_read)+len(triple_expansion_sequences_three_read)

#Write a excel file with all the relevant data
workbook = xlsxwriter.Workbook('%s/%s_Results/SPCR_analysis.xlsx' % (Data_Path,condition))
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,1,condition)
worksheet.write(1,0,'Unexpanded Reads')
worksheet.write(2,0,'Single Expansion Reads')
worksheet.write(3,0,'Double Expansion Reads')
worksheet.write(4,0,'Triple Expansion Reads')
worksheet.write(5,0,'Unexpanded (percentage)')
worksheet.write(6,0,'Single Expanded (percentage)')
worksheet.write(7,0,'Double Expanded (percentage)')
worksheet.write(8,0,'Triple Expanded (percentage)')
worksheet.write(9,0,'Percentage Expanded')
worksheet.write(11,0,'Number of New SPCRs in Position 1:')
worksheet.write(12,0,'Number of New SPCRs in Position 2:')
worksheet.write(13,0,'Number of New SPCRs in Position 1+2:')
worksheet.write(14,0,'Binned Counts (length, raw):', bold)
row = 15
col = 0
for i in range(28,41):
    worksheet.write(row,col,i)
    row +=1
#Add data
worksheet.write(1,1, len(wildtype_sequences_one_read)+len(wildtype_sequences_two_read)+len(wildtype_sequences_three_read))
worksheet.write(2,1, len(single_expansion_sequences_one_read)+len(single_expansion_sequences_two_read_pos_one)+len(single_expansion_sequences_two_read_pos_two)+len(single_replacement_sequences_two_read_pos_one)+len(single_expansion_sequences_three_read_pos_one)+len(single_expansion_sequences_three_read_pos_two))
worksheet.write(3,1, len(double_expansion_sequences_two_read)+len(double_expansion_sequences_three_read))
worksheet.write(4,1, len(triple_expansion_sequences_three_read))
worksheet.write(5,1, (float(len(wildtype_sequences_one_read)+len(wildtype_sequences_two_read))/number_nonwierd_reads)*100)
worksheet.write(6,1, (float(len(single_expansion_sequences_one_read)+len(single_expansion_sequences_two_read_pos_one)+len(single_expansion_sequences_two_read_pos_two)+len(single_replacement_sequences_two_read_pos_one)+len(single_expansion_sequences_three_read_pos_one)+len(single_expansion_sequences_three_read_pos_two))/number_nonwierd_reads)*100)
worksheet.write(7,1, (float(len(double_expansion_sequences_two_read)+len(double_expansion_sequences_three_read))/number_nonwierd_reads)*100)
worksheet.write(8,1, (float(len(triple_expansion_sequences_three_read))/number_nonwierd_reads)*100)
worksheet.write(9,1, (float(len(single_expansion_sequences_one_read)+len(single_expansion_sequences_two_read_pos_one)+len(single_expansion_sequences_two_read_pos_two)+len(single_replacement_sequences_two_read_pos_one)+len(single_expansion_sequences_three_read_pos_one)+len(single_expansion_sequences_three_read_pos_two)+len(double_expansion_sequences_two_read)+len(double_expansion_sequences_three_read)+len(triple_expansion_sequences_three_read))/number_nonwierd_reads)*100)
worksheet.write(11,1, len(single_expansion_sequences_two_read_pos_one)+len(single_expansion_sequences_three_read_pos_one))
worksheet.write(12,1, len(single_expansion_sequences_two_read_pos_two)+len(single_expansion_sequences_three_read_pos_two))
worksheet.write(13,1, len(double_expansion_sequences_two_read)+len(double_expansion_sequences_three_read))
row = 15
col = 1
for i in range(28,41):
    worksheet.write(row,col,new_SPCR_lengths.count(i))
    row +=1
workbook.close()