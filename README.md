# crispr-images
Code for Shipman SL, Nivala J, Macklis JD, Church GM. (2017). **CRISPR-Cas encoding of a digital movie into the genomes of a population of living bacteria**. *Nature*.

*These Python scripts generate nucleotide sequences in a protospacer format for encoding images into live cells using CRISPR adaptation, analyze the spacer content of bacteria to identify newly acquired spacers, and reconstruct images that have been encoded in newly acquired spacers.*

# Contents of this repository:
* **seq_from4color_56x56.py**
  
  Generates protospacer sequences to encode a 4 color, 56x56 pixel image  
  *Input*: .tif image (see Image_set_hand4c)  
  *Output*: Excel file with protospacer sequences

* **seq_from21color_30x30.py**
  
  Generates protospacer sequences to encode a 21 color, 30x30 pixel image  
  *Input*: .tif image (see Image_set_hand21c)  
  *Output*: Excel file with protospacer sequences

* **seq_from21colorMovie_36x26.py**
  
  Generates protospacer sequences to encode a 21 color, 36x26 pixel GIF of five frames  
  *Input*: five .tif images (see Image_set_horse_movie_21c_36_26)  
  *Output*: Excel file with protospacer sequences
  
* **Extract_Spacers.py**
  
  Identifies newly acquired spacers from a population of bacteria  
  *Input*: fastq file of sequencing reads containing the CRISPR array  
  *Output*: 1. fasta file of newly acquired sequences  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2. fastq files of sequenced binned by how many new spacers are contained in the array  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3. excel file describing acquisition events (percentage of arrays expanded, length of new spacers, etc.)

* **SPCR_blast.py**

  Identifies the origin of newly acquired spacers from a population of bacteria: genome, plasmid, or unaligned  
  *Input*: fasta file of newly acquired sequences (additionally requires local blast database for genome and plasmid)  
  *Output*: 1. fasta files of sequences by origin with surrounding nucleotide context (for genome/plasmid)  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2. excel file with spacer numbers by origin and a list of the most commonly acquired unaligned spacers

* **Analyze_Image_4c_56x56.py**

  Reconstructs a 4 color, 56x56 pixel image from newly acquired spacers  
  *Input*: fasta file of newly acquired sequences (additionally requires local blast database for genome and plasmid)  
  *Output*: .png image

* **Analyze_Image_21c_30x30.py**

  Reconstructs a 21 color, 30x30 pixel image from newly acquired spacers  
  *Input*: fasta file of newly acquired sequences (additionally requires local blast database for genome and plasmid)  
  *Output*: .png image

* **Analyze_GIF.py**

  Reconstructs a 21 color, 36x26 pixel, 5 frame GIF from newly acquired spacers  
  *Input*: fasta file of newly acquired sequences (additionally requires local blast database for genome and plasmid)  
  *Output*: 1. five .png images (by each frame)  
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2. excel file comparing the recovered information to the intended information 
