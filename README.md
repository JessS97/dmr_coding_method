## Software: DNA Data Storage: Error Correction with Dynamic Mapping Rule

For more information on the purpose of the programme, please refer to the paper.md file. 

### Installation: 

1. Install anaconda
   - Download anaconda: https://www.anaconda.com/data-science-platform

2. Create the environment with anaconda
   - Use this command in your terminal: conda env create -f environment.yml
   
### Programm structure
The following section describes the programme files and their functions. 

#### run.py
This option parser can be used to run the programm which will do the encoding, error simulation and decoding based on the user inputs.
Example to run the file: 
python run.py  -o "output//" -r 32 -c 1 -s 10,15,20,25,30,35,40,45,50 -l short -e "DMR" -f "subs" -t 10 -p 0.01,0.02

#### Input images
The fraunhofer logo was used for the encoding. The logo was additionally reduced to 47x47 pixels to shorten the time for testing. 
Example images:
Logo_long__binary_.png
Logo_short__binary_.png

#### Mapping Tables
Mapping tables to translate binary strings into DNA and back. 
mapping_table_dmr.json
mapping_table_two_bit.json

#### main_functions.py
This Python script contains the encoding and decoding functions for the encoding methods RS, RS bytearray segmentation, RS bitarray segmentation, DMR and DMR 
bitarray segmentation.  
##### Single RS fragment
With this method, the RS code is applied to a long segment. The image is translated into a complete byte sequence and encoded with RS. When coding with RS, the
Galois field must be adapted to the length of the resulting overall sequence of data and RS symbols. A Galois field with a size of 8 is selected as standard,
which means that a total segment size of 2^8 - 1 = 255 can be used for encoding. If RS(255,223) is used, the segment consists of 223 data symbols and 32 RS 
symbols. If a larger Galois field of e.g. 12 is used, a total of 4095 (2^12 - 1 = 4095) symbols can be used for encoding the data and RS symbols.
##### RS/DMR with byte array segmentation
In the second method, the RS code is applied to individual subsections, referred to here as segments. The compressed binary string is translated into bytes 
and divided into segments of a certain size, to which the RS code is applied separately. The segments can later be separated from each other by adding spacer 
sequences. Depending on how large the spacer sequences are selected for later use, the length of the DNA strands to be stored varies. The advantage of this
method is that the spacer sequences are located at a certain distance from each other and should therefore be easy to find.
##### RS/DMR with bitarray segmentation
In the third encoding scheme RS with Bitarray segmentation the input data was first converted into bits as described below. These Bitarray was separated into
segments with a specific length so that a fixed number of bits or pixels could be assigned to the segments (Figure 4 c). Each segment got compressed using 
PackBits, which results in compressed segments with different lengths. Due to the different lengths of the segments, they were coded with RS symbols on the 
basis of proportions calculated from RS(255,223). The merging of the segments using spacers results in spacer located at different distances on the resulting 
segment.

#### functions.py
Python script containing functions used in main_functions.py. These functions are for reading and writing images, bit to byte translations, insertion and 
removing of a random mask, segmentation, encoding and decoding with the DMR scheme, DNA-base checks of the DMR scheme, encoding and decoding with RS, 
generation if random errors and the calculation of the edit distance. 

#### dmr_rs_coder.py 
This file contains the clas of the DMR encoder/decoder, which is designed similar to the Reed-Solomon encoder/decoder.
The class contains methods to encode binary to DNA and decode DNA to binary. The used DMR (Dynamic Mapping Rule) mapping-rule is an attempt to utilize 
a self-correcting DNA code. This mapping rule is combined with the Reed-Solomon code to further enhance the correctability of the sequence. Reed-Solomon and
DMR interact with each other so that in synergy a better correction should be possible.

#### dmr_level_master.py
During correction, the sequences are first checked to see whether they are already correct or can already be corrected using the attached RS characters. 
Sorting out the sequences serves to speed up the correction and only carry out further correction on segments where this is necessary. The remaining sequences 
are corrected by running them through different levels, which have different assumptions and test correction options. The first step is to check whether the 
bases of the segment fit into the mapping scheme or whether conspicuous areas (inconsistencies/inconsistencies) are already recognisable here. A total of 3 
levels were designed to correct substitution errors.

#### output folder
This folder saves the results of the coding and decoding. This includes log files, the edit distances, the decoded pictures, the decoded bit sequences and the 
mutated dna sequences.
