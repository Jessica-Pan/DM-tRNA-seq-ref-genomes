# DM-tRNA-seq-ref-genomes
Jessica Pan, 8/8/2019

In order to do the run described in the Word document tRNA-seq outputs 080719, first run "Reference Sequence Alignment.py" then "creating output from sequence_alignments.py."

Description of included files:
-hg19CCA+mito.fasta is a list of the seed sequences used for the initial run of this code. It may or may not still be relevant to you. The seed sequence fileName is specified in "Reference Sequence Alignment.py" if you would like to change it. The raw data files are also specified in this file.

-Reference Sequence Alignment.py is the main code for this run. It easily takes up the bulk of the runtime. This file aligns sequences in the raw data to the seed sequences, with each sequence aligning to only one seed. It creates a folder called sequence_alignments, where it creates a fasta file for each seed sequence, containing its mapped reads. The three variables easily accessible at the top of the file are the paths to the seed sequence list(s), the paths to the raw data files, and the toggle boolean that determines whether or not seed sequences without reads output empty files for completion.

-creating output from sequence_alignments.py takes the folder sequence_alignments (created by the file Reference Sequence Alignments.py) and creates all of the outputs described in the word document that describes the final outputs of this program. These go into the folder tRNA-seq-outputs. These include: a file detailing the data quality in each raw data file, a breakdown of the data by isodecoder, a breakdown of the data by anticodon, a breakdown of the charging proportion by isodecoder, and an output for each sequence fulfilling certain criteria that describes mutation frequency by position in the seed sequence and various other stats. 

-reformatting dataset.py attempts to reformat your seed sequence database. This database must be in fasta format, preferably with the type of header described in the word document. If the header type differs, the code will likely be unable to identity a proper isodecoder and will use the entire ID line instead. 

-Within the folder "Alternate Version," lives a slightly older version of the code that was repurposed for the ALKb dataset. The main difference is the addition of running each seed sequence through the code twice, the second time containing a 20 nucleotide sequence that could appear in these sequences before the acceptor CCA. 
-The file output from sequence_alignments alkb.py is very similar to its alternate version outside of this folder. It may vary slightly. 