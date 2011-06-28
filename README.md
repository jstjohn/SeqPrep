SeqPrep is a program to merge paired end Illumina reads that are overlapping into a single longer read. It may also just be used for its adapter trimming feature without doing any paired end overlap. When an adapter sequence is present, that means that the two reads must overlap (in most cases) so they are forcefully merged. When reads do not have adapter sequence they must be treated with care when doing the merging, so a much more sensitive approach is taken. The default parameters were chosen with sensitivity in mind, so that they could be ran on libraries where very few reads are expected to overlap. It is always safest though to save the overlapping procedure for libraries where you have some prior knowledge that a significant portion of the reads will have some overlap. 


    Usage:
    ./SeqPrep [Required Args] [Options]
    NOTE 1: The output is always gziped compressed.
    NOTE 2: If the quality strings in the output contain characters less than ascii 33 on an ascii table (they look like lines from a binary file), try running again with or without the -6 option.
    Required Arguments:
    	-f <first read input fastq filename>
    	-r <second read input fastq filename>
    	-1 <first read output fastq filename>
    	-2 <second read output fastq filename>
    General Arguments (Optional):
    	-h Display this help message and exit (also works with no args) 
    	-6 Input sequence is in phred+64 rather than phred+33 format, the output will still be phred+33 
    	-q <Quality score cutoff for mismatches to be counted in overlap; default = 10>
    	-L <Minimum length of a trimmed or merged read to print it; default = 30>
    Arguments for Adapter/Primer Trimming (Optional):
    	-A <forward read primer/adapter sequence to trim as it would appear at the end of a read
    		 (should validate by grepping a file); default = AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG>
    	-B <reverse read primer/adapter sequence to trim as it would appear at the end of a read
    		 (should validate by grepping a file); default = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT>
    	-O <minimum overall base pair overlap with adapter sequence to trim; default = 10>
    	-M <maximum fraction of good quality mismatching bases for primer/adapter overlap; default = 0.130000>
    	-N <minimum fraction of matching bases for primer/adapter overlap; default = 0.700000>
    	-b <adapter alignment band-width; default = 75>
    	-Q <adapter alignment gap-open; default = 8>
    	-t <adapter alignment gap-extension; default = 2>
    	-e <adapter alignment gap-end; default = 2>
    	-Z <adapter local alignment cutoff score ((2*num_matches) - (gap_open*num_gaps) - (gap_close*num_gaps) - (gap_ext*gap_len)) ; default = 18>
    Optional Arguments for Merging:
    	-g <print overhang when adapters are present and stripped (use this if reads are different length)>
    	-s <perform merging and output the merged reads to this file>
    	-E <write pretty alignments to this file for visual Examination>
    	-x <max number of pretty alignments to write (if -E provided); default = 10000>
    	-o <minimum overall base pair overlap to merge two reads; default = 20>
    	-m <maximum fraction of good quality mismatching bases to overlap reads; default = 0.020000>
    	-n <minimum fraction of matching bases to overlap reads; default = 0.750000>


My current strategy to deal with ambiguous alignments to low complexity regions is as follows:

1. I have some minimum requirements for an overlap to be accepted
2. After the first one is found (ie the one with the maximal overlap between the two sequences), if low complexity filtering is enabled, I keep searching
3. if a second viable hit is found, I give up and say that it is not a good idea to merge the two reads.

I check for ambiguous alignments in read overlapping, but not in adapter trimming where the most conservative thing to do is strip the most aggressively aligned adapter (The closest to the beginning of the read).

To accept an alignment I allow some fraction of mismatches (currently the floor of 0.06 of the alignment length for adapter and 0.02 of the alignment length for two reads). That means that in most cases for overlapping two reads I don't allow any mismatches between adjacent reads, but if there is a 50bp potential overlap with 1 mismatch over q20 for example, I allow it. Anything below 50 needs to be perfect other than with low quality bases.

Since we ignore poor quality bases, we could have the case where a single real match followed by a long string of poor quality bases to the end of the read would result in a called overlap. That seemed like a bad idea. To get around that I require that at least some fraction of the overlapping length be matches. Right now I have that parameter set at 0.7 for adapter trimming and 0.75 for read merging, so for a case where only the last 10 bases overlap, at least 7 of those must be matches. 

Since doing that many floating point multiplications seems like a bad idea, I just have a table that pre-calculates all of those min matches and max mismatch numbers for every overlap length up to the maximum allowed read length.

Finally I have a parameter you can set which specifies a minimum resulting read length after adapter trimming and/or merging so that ultra short trimmed reads aren't output.

Following are results from hand testing the three main merge cases. Now to generate similar output automatically just supply the `-E readable_alignment.txt.gz` argument to the program (the output is gzip compressed into the file name specified).


Sequence Merge No Adapter Present:

    QUER: NCCTGCTACTACCACCCGTTCCGTGCCTGGAGCCTGCATGTTGGGCAGATACGTGCTGCCACAGCCTGTCTCTGCTGGTGCCTGGGCCTC
                                            ||  |||||||||||| || |  |||||||||||||||||||||||||||||||||
    SUBJ:                                   TGTGTGTTGGGCAGATGCGGGGGGCCACAGCCTGTCTCTGCTGGTGCCTGGGCCTCTCCTGTTCCTTGCCCACGTCTCCGTCTCCTGTTG
    RESU: NCCTGCTACTACCACCCGTTCCGTGCCTGGAGCCTGCATGTTGGGCAGATACGTGCTGCCACAGCCTGTCTCTGCTGGTGCCTGGGCCTCTCCTGTTCCTTGCCCACGTCTCCGTCTCCTGTTG
    Quality Merge:
    QUER: !223387787@@@CCC22C@@@@@@@@@@@@@@@@@@@@@@@@@@@@?@@89887:::::.2125@@:@@:::::@@@@@<<::8@@@@@
    SUBJ:                                   !!!!!!!!!!!!!!!!!!!!!!!!!!!@@@8DEGE@EDDBB2<BBE@EHBFE@EE>D8@DBE>BFIDH@IIEEIIBEIEIIGBIIGIFII
    RESU: !223387787@@@CCC22C@@@@@@@@@@@@@@@@@@@@@@@@@@@@?@@89887:::::.QPQLSSSSSSSSSSQSSSSSSSSSSSSSSD8@DBE>BFIDH@IIEEIIBEIEIIGBIIGIFII


Sequence Merge Adapter Present, Easy Peezy Mode (same lengths):

    SUBJ: NGATATGATTCCCAATCTAAGCAAACTGTCATGGAAAC
           |||||||||||||||||||||||||||||||||||||
    QUER: GGATATGATTCCCAATCTAAGCAAACTGTCATGGAAAC
    RESU: GGATATGATTCCCAATCTAAGCAAACTGTCATGGAAAC
    Quality Merge:
    SUBJ: !.-/.53444@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    QUER: IHGIIIDIIHGEHIGHIFHIFIIIIHIIIIIIIIIHII
    RESU: ISSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS


Sequence merge Adapter but lengths differ:

    SUBJ: AATTGATGGGTGCCCACCCACGGGCCAGACAAAATCATCTGGCAAGCTGGATGCAGCCTACAAGCTGTAAGATTGGA
          |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    QUER: AATTGATGGGTGCCCACCCACGGGCCAGACAAAATCATCTGGCAAGCTGGATGCAGCCTACAAGCTGTA
    RESU: AATTGATGGGTGCCCACCCACGGGCCAGACAAAATCATCTGGCAAGCTGGATGCAGCCTACAAGCTGTAAGATTGGA
    Quality Merge:
    SUBJ: =DEC??DDBD?4B=BEE@@@GB>GEE:DE8=2::6GDGBGEGDD<=;A?=AGGGG=5.=<BD?B?DDB>B4725:E>
    QUER: GDDBBFBGGFBHFIEDGGGBDGGG<GGDDG@IIIEIHDIHGIIIDDGDGDFDIFIHGIDEGGGDIIIGI
    RESU: SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSB4725:E>

