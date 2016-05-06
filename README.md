SeqPrep is a program to merge paired end Illumina reads that are overlapping into a single longer read. It may also just be used for its adapter trimming feature without doing any paired end overlap. When an adapter sequence is present, that means that the two reads must overlap (in most cases) so they are forcefully merged. When reads do not have adapter sequence they must be treated with care when doing the merging, so a much more specific approach is taken. The default parameters were chosen with specificity in mind, so that they could be ran on libraries where very few reads are expected to overlap. It is always safest though to save the overlapping procedure for libraries where you have some prior knowledge that a significant portion of the reads will have some overlap. 

Before running SeqPrep make sure to check that the program's defaults are indeed the adapters you are looking for. Try copying the default forward adapter from this file and grep it against your reads doing a word count, also try the same with the reverse adapter with grep. You should see some hits. You can also try using (and validating with grep) `-A GATCGGAAGAGCACACG -B AGATCGGAAGAGCGTCGT` as parameters. To find a list of Illumina adapter sequences you should write to Illumina tech support TechSupport@illumina.com (they do not like people to share the list of sequences outside of their institution).

Chose about 20bp of an adapter sequence where:

1.    You see the most hits with grep
2.    When you run a command like `zcat Lane2_0d_2.fastq.gz | head -n 1000000 |grep "INSERT ADAPTER HERE" | head` you see the adapter sequence show up at the beginning of a few reads. Also the -A and -B arguments should be as they show up in your data, SeqPrep searches directly for these sequences without doing reverse complementing.
3.    Check the forward and reverse and make sure that you have roughly the same number of hits via a command to count hits like: `zcat Lane2_0d_2.fastq.gz | head -n 1000000 |grep "INSERT ADAPTER HERE" | wc -l`

As an additional precaution, the program checks for good read overlap once the adapters are trimmed. If the adapter is trimmed and the reads do not have a reasonable adapter overlap (you can modify this setting with -X) then the reads aren't printed or merged. 

See `Test/README.md` for some information on testing out other parameters. `Test/SimTest` has some particularly cool test data which you can use to check out sensitivity and specificity of adapter trimming using different parameters. The results of the test are displayed in `results.html` which uses the google charts API so that the points are interactive and you can easily determine which settings made which points.


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

	-3 <first read discarded fastq filename>
	-4 <second read discarded fastq filename>
	-h Display this help message and exit (also works with no args) 
	-6 Input sequence is in phred+64 rather than phred+33 format, the output will still be phred+33 
	-q <Quality score cutoff for mismatches to be counted in overlap; default = 13>
	-L <Minimum length of a trimmed or merged read to print it; default = 30>

Arguments for Adapter/Primer Trimming (Optional):

	-A <forward read primer/adapter sequence to trim as it would appear at the end of a read (recommend about 20bp of this)
		 (should validate by grepping a file); default (genomic non-multiplexed adapter1) = AGATCGGAAGAGCGGTTCAG>
	-B <reverse read primer/adapter sequence to trim as it would appear at the end of a read (recommend about 20bp of this)
		 (should validate by grepping a file); default (genomic non-multiplexed adapter2) = AGATCGGAAGAGCGTCGTGT>
	-O <minimum overall base pair overlap with adapter sequence to trim; default = 10>
	-M <maximum fraction of good quality mismatching bases for primer/adapter overlap; default = 0.020000>
	-N <minimum fraction of matching bases for primer/adapter overlap; default = 0.870000>
	-b <adapter alignment band-width; default = 50>
	-Q <adapter alignment gap-open; default = 8>
	-t <adapter alignment gap-extension; default = 2>
	-e <adapter alignment gap-end; default = 2>
	-Z <adapter alignment minimum local alignment score cutoff [roughly (2*num_hits) - (num_gaps*gap_open) - (num_gaps*gap_close) - (gap_len*gap_extend) - (2*num_mismatches)]; default = 26>
	-w <read alignment band-width; default = 50>
	-W <read alignment gap-open; default = 26>
	-p <read alignment gap-extension; default = 9>
	-P <read alignment gap-end; default = 5>
	-X <read alignment maximum fraction gap cutoff; default = 0.125000>
	-z <use mask; N will replace adapters>

Optional Arguments for Merging:

	-y <maximum quality score in output ((phred 33) default = ']' )>
	-g <print overhang when adapters are present and stripped (use this if reads are different length)>
	-s <perform merging and output the merged reads to this file>
	-E <write pretty alignments to this file for visual Examination>
	-x <max number of pretty alignments to write (if -E provided); default = 10000>
	-o <minimum overall base pair overlap to merge two reads; default = 15>
	-m <maximum fraction of good quality mismatching bases to overlap reads; default = 0.020000>
	-n <minimum fraction of matching bases to overlap reads; default = 0.900000>



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


If interested there is a website where I post my tests of different parameters for SeqPrep on simulated data. There are also a few comparison stats of different programs to trim adapters. The website can be accessed here:
http://hgwdev.cse.ucsc.edu/~jstjohn/seqprep/

where the pages are named `result(date).html`. The latest ones (as of when I have gotten around to edit this) can be found here:

http://hgwdev.cse.ucsc.edu/~jstjohn/seqprep/results2011-09-15.html

Note that although my program is more sensitive and specific than fastq-clipper, I optimized my default parameters based on this test. Results on real data may be
different, although I believe my method takes advantage of a more realistic adapter model than other software does. For example, even though my program requires
10bp of adapter to be present at the end of a read to trim it off (by default) there is a backup adapter trimming function that trimms based on strong and
unambiguous read overlap. Because of this my program can trim the adapter even if it is only present in the last few bases of the read.

Also note that fastq-mcf appears to do a little better at sensitivity (0.992 vs 0.985) at a very large cost to specificity (0.497 vs 0.994).
