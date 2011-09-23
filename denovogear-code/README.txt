#### Rough Draft Instructions for compiling denovogear #########
Authors: Don Conrad & Avinash Ramu
09/23/11


***** INSTALLATION ******
Requires newmat11 library, http://www.robertnz.net/nm11.htm

cd to newmat, if you have gnu compiler you can compile with
'make -f nm_gnu.mak'

after compiling newmat, change back to this directory and 
compile with "make". If you chose to keep the newmat 
library elsewhere, be sure to edit the denovogear Makefile 
NEWMAT= definition accordingly. 

***** RUNNING THE CODE *****

Input requires a PED file, a BCF file and 2 lookup files called "lookup.txt" and "lookup_indel.txt"

usage:
	./denovogear sample.ped sample.bcf

about sample.bcf:
BCF files can be generated from the alignment files using the samtools mpileup command. The -D option of the mpileup command retains the per-sample read depth which is preferred by denovogear (but note that DNG will work without per-sample RD information). 

For example to generate a bcf file from sample.bam the command is:
	samtools mpileup -gDf reference.fa sample.bam > sample.bcf

Note that the -g option specifies a compressed output and the -f option is used to indicate the reference.

about sample.ped:
The PED file contains information about the trios present in the BCF file. Please make sure that all the members of the trios specified in the PED file are present in the BCF file. The PED file can be used to specify a subset of individuals for analysis from the BCF (that is not every sample in the BCF need be represented in the PED file). An example PED file, CEU.ped, is included in the distribution directory. 

about "lookup.txt" and "lookup_indel.txt":
This is a table with precomputed priors (and other useful numbers) for all possible 
trio configurations, under the null (no mutation present) and alternative (true de novo). 
This file needs to be in the working directory. The default tables were generated using a prior of 1 x 10 ^-8 /bp/generation on the haploid germline point mutation rate, and 1 x 10 ^-9 /bp/generation on the haploid germline indel mutation rate. 

New versions of the table, for instance using a different prior on the mutation rate, can be generated using the enclosed utility program "make_lookup.pl". Run the script without arguments to see usage. A similar lookup table for indels is stored in the "lookup_indel.txt" file.


***** OUTPUT FORMAT *******

The output format is a single row for each putative de novo mutation (DNM), with the following fields:
1. Event type (POINT MUTATION or INDEL)
2. Sample ID of offspring with the DNM
3. Chromosome 
4. Physical Position 
5. Base present in reference sequence at this position
6. maxlike_null  - likelihood of the most likely mendelian-compatible config.
7. pp_null - posterior probability of most likely mendelian configuration        
8. tgt  - genotypes of the most likely mendelian configuration
9. Code that indicates whether the configuration shown in field 6 is monomorphic (1) or contains variation (2)
10. This field seems to be redundant to field 7, except the codes are (6) and (9).

11. maxlike_DNM  -11, 12 and 13 are analogous to 6,7,8, but for a de novo mutation
12. posterior_probability_DNM
13. tgt: DNM_configuration

14. Code that indicates if the most likely DNM is a transition (4) or transversion (5)
15. This is a flag that indicates wether the data for the site passed internal QC thresholds (for development use only).

16-18. Read depth of child, parent 1 and parent 2. 
19-21. Root mean square of the mapping qualities of all reads mapping to the site for child, parent 1 and parent2. Currently these values are the same for all samples when using BCF as the input format.

Fields 16-21 are meant for filtering out low quality sites. 


############
