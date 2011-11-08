#### Rough Draft Instructions for compiling denovogear #########
Authors: Don Conrad & Avinash Ramu
11/2/11


***** INSTALLATION ******
Requires newmat11 library, http://www.robertnz.net/nm11.htm

cd to newmat, if you have gnu compiler you can compile with
'make -f nm_gnu.mak'

after compiling newmat, change back to this directory and 
compile with "make". If you chose to keep the newmat 
library elsewhere, be sure to edit the denovogear Makefile 
NEWMAT= definition accordingly. 

***** RUNNING THE CODE *****

Input requires a PED file and a BCF file.

usage:
	./denovogear --ped sample.ped --bcf sample.bcf

about sample.bcf:
BCF files can be generated from the alignment files using the samtools mpileup 
command.  

For example the command to generate a bcf file from sample.bam is:
	samtools mpileup -gDf reference.fa sample.bam > sample.bcf

The -D option of the samtools mpileup command retains the per-sample read depth 
which is preferred by denovogear (but note that DNG will work without per-sample 
RD information). The -g option specifies a compressed output and the -f option 
is used to indicate the reference. 

about sample.ped:
The PED file contains information about the trios present in the BCF file. 
Please make sure that all the members of the trios specified in the PED file 
are present in the BCF file. The PED file can be used to specify a subset of 
individuals for analysis from the BCF (that is not every sample in the BCF need 
be represented in the PED file).

The PED file is a tab delimited file. The first six columns of the PED file are 
mandatory, these are Family ID, Individual ID, Paternal ID, Maternal ID, 
Sex (1 = male; 2 = female; other = unknown) and Phenotype. The sample ID's in 
the PED file need to be exactly the same as they appear in the BCF file header. 
Sample order within the PED file does not matter, as family relationships are 
completely specified by the value of the child/mother/father fields in each row.
 
For example, a single line in the PED file that specifies a trio looks like:
CEU	NA12878_vald-sorted.bam.bam	NA12891_vald-sorted.bam.bam	NA12892_vald-sorted.bam.bam	2	2

An example PED file, CEU.ped, is included in the distribution directory. 

about "snp_lookup.txt" and "indel_lookup.txt":
These are tables with precomputed priors (and other useful numbers) for all possible 
trio configurations, under the null (no mutation present) and alternative (true de novo). 
The default tables are generated during each program run using a prior of 
1 x 10 ^-8 /bp/generation on the haploid germline point mutation rate, 
and 1 x 10 ^-9 /bp/generation on the haploid germline indel mutation rate. 

If you wish to change the default point or indel mutation rates use the --snp_mrate 
or --indel_mrate switches respectively. 

For example
	./denovogear --ped sample.ped --bcf sample.bcf --snp_mrate 2e-10 --indel_mrate 1e-11

The indel mutation rate varies according to the length of the insertion or deletion, 
separate models are used for insertions and deletions. The two models were calibrated
based on the indel observations from the 1000Genomes phase 1 data.

The insertion mutation rate is modeled using the function
 	log (mrate) = mu_scale * (-22.8689 - (0.2994 * insertionLength))

The deletion mutation rate is modeled using the function
	log (mrate) = mu_scale * (-21.9313 - (0.2856 * deletionLength))

Note that a constant factor is used to scale the mutation rate, it is set to 1.0 
by default and can be set using the switch --mu_scale. 

For example, 
	./denovogear --ped sample.ped --bcf sample.bcf --mu_scale 3


***** OUTPUT FORMAT *******

The output format is a single row for each putative de novo mutation (DNM), with the following fields:
1. Event type (POINT MUTATION or INDEL)
2. Sample ID of offspring with the DNM
3. Chromosome 
4. Physical Position 
5. Base present in reference sequence at this position
6. ALT - Comma separated list of alternate non-reference alleles called on at-least one sample.
7. maxlike_null  - likelihood of the most likely mendelian-compatible config.
8. pp_null - posterior probability of most likely mendelian configuration        
9. tgt  - genotypes of the most likely mendelian configuration
10. Code that indicates whether the configuration shown in field 6 is monomorphic (1) or contains variation (2)
11. This field seems to be redundant to field 7, except the codes are (6) and (9).

12. maxlike_DNM  -11, 12 and 13 are analogous to 6,7,8, but for a de novo mutation
13. posterior_probability_DNM
14. tgt: DNM_configuration

15. Code that indicates if the most likely DNM is a transition (4) or transversion (5)
16. This is a flag that indicates whether the data for the site passed internal QC thresholds (for development use only).

17-19. Read depth of child, parent 1 and parent 2. 
20-22. Root mean square of the mapping qualities of all reads mapping to the site for child, parent 1 and parent 2. Currently these values are the same for all samples when using BCF as the input format.

Fields 17-22 are meant for filtering out low quality sites. 


############
