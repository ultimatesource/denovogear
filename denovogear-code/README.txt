#### Rough Draft Instructions for compiling denovogear #########
Author: Don Conrad
06/12/10


***** INSTALLATION ******
Requires newmat11 library, http://www.robertnz.net/nm11.htm

cd to newmat, if you have gnu compiler you can compile with
'make -f nm_gnu.mak'

after compiling newmat, change back to this directory and 
compile with "make". If you chose to keep the newmat 
library elsewhere, be sure to edit the denovogear Makefile 
NEWMAT= definition accordingly.

***** RUNNING THE CODE *****

Input requires three GLF files and a lookup file called "lookup.txt".


usage is:
./denovogear child.glf parent1.glf parent2.glf > outfile

about "lookup.txt":
This is a table with precomputed priors (and other useful numbers) for all possible 
trio configurations, under the null (no mutation present) and alternative (true de novo). 
This file needs to be in the working directory. The default table was generated using a prior 
of 1 x 10 ^-8 /bp/generation on the haploid germline mutation rate. 

new versions of the table, for instance using a different prior on the mutation rate, using 
the enclosed utility program "make_lookup.pl"


***** OUTPUT FORMAT *******

The output format is a single row for each putative de novo mutation (DNM), with the following fields:

1.Chromosome 
2. Position 
3. Reference_base  base present in ucsc18/NCBI36 at this position

4. likelihood_mendelian  - likelihood of the most likely mendelian-compatible config.
5. posterior_probability_mendelian  posterior probability         
6. mendelian_configuration  - genotypes of the most likely mendelian configuration

7. Code that indicates whether the configuration shown in field 6 is monomorphic (1) or contains variation (2)
8. I believe this field is completely redundant to field 7, except the codes are (6) and (9).

9. likelihood_DNM  - 7, 8 and 9 are analogous to 4,5,6, but for a de novo mutation
10. posterior_probability_DNM
11. DNM_configuration

12. Code that indicates if the most likely DNM is a transition (4) or transversion (5)
13. I dont know what this code is, right now!

15-17. Read depth of child, parent 1 and parent 2. 
19-20. Root mean square of the mapping qualities of all reads mapping to the site for child, parent 1 and parent2. 

Fields 15-20 are meant for filtering out low quality sites. 

Example:

chr1 16347710 C 9.53194e-10 0.131685 CT/CT/CC 2 9 9.98e-09 0.868315 CT/CC/CC 4 0 depth 14       21      17       Qs 60  60      60

############
