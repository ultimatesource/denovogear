Columns:

1-SNP status: total number of unique alleles observed at the site across all three samples in trio. Not used in calculations or algorithm

2-Informative flag: flag that categorizes site based on trio genotype calls, not used for calculations, but it is used in the algorithm.

3-Transmission probability, corresponds to L(G_c|G_m,G_f), used in calculations.
 -values: if snp_status >=4 then tprob=0.
          if snp_status == 3, then tprob is assigned using the following rules:
          -if both child alleles are in parents, then assign standard mendelian
           tprobs (1/4 or 1/2). 
          -if site has a DNM, then set tprob=0.
          if snp_status == 2, then use these rules:
          -if mendelian, then set standard tprobs
          -if de novo, _AND_ trio configuration fits one of these models for CC/MM/DD
:           AA/BB/AA, AA/AA/BB then set tprob=0. 

6-Denovo flag: should roughly be 0 if there is no denovo, 1 if there is. This is used in the DNG executable for picking which matrix entries correspond to a de novo state 

7-Normal flag: analogous to column 6, but picking out normal mendelian sites instead.
