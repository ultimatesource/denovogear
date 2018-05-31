# DeNovoGear - Estimating *de novo* mutations from related individuals and cells

DeNovoGear is a software package to detect somatic and germline *de novo* mutations using next-generation sequencing data. It uses advanced statistical models to reduce the false positive rate and supports the analysis of many differential experimental designs.

**Citation:** Ramu et al. (2013) DeNovoGear: de novo indel and point mutation discovery and phasing. Nature Methods 10:985--987 doi:[10.1038/nmeth.2611](https://www.nature.com/articles/nmeth.2611)

## Table of Contents
* [Installation](#installation)
* [dng call](#dng-call-finding-mutations-in-general-pedigrees)
* [dng dnm](#dng-dnm-finding-denovo-mutations-in-trios-and-pairs)
* [dng treecall](#dng-treecall-estimating-somatic-phylogenies)
* [Release Notes](#release-notes)
* [Dependencies](#dependencies)
* [Acknowledgments](#acknowledgements)

## Installation

### Download

Source code and binaries for release versions is available at <https://github.com/denovogear/denovogear/releases>.

Source code for the most recent beta versions is available at <https://github.com/denovogear/denovogear/archive/develop.tar.gz>.

### Compiling

Compilation of DeNovoGear requires CMake. Most Linux distributions allow you to install CMake using their package software.

#### Optimized Build on Unix from a Command Line Terminal
```
    tar -xvzf denovogear*.tar.gz
    cd denovogear*/build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j4
```

### Dependencies

* Recent C++ compiler, supporting C++11 (eg. gcc 4.8.1+ or clang 3.3+)
* CMake 3.1+ when compiling <http://www.cmake.org/download/#latest>
* HTSlib 1.3+ <http://www.htslib.org/>
* Eigen 3+ <http://eigen.tuxfamily.org/>
* Boost 1.47+ <http://www.boost.org/>

Most Unix distributions contain package software that will install these dependencies for you.
DNG also can download and build missing dependencies:

```
    tar -xvzf denovogear*.tar.gz
    cd denovogear*/build
    cmake -DBUILD_EXTERNAL_PROJECTS=1 -DCMAKE_BUILD_TYPE=Release  ..
    make -j4
```

### Local Install
```
    cd denovogear*/build
    cmake -DCMAKE_INSTALL_PREFIX="${HOME}/dng" -DCMAKE_BUILD_TYPE=Release ..
    make -j4
    make install
```

### Global Install (requires root access)
```
    cd denovogear*/build
    cmake -DCMAKE_BUILD_TYPE=Release ..    
    make -j4
    sudo make install
```

### Testing

Denovogear comes with unit tests as well as full-suite test data available at <https://github.com/denovogear/testdata>. After running the build commands the tests can be downloaded and run using:
```
    make testdata
    make test  
```


## dng call: Finding Mutations in General Pedigrees

`dng call` is a module that utilizes advanced peeling algorithms to identify mutations and other statistics on zero-loop pedigrees and somatic phylogenies.
`dng call` produces valid variant call format files (vcf or bcf) that can be used with other applications.

`dng call` is capable of calculating:

 - The probability of there being exactly one, or at least one, mutation at a site.
 - The most likely genotype at a variant site, including quality score and posterior probabilities.
 - The expected number of de novo mutations.
 - The likelihood of observed data at a site.
 - The Anderson-Darling test statistic for reference vs. alternative alleles.
 - The likely gametic genotype of the parents in a trio.

### Usage

Examples:

    dng call --ped family_1.ped family_1.bam
    dng call --model=w-linked --ped family_2.ped family_2.bcf

Print Usage: `dng help call`

### Pedigree file format

`dng call` uses a 5-column, whitespace (space and tab) separated file format to specify the pedigree and somatic samples. (Note that this format is new as of June 2018.)

```
##PEDNG v1.0
#Individual Father Mother Sex Samples
NA12891 . . male =
NA12892 . . female =
NA12878 NA12891 NA12892 female (NA12878a,NA12878b);
```

Features:

 * The file must start with a `##PEDNG` header. This requirement ensures that the user is using a properly curated pedigree file.
 * Comments are lines that begin with a `#`.
 * Individual ID (first column) contains an identifier of each individual in the pedigree, along with meta information (see below).
 * Father ID and Mother ID (second and third columns) contain the identifier of the parents of each sample. The format supported is `id:length`, allowing users to scale germline mutation rates if desired. If length is not specified, it defaults to `1.0`. An id of `.` specifies that the parent is missing or invalid. Founder nodes are specified by indicating that both parents are invalid.
 * Sex (forth column) contains the information about the sex and sex chromosomes of an individuals. Both numeric and descriptive values are supported: 1=male, 2=female, and 0=autosomal.
 * Sample IDs (fifth column and beyond) specific the samples that came from the individual. Multiple samples can be specified for each individual via separate columns and/or a newick-formated tree with branch lengths. Users can use the newick tree if something is known about the somatic relationship of the samples. If a sample id is `=`, the value of the individual id will be copied. A value of `.` indicates that the sample is missing/invalid.
 * Meta information can be specified in the first column as `id@meta1@meta2@...` and can be used to specify information that affects the processing of the pedigree. Supported values:
     * Individuals default to being diploid, but this can be modified via tags: `@haploid`, `@diploid`, `@ploidy=1`, `@ploidy=2`, `@p=1`, `@p=2`.
     * Individuals can be specified as founders using the `@founder` tag. Any values in the parental columns will be ignored.
     * Individuals can be specified as clones/somatic samples of other individuals using `@clone`. Only one parent should be specified, and its ploidy will be copied to the clone.
     * Individuals can be specified as a gamete of another individual using `@gamete`. Only the appropriate parent should be specified. This tag will also set the ploidy of the individual to 1.

#### Mixed Ploidy and Mutation Accumulation Example

```
##PEDNG v1.0
Anc@diploid . . autosomal Anc-1 Anc-2 Anc-3
MA1@clone Anc:3.0 . autosomal .
MA2@clone Anc:2.0 . autosomal =
A1a@gamete MA1 . autosomal = 
A1b@gamete MA1 . autosomal =
```

#### Haplo-Diploidy Example

```
##PEDNG v1.0
drone@founder@haploid . . 1 =
queen@founder . . 2 =
bee1 drone queen 2 =
bee2@gamete . queen 1 =
```

Note that using x-linked inheritance will give you the same result as using autosomal inheritance with specified ploidies in this case

### Controlling output

`--min-quality`: threshold for reporting a mutation or variant. Based on the quality of at least one mutation.
`--all`: include sites with segregating germline variation. Based on the the quality of at least one alt allele at the site.

### Model parameters

 * Population Parameters
     - `--theta`: the population diversity parameter, i.e. 4Nu.
 * Mutation Parameters 
     - `--mu`: germline mutation rate.
     - `--mu-somatic`: somatic mutation rate. Branch lengths in each somatic phylogeny (normalized or unnormalized) are scaled by this parameter.
     - `--mu-library`: library mutation rate. Branch length connecting sequenced libraries to collected somatic samples.
 * Sequencing Parameters
     - `--lib-error`: error-rate per base-call.
     - `--lib-bias`:  reference bias in heterozygotes (ref/alt ratio).
     - `--lib-overdisp-hom`: amount of overdispersion in sequencing homozygous genotypes.
     - `--lib-overdisp-het`: amount of overdispersion in sequencing heterozygous genotypes.

To see a complete list of parameters run `dng call --help`.

### Sex-linked inheritance

`dng call` supports the analysis of data based on sex-linked chromosomal inheritance (via the `--model` flag). The supported models are `autosomal`, `x-linked`, `y-linked`, `w-linked`, `z-linked`, `mitochondrial`, and `paternal`.


## dng dnm: Finding Denovo Mutations in trios and pairs.

`dng dnm` takes in a PED file and a BCF file as input. The PED file describes the relationship between the samples and the BCF file contains the sequencing information for every locus.

### Usage

`dng dnm auto --ped testdata/sample_CEU/sample_CEU.ped --bcf testdata/sample_CEU/sample_CEU.bcf`

If you would rather start with the BAM files and have samtools installed, this would work,

`samtools mpileup -gDf hg19.fa s1.bam s2.bam s3.bam | dng dnm auto --ped sample.ped --bcf -`

### About sample.bcf

BCF files can be generated from the alignment using the samtools mpileup command. The command to generate a bcf file from sample.bam is: `samtools mpileup -gDf reference.fa sample.bam > sample.bcf` The -D option of the samtools mpileup command retains the per-sample read depth which is preferred by denovogear as it helps to filter out sites without a minimum number of reads(but note that DNG will work without per-sample RD information, in which case the RD tag encodes the average read depth information). The -g option computes genotype likelihoods and produces a compressed bcf output and the -f option is used to indicate the reference fasta file against which the alignment was built. A sample BCF file 'sample_CEU.bcf' is included in the distribution.

### About sample.ped

The PED file contains information about the trios present in the BCF file. Please make sure that all the members of the trios specified in the PED file are present in the BCF file. The PED file can be used to specify a subset of individuals for analysis from the BCF (in other words not every sample in the BCF needs to be represented in the PED file).

The PED file is a tab delimited file. The first six columns of the PED file are mandatory, these are Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1 = male; 2 = female; other = unknown) and Phenotype. `dng dnm` makes use of the first four columns. The sample IDs in the PED file need to be exactly as they appear in the BCF file header. Sample order within the PED file does not matter, as family relationships are completely specified by the value of the child/mother/father fields in each row.

For example, a single line in the PED file that specifies a trio looks like:

```
CEU	NA12878	NA12891	NA12892	2	2
```

### About "snp_lookup.txt" and "indel_lookup.txt"

These are tables with precomputed priors (and other useful numbers) for all possible trio configurations, under the null (no mutation present) and alternative (true de novo). The default tables are generated during each program run using a prior of 1 x 10 ^-8 /bp/generation on the haploid germline point mutation rate, and 1 x 10 ^-9 /bp/generation on the haploid germline indel mutation rate.

If you wish to change the default point or indel mutation rates use the --snp_mrate or --indel_mrate switches respectively.

For example, `dng dnm auto --ped testdata/sample_CEU/sample.ped --bcf testdata/sample_CEU/sample.bcf --snp_mrate 2e-10 --indel_mrate 1e-11`

The indel mutation rate prior is calculated based on the length of the insertion or deletion event, separate models are used for insertions and deletions. The two models are based on the indel observations from the 1000Genomes phase 1 data.

The insertion mutation rate is modeled using the function
 	log (mrate) = mu_scale * (-22.8689 - (0.2994 * insertionLength))

The deletion mutation rate is modeled using the function
	log (mrate) = mu_scale * (-21.9313 - (0.2856 * deletionLength))

Note that a constant factor is used to scale the mutation rate, it is set to 1.0 by default and can be set using the switch `--mu_scale`. This provides the users to scale the mutation rate prior according to their data-set.

For example, `dng dnm auto --ped testdata/sample_CEU/sample.ped --bcf testdata/sample_CEU/sample.bcf --mu_scale 3`


### Output format for trios

The output format is a single row for each putative de novo mutation (DNM), with the following fields

    1. Event type (POINT MUTATION or INDEL).
    2. CHILD_ID - sample ID of trio-offspring with the DNM.
    3. chr - chromosome.
    4. pos - physical Position.
    5. ref - base present in reference sequence at this position (from BCF).
    6. alt - comma separated list of alternate non-reference alleles called on at-least one sample (from BCF).
    7. maxlike_null  - likelihood of the most likely mendelian-compatible genotype configuration.
    8. pp_null - posterior probability of the most likely mendelian-compatible genotype configuration.
    9. tgt_null(child/mom/dad)  - genotypes of the most likely mendelian-compatible configuration.
    10. snpcode - code that indicates whether the configuration shown in field 6 is monomorphic (1) or contains variation (2)(internal filters, can be ignored).
    11. code - This field seems to be redundant to field 7, except the codes are (6) and (9).(internal filters, can be ignored).
    12. maxlike_dnm - likelihood of the most likely DENOVO genotype configuration.
    13. pp_dnm - posterior probability of the most likely DENOVO genotype configuration. Range = [0-1], a value closer to 1 indicates higher probability of observing a denovo event at this position. This is the field that is used to rank the calls.
    14. tgt_dnm(child/mom/dad)  - genotypes of the most likely mendelian-compatible configuration.
    15. lookup - Code that indicates if the most likely DNM is a transition (4) or transversion (5) (for development use only).
    16. flag - Flag that indicates whether the data for the site passed internal QC thresholds (for development use only).
    17-19. Read depth of child, parent 1 and parent 2.
    20-22. Root mean square of the mapping qualities of all reads mapping to the site for child, parent 1 and parent 2. Currently these values are the same for all samples when using BCF as the input format.

Fields 17--22 are meant for filtering out low quality sites.

### Separate models for the X chromosome

`dng dnm` has separate models for autosomes, X chromosome in male offspring and X chromosome in female offspring. To use this model, create separate BCFs for the X chromosome.

#### Autosomes model usage

`dng dnm auto --ped sample.ped --bcf sample.bcf`

#### X in male offspring model usage

`dng dnm XS --ped sample.ped --bcf sample.bcf --region X`

#### X in female offspring model usage

`dng dnm XD --ped sample.ped --bcf sample.bcf --region X`

### Paired-Sample Analysis

DNG can be used to analyze paired samples for example to call somatic mutations between tumor/matched-normal pairs, the main difference in how to run the program is the way samples are specified in the PED file(see below),

#### Usage

`dng dnm auto --ped paired.ped --bcf sample.bcf`

About the arguments,

`paired.ped` is a ped file containing the family-name and the name of the two samples in the bcf file. The last three columns are mandated by the 	PED format but are ignored by the program. An example line looks like

```
YRI	NA19240_blood	NA19240_vald	0	0	0
```

`sample.bcf` is a BCF/VCF file containing both the samples.

#### Output format for paired sample analysis

The output format is a single row for each putative paired denovo mutation(DNM), with the following fields

    1. Event type (POINT MUTATION or INDEL).
    2. TUMOR_ID - sample ID of the 'TUMOR' sample.
    3. NORMAL_ID - sample ID of the 'NORMAL' sample.
    4. chr - chromosome.
    5. pos - physical Position.
    6. ref - base present in reference sequence at this position (from BCF).
    7. alt - comma separated list of alternate non-reference alleles called on at-least one sample (from BCF).
    8. maxlike_null  - likelihood of the most likely compatible genotype configuration.
    9. pp_null - posterior probability of the most likely compatible genotype configuration.
    10. tgt_null(normal/tumor)  - genotypes of the most likely compatible configuration.
    11. maxlike_dnm  - likelihood of the most likely denovo genotype configuration.
    12. pp_dnm - posterior probability of the most likely denovo genotype configuration, a value closer to 1 indicates strong evidence of a denovo event. This is the field that is used to rank the calls.
    13. tgt_dnm(normal/tumor)  - genotypes of the most likely denovo configuration.
    14-15. Read depth of tumor, normal samples.
    16-17. Root mean square of the mapping qualities of all reads mapping to the site for tumor, normal samples.
    18-19. null_snpcode, dnm_snpcode - snpcode is a field used to classify the genotype configurations for the null and alternate case. 1 stands for hom/hom, 2 stands for het/hom, 3 stands for hom/het and 4 stands for het/het.


## dng phaser: read-based phasing if de novo mutations

`dng phaser` can be used to obtain parental phasing information for de novo Mutations where phase informative sites are present. This is done by looking at reads which cover both the de novo base and a phase informative positions. Phase informative positions are SNP positions that lie within a certain window from the de novo site, the default window size is 1000 bp but the window-size can be set by the user.

For each DNM, a list of phase informative sites from the genotypes file is obtained, these are the loci which lie within the phasing window. Also, these sites do not have a het/het GT configuration for the parents and where the child is het. The number of DNM and phasing  allele combinations seen in the read level data is output by the program. For the phasing sites, the inferred parental origin is displayed.

For example if the base at the phasing site is T and the parental genotypes are P1:CC and p2:TC at this site, then the parent of origin of this base is p2. By looking at the base in the de novo position on this read it is possible to infer the parent of origin of the denovo mutation.

### Usage

`dng phaser --dnm dnm_f --pgt gts_file --bam alignment.bam --window 1000`

About the arguments,

 1. "dnm" is the list of DNMs whose parental origin is to be determined. It is a tab delimited file of the format
  	chr pos inherited_base mutant_base
 2. "pgt" contains the genotypes of the child and the parents at SNP sites. It is a tab delimited file of the format
	chr pos child_GT parent1_GT parent2_GT
 3. "bam" is the alignment file (.bam) of the child.
 4. "window" is an optional argument which is the maximum distance between the DNM and a phasing site. The default value is 1000.

### Output

	DNM_pos 1:182974758     INHERITED G     VARIANT A
        HAP POS 182974328 p1: CC p2: TC Number of denovo-phasing pairs found: 0
        HAP POS 182974572 p1: CC p2: TC Number of denovo-phasing pairs found: 2
                Base at DNM position: A Base at phasing position: C      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 3
                Base at DNM position: G Base at phasing position: T      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 3
        HAP POS 182974598 p1: GG p2: CG Number of denovo-phasing pairs found: 2
                Base at DNM position: A Base at phasing position: G      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 1
                Base at DNM position: G Base at phasing position: C      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 2
        HAP POS 182974602 p1: TT p2: CT Number of denovo-phasing pairs found: 2
                Base at DNM position: A Base at phasing position: T      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 1
                Base at DNM position: G Base at phasing position: C      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 1
        HAP POS 182974707 p1: TT p2: CT Number of denovo-phasing pairs found: 2
                Base at DNM position: A Base at phasing position: T      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 25
                Base at DNM position: G Base at phasing position: C      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 24
        HAP POS 182974750 p1: AA p2: TA Number of denovo-phasing pairs found: 2
                Base at DNM position: A Base at phasing position: A      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 45
                Base at DNM position: G Base at phasing position: T      INFERRED PARENT OF ORIGIN for DNM: p1 SUPPORTING READ COUNT: 47

Please feel free to contact the authors about any concerns/comments.

### General options for trios and paired sample calling

`--snp_mrate`:     Mutation rate prior for SNPs. [1e-8]
`--indel_mrate`:   Mutation rate prior for INDELs. [1e-9]
`--pair_mrate`:    Mutation rate prior for paired sample analysis. [1e-9]
`--indel_mu_scale`:        Scaling factor for indel mutation rate. [1]
`--pp_cutoff`:     Posterior probability threshold. [0.0001]
`--rd_cutoff`:     Read depth filter, sites where either one of the sample have read depth less than this threshold are filtered out. [10
`--region`: region of the BCF file over which to perform denovo calling. [string of the form "chr:start-end"]


## dng treecall: Estimating somatic phylogenies

`dng treecall` is a program to estimate somatic phylogenies from VCF/BCF files with annotated PL fields.

### Usage Example

```{bash}
# Search treespace
dng treecall search samples.vcf output

# Genotype samples based on best scoring tree
dng treecall genotype -t output.best.tree samples.vcf output.gtcall

# Annotate tree with number of mutations that occur on each branch
dng treecall annotate -t output.best.tree samples.vcf output.gtcall output.best.annotated.tree
```

### dng treecall: Finding Mutations and Trees from Single Cell Sequencing

`dng treecall` is an experimental module to identify mutations and trees from whole-genome sequencing of single cells

#### Usage:

Print Usage: `dng treecall -h`

see https://github.com/nh3/treecall/

## ACKNOWLEDGEMENTS
DeNovoGear uses the Samtools libraries to parse BCF files, we thank the Heng Li and the rest of the authors of Samtools for making this resource available to developers.
