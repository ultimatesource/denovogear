## Release Notes

### develop

* Lots of new features and bug fixes are added to the develop branch, which is generally stable. Consider trying it out.

#### latest changes in develop

* CHANGE: HTSLIB 1.3.1+ is now required
* FEATURE: added `dng treecall` v2
* FEATURE: new k-alleles model for `dng call` that improves mutation calling
* BUGFIX: memory leak fixed with when using BCFtools parsing functions
* FEATURE: improved `dng pileup` now supports slurping regions and filter in number of alleles
* FEATURE: added region specification with BED files when using VCF on `dng call` 
* FEATURE: added support for AD/TAD format in `dng call` and `dng loglike`
* FEATURE: updated version of external dependencies
* FEATURE: improved support for regions in BAM files
* FEATURE: `dng dnm` can now handle indels 
* BUGFIX: fixed segfault found in some cases when using `GCC 4.8` 
* FEATURE: improved compiling on Apple
* FEATURE: optimized peeler implementation and interface

### v1.1.1

* FEATURE: When building with cmake, users have the option to download and install dependencies.
Just run `cmake -DBUILD_EXTERNAL_PROJECTS=1 ..`
* FEATURE: `dng help` now fully implemented
* FEATURE: more information added to `dng call` output
* BUGFIX: When running `dng call` with HTSLIB 1.2.1, the following error message was being
emitted: "FIXME: dirty header not synced".  This has been fixed.
* BUGFIX: `dng call` now outputs correct 1-based site locations.
* BUGFIX: a segfault was fixed in `dng dnm` and `dng phaser` caused by invalid commandline arguments
* CHANGE: HTSLIB 1.2+ is now required.
* Miscellaneous improvements to the build system

### v1.1

* Main program now called 'dng'
* Added experimental 'dng call' module.
* DeNovoGear now requires HTSLIB 1+, CMake 3.1+, Boost 1.47+, and Eigen 3+.

### v1.0

* made changes to indel_mrate parameter
* better indenting
* mu_scale scales indel mutation rate linearly

### v0.5.4

* added GPL v3
* updated output fields for indels, snps to be the same

### v0.5.3

* removed 'X' allele in VCF op. VCF can be indexed by Tabix, IGVTools and used in Annovar.
* added region based denovo calling on BCF files, invoked with --region flag
* added vcf parser for denovo calling, invoked with --vcf flag

### v0.5.2

* Added read-depth, posterior-probability filters.
* Output number of sites in the BCF and number of sites passing filters.
* Modified paired caller output.

### v0.5.1

* Fixed bug in triallelic configuration.
* Some trialleic denovo configurations were being called incorrectly.

### v0.5

* switched to cmake
* separate models for X chromosome calling in male offspring, female offspring

### v0.4

* paired sample analysis
* parental phaser

### v0.3

* INDEL length based mutation prior
* look at additional trio configs compared to initial version

### v0.2

* incorporated BCF parser
* denovo INDEL calling
