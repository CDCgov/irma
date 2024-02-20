# IRMA Changelog

## v1.1.4 : 2024-02

### Fixes

- Fixes a memory bottleneck for segmented modules (FLU) processing a large number of reads. Many thanks to T. Stark, K. Lacek, and M. Scotch for helping to diagnose the issue.

## v1.1.3 : 2023-11

### Fixes

- Fixes regression (fixed in v1.1.1 but reintroduced in v1.1.2) of the FASTQ header length issue for single-ended assemblies that can cause BAM files to not be written.

## v1.1.2 : 2023-10

### Fixes

- SRA format (with or without read side identifiers in the ID) now processes correctly without user FASTQ SEQNAME manipulation. IRMA still expects paired-end FASTQ files to be in split format and processes Illumina format as usual. Note: SRA read-side identifiers are order-based and their distribution may not be preserved between the original and SRA "standardized" sample.
- IRMA now correctly processes legacy Illumina format (FASTQ SEQNAME ends in '/1' and '/2') with respect to paired-end reads.
- IRMA no longer underlines FASTQ headers in the 'intermediate/4-ASSEMBLY-*/' sub-folder, which is sometimes used for submitting to SRA. Sorting is still by gene segment as needed and paired-end reads are still interleaved. The read set is also still adapter-trimmed as needed.

## v1.1.1 : 2023-09

### Config Changes

- Config file for "CoV-minion-long-reads" now uses `BLAT` instead of `SAM` for rough alignment during the read gathering phase. Results were similar in recent tests and performance is improved. If you need more sensitivity, use the "CoV-minion-sensitive" config file.

### Fixes

- Clarify semantic versioning standards established in v1.0.6 for program/module (minor) vs config file levels (patch)
- Single-ended assemblies no longer underscore spaces in the FASTQ SEQNAME as it can cause failures in SAM to BAM conversion due to excess length (aligners typically remove SEQNAME descriptions from the QNAME). The BAM specification limits the QNAME field to 254 bytes. If needed for paired-end, IRMA now truncates excess FASTQ header (SEQNAME) length.
- IRMA adds warning messages upon SAMTOOLS failure and for when FASTQ headers need to be truncated pre-emptively.
- The final BAM files were not being sorted and indexed as expected due to usage errors on Linux; this has been corrected.

## v1.1.0 : 2023-08

### Features

- Adds checked external configuration parsing using the `--external-config <FILE>` or `-c <FILE>` option. Example:

```bash
## Takes -c or --external-config at the beginning or end of the argument list
## Note that process substitution will not work with Docker!
./IRMA FLU ../flu-amd/tests/test2.fastq.gz S1 -c <(echo 'PARAM_FILE_AUTHOR="Wile E. Coyote"')
./IRMA FLU ../flu-amd/tests/test2.fastq.gz S1 --external-config config.sh
```

- Adds **opt-in** configuration `USE_IRMA_CORE` to run the *mergeSAMpairs* and *fastQ_converter* tasks more efficiently. This is currently most helpful for paired-end samples with a higher quantity of reads. The `irma-core_<OS>` binary must be in the *IRMA_RES/scripts* folder or in the `PATH`, otherwise the configuration will be turned back off. Example:

```bash
./IRMA FLU R1.tar.gz R2.tar.gz S1 -c <(echo 'USE_IRMA_CORE=1')
```

### Removed

- Deletes the obsolete defaults back-up file
- Removes superfluous configs in the defaults file and from loggers

### Fixes

- Removes an erroneous tab character from the line ends of the insertion tables. Can cause issues downstream when `ALIGN_AMENDED` is enabled (Thanks to K. Lacek).
- Fixed bug related to observed read-pair insertion discrepancies being ignored when calculating estimated insertion error.

## v1.0.6 : 2023-05

- When `ALLOW_TMP` is set but `TMP` is unset, adds an error and aborts to avoid writing to unintended directories.
- With IRMA v1.0.6 and beyond attempt we will attempt to follow *semantic versioning* (v Major.Minor.Patch) more closely:
  - **Major**: breaking changes to the CLI, output, configurations, and supported inputs
  - **Minor**: changes that break undocumented/unintended behaviors for generating data without altering format, adding new features, and changes in module or program configuration defaults
  - **Patch**: fix bugs, minor performance improvements, non-breaking updates to dependencies, log or error improvements, changes in config file defaults

## v1.0.5 : 2022-10

- Performs more efficient fuzzy adapter matching within IRMA for Illumina configurations.
- Provides a fix to actually allow turning off the default median read quality filter and instead use the average read quality filter as expected.  Configuration logs now correctly display USE_MEDIAN status.
- Adds more uniform and updated defaults to the CoV module configurations (thanks to K. Lacek).

## v1.0.4 : 2022-09

- Added new option "PACKAGED_FASTQ" (default on). One can now optionally turn off packaging together sorted assembly fastq as a "tar.gz" and instead gzip up each fastq separately.

## v1.0.3 : 2022-09 (thanks for contributions from K. Lacek & B. Rambo-Martin)

- Configuration tweaks for CoV
- Added S-gene ONLY configurations + references (v1.0.2p1; thanks to K. Lacek)
- Added unstable feature MIN_DROPOUT_EDGE_DEPTH (v1.0.2p1; off and may be removed, not recommended for regular use)
- Resolved a major memory issue with CoV assembly and Illumina paired-end read merging (v1.0.2p2)
- Added recombinant configurations + references for CoV (v1.0.2p3; thanks to K. Lacek)
- Updated samtools to v1.14 (Linux) that does not rely on ncurses
- Updated SSW to v0.1.5M (Linux+MacOS) and provided build script to repo fork: adds memory safety improvements
- Updated BLAT to v35 (Mac) and re-compiled for 64-bit

## v1.0.2 : 2021-03-11

- When available, left join HMM coordinates and alignment states from the coverage.a2m.txt file to the variants, allAlleles, insertions, and deletions table.  Many thanks to K. Lacek for assistance.

## v1.0.1 : 2021-03-03

- Refined amplicon dropout padding for greater accuracy; added plurality global alignment and coverge tables for A2M and padded data. Many thanks to K. Lacek for assistance.

## v1.0.0 : 2021-02-16 (many thanks for suggestions and feedback to T. Stark, K. Lacek, B. Rambo-Martin, C. Paden, and J. Tognarelli)

- Added SARS-CoV-2 to the experimental CoV module and performed optimization for MiSeq 2x150 PE reads.
- Added ENFORCE_CLIPPED_LENGTH to enforce minimum read length during QC post IRMA adapter trimmed Illumina PE reads.
- Added MIN_BLAT_MATCH to filter the BLAT matches during read gathering to the minimum length. Default progrm settings result in about 30bp even if set lower.
- Added MIN_CONS_SUPPORT to control total ambiguation ('N') of the amended consensus based on the plurality allele count.
- Added MIN_CONS_QUALITY to control total ambiguation ('N') of the amended consensus based on the plurality allele average quality.
- Added INS_T_DEPTH, the minimum insertion coverage depth, that allows insertion editing given INS_T (frequency minimum) isalso satisfied.
- Added DEL_T_DEPTH, the minimum deletion coverage depth, that allows deletion editing given DEL_T (frequency minimum) is also satisfied. The deletion must also be the plurality allele.
- Added more control for deletion editing (see DEL_TYPE)
- Added ALIGN_AMENDED option to do a global sequence of the amended consensus using the HMM profile.
- Added experimental support for minimap2 for the final assembly (ASSEM_PROG="MINIMAP2"). Helpful for long genomes + reads that Smith-Waterman can't handle so well. Split chains not yet supported.
- Added output directory (relative to working directory or absolute path) specification in the RUN argument.
- Updated version of GNU Parallel to 20200422 (more stable)
- Fixed bug with annotated (curly brackets) fastq, fuzzy adapter trimming in rare cases, and an issue affecting VCF file underreporting when AUTO_F=1 (variant tables remain correct).
- Fixed bug that could affect the performance and accuracy of BLAT alignment (though retroactively corrected in following rounds of read-gathering)
- Other bug fixes and stability improvements.

## v0.9.4 : 2019-12-04

- Heatmap defaults are now statically determined for consistency with later version of R (3.6+)

## v0.9.3 : 2019-09-12

- Increased local alignment end repair (read-gathering phase) from 7 to 9 bases (SAM/BLAT).

## v0.9.2 : 2019-09-06

- For direct RNA sequencing, always convert Uracil to Thymine to avoid unexpected side-effects
- Prevent whitespace in paths and inputs from crashing IRMA, improved error handling

## v0.9.1 : 2019-05-08

- Added the "FLU_AD" module, capable of assembling influenza C and D in addition to A and B
- Allow for more flexible LABEL module and SORT_GROUP specification; bug fixes

## v0.9.0 : 2019-03-20

- Added new options for "residual" and "secondary" assembly (currently on type and subtype level). Please read more.
- improved clean up code

## v0.8.4 : 2018-09-25

- Fixed Perl regex deprecation issue

## v0.8.3 : 2018-09-18

- Added variable "FUZZY_ADAPTER" (default ON) to trim library adapters with up to 1 mismatch

## v0.8.2 : 2018-07-09

- Simplified temporary directory randomization TOKEN to Perl only

## v0.8.1 : 2018-06-05

- Added variable "SILENCE_COMPLEX_INDELS" to silence reads with 4 or more indels within the final assembly

## v0.8.0 : 2018-03-23

- Adds new variable "GRID_PATH" for custom grid working directory

## v0.7.2 : 2018-03-23

- Fixes issues related to phase assignment, MATCH phase

## v0.7.1 : 2018-02-12

- IRMA can run in /tmp even if noexec set

## v0.7.0 : 2018-02-02

- Revised the EXPENR phase association measure to be more robust
- Added SNV phase assignment numbers to the variants files

## v0.6.8 : 2017-04-27

- Added disk space check based on fastq size to avoid running on a small disk
- Added fail-over from /tmp to the project directory

## v0.6.7 : 2017-04-14

- added parallelization for zipping files via pigz. If not compatible, defaults to gzip.
- updated the Ebola and CoV modules to use MERGE_SECONDARY=1

## v0.6.6 : 2017-04-13

- added "MERGE_SECONDARY" option, puts secondary data into the unmatched read pool for round 2. Recommended if lineage references are close, and for non-segmented viruses where co-infection is unlikely.
- added "GENE_GROUP" variable for customizable two-stage sorting with BLAT performing rough sorting during MATCH step and before the SORT step. Example: influenza "LABEL-trio" in the paper

## v0.6.5 : 2016-12-02

- added experimental configuration file "FLU-minion" for MinION R7/R9 chemistries
- fixed VCF bug that prevented correct VCF creation

## v0.6.4 : 2016-11-04

- added plurality consensus script to IRMA_RES/scripts folder to aid in new module creation

## v0.6.3 : 2016-08-12

- added failovers to local node computation when grid computation doesn't work

## v0.6.2 : 2016-06-01

- optimized indel calling algorithm for performance
- added flexible segment numbering option for module init.sh (see SEG_NUMBERS in configuration, and SORT_GROUPS in primary vs. secondary)
- added NO_SORT_REFS option for custom reference sorting into LABEL classes
- added experimental MERS CoV module
- changed NextTera adapter masking back to adapter clipping but without length penalization (see ADAPTER in configuration)
- fixed bug with flexible parameter file naming

## v0.6.1 : 2016-04-01

- added name + hashing to allow for multiple concurrent projects
- tweaked reference elongation algorithm
- changed ADAPTER clipping to masking

## v0.6.0 : 2016-03-09  (test release)
