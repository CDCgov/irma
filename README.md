# IRMA - the Iterative Refinement Meta-Assembler

IRMA is a highly configurable next-generation sequencing assembly pipeline for virus genomes. Seed references are refined and edited as reads are matched, sorted, and aligned across virus genomes or gene segments. More information can be found on the IRMA website: https://wonder.cdc.gov/amd/flu/irma/

## Usage

IRMA takes a module-configuration name, one or two fastq, and a project name. The project name can be a full path. For example:

```{bash}
IRMA <MODULE|MODULE-CONFIG> <R1.fastq.gz|R1.fastq> <R2.fastq.gz|R2.fastq> [path/to/]<sample_name>
IRMA <MODULE|MODULE-CONFIG> <fastq|fastq.gz> [path/to/]<sample_name>
```

More usage information: https://wonder.cdc.gov/amd/flu/irma/run.html

## Attributions

Please cite IRMA if you use it in your manuscript:

> Shepard SS, Meno S, Bahl J, Wilson MM, Barnes J, Neuhaus E. Viral deep sequencing needs an adaptive approach: IRMA, the iterative refinement meta-assembler. BMC Genomics. 2016;17(1). doi:10.1186/s12864-016-3030-6

Read the paper: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6

## License

As a complete package, IRMA is under GPLv3 with non-commercial use clauses added (due to use of BLAT in IRMA and SAM in LABEL). See also: https://wonder.cdc.gov/amd/flu/irma/disclaimer.html

## System requirements

Perl 5.16.1, R 3.6+, BASH 3+; CentOS 7+, Ubuntu Linux, or MacOS 10.14+