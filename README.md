# IRMA - the Iterative Refinement Meta-Assembler

IRMA is a highly configurable next-generation sequencing assembly pipeline for virus genomes. Seed references are refined and edited as reads are matched, sorted, and aligned across virus genomes or gene segments. More information can be found on the [IRMA website](https://wonder.cdc.gov/amd/flu/irma/) or you can read the [manuscript]. If you use IRMA in a paper, please [cite](CITATION.bib) us.

## Usage

IRMA takes a `module`-`configuration` name, one or two fastq, and a project name. The project name can be a full path. For example:

```bash
Usage:
(PAIRED-END):   IRMA <MODULE|MODULE-CONFIG> <R1.fastq.gz|R1.fastq> <R2.fastq.gz|R2.fastq> [path/to/]<sample_name> [options]
(SINGLE-END):   IRMA <MODULE|MODULE-CONFIG> <fastq|fastq.gz> [path/to/]<sample_name> [options]

Options:
        --external-config|-c <VALID_CONFIG_PATH>
```

More usage information: <https://wonder.cdc.gov/amd/flu/irma/run.html>

## INSTALLATION & REQUIREMENTS

We recommend a multi-core machine with no fewer than 8 cores (16 or more threads work best) and at least 16 GB of RAM.  IRMA runtime is impacted by the number of cores available on a machine. In addition software requirements include:

- Perl 5.16
- R 3.6+
- BASH 3+
- Linux: RHEL 8+ or any recent Ubuntu like Bookworm
- *OR* macOS 10.14 on intel, macOS 11 on arm64 (Rosetta 2 required for LABEL)

### Via Archive

Download the latest archive via our [releases page](https://github.com/CDCgov/irma/releases). Use of `wget` or `curl` for downloads is *recommended for macOS to preserve functionality*.

1) Unzip the archive containing IRMA.
2) Move the package to your desired location and add the folder to your `PATH`
   - Note: IRMA_RES and IRMA must be in the same folder.
3) IRMA is now installed.  To test it from the package folder, execute:

   ```bash
   ./IRMA FLU tests/test2.fastq.gz test_project
   ```

### Via Docker

Simply run:

```bash
## From Github
docker run --rm -itv $(pwd):/data ghcr.io/cdcgov/irma:latest IRMA # more args

## From Dockerhub
docker run --rm -itv $(pwd):/data docker/cdcgov/irma:latest LABEL # more args
```

### Via MIRA

For a GUI interface, please visit [MIRA](https://github.com/CDCgov/MIRA). For usage with [nextflow](https://nextflow.io) consider [MIRA-NF](https://github.com/CDCgov/MIRA-NF). MIRA provides more end-to-end functionality out-of-the-box by integrating other components and providing extra summaries.

## Third Party Software

We aggregate and provide [builds of 3rd party software](IRMA_RES/third_party/) for execution at runtime with IRMA. Similarly, still more third-party software is [co-distributed with LABEL][label-codistribute] for runtime with IRMA. You may install or obtain your own copies and IRMA will detect them, but the user will be required to test for compatibility.

### In IRMA proper

- [blat] v35
  - Artifacts: `blat`
  - Purpose: match, sort, and rough align phase
  - License: [noncommerical license][blat-license] (licensed for personal,
    academic, and non-profit) with [permission to redistribute][blat-permission]
    in IRMA
- [minimap2] v2.29
  - Artifacts: `minimap2`
  - Purpose: final assembly phase
  - License: [MIT]
- [GNU Parallel] v20200422
  - Artifacts: `parallel`
  - Requires: system Perl
  - Purpose: local parallelization
  - License: [GPL v3]
- [pigz] v2.8
  - Artifacts: `pigz`
  - Purpose: zipping FASTQ, but can also use `gzip`
  - License: [free to redistribute][pigz-license]
- [samtools] v1.21
  - Artifacts: `samtools`
  - Purpose: converting SAM to BAM, indexing BAM, sorting BAM
  - License: [MIT/Expat][samtools-license]
- [SSW] v1.2.5M
  - Artifacts: `ssw`
  - Custom modifications:
    <https://github.com/sammysheep/Complete-Striped-Smith-Waterman-Library/tree/IRMA%40v1.3>
  - Purpose: final assembly phase
  - License: [MIT]

### Co-distributed under LABEL

- [GNU Parallel]
- [SHOGUN] version 1.1.0 (2.1+ is not compatible)
  - Artifacts: `shogun` (cmdline_static)
  - Exception: apple/arm64 requires **Rosetta2**
  - Purpose: executes the SVM decision phase.
  - License: [GPL v3]
- [SAM] version 3.5
  - Artifacts: `align2model`, `hmmscore`, `modelfromalign`
  - Purpose: build HMM profiles, score sequences for evaluation
  - License: [Custom][sam-license] academic/government, not-for-profit,
    redistributed [with permission][sam-permission]

> [!WARNING]
> Note that [blat] and [SAM] (via LABEL) is redistributed with
> permission but their terms exclude commerical use without a license. If you
> are a commercial entity, you might need to reach out to the respective authors
> (see example [blat-license] and [sam-license]) to obtain a license for commercial use.

## Contributing

New feature integration is mainly being moved to the [irma-core](https://github.com/CDCgov/irma-core) project (support by [Zoe](https://github.com/CDCgov/zoe)), which is a mandatory component of IRMA as of v1.3.0. The IRMA repo itself will continue to address fixes needed for the Perl, R, and BASH code base, albeit some pieces may be obviated in future versions by IRMA-core. This repo will also continue to be responsible for data updates and other module artifacts and for the overall build and release.

## Notices

### Contact Info

For direct correspondence on the project, feel free to contact: [Samuel S. Shepard](mailto:sshepard@cdc.gov), Centers for Disease Control and Prevention or reach out to the [IRMA contributors](CONTRIBUTORS.md).

### Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).  All contributions to this repository will be released under the CC0 dedication.  By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

### License Standard Notice

The repository utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later. This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version. This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details. You should have received a copy of the Apache Software License along with this program. If not, see: <http://www.apache.org/licenses/LICENSE-2.0.html>. The source code forked from other open source projects will inherit its license.

### Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the [Disclaimer](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md). For more information about CDC's privacy policy, please visit <http://www.cdc.gov/other/privacy.html>.

### Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo) and submitting a pull request. (If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git).) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the [Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).

[label-codistribute]: https://github.com/CDCgov/label/blob/master/LABEL_RES/third_party/MANIFEST.md
[GNU Parallel]: https://www.gnu.org/software/parallel/
[sam-permission]: https://github.com/CDCgov/label/tree/master/LABEL_RES/third_party/copyright_and_licenses/sam3.5/SAM%20Redistribution%20Special%20Permissions.pdf
[GPL v3]: https://www.gnu.org/licenses/gpl-3.0.txt
[SHOGUN]: https://github.com/shogun-toolbox/
[sam-license]: https://users.soe.ucsc.edu/~karplus/projects-compbio-html/sam-lic/obj.0
[SAM]: https://users.soe.ucsc.edu/~karplus/projects-compbio-html/sam2src/
[MIT]: https://opensource.org/license/mit
[SSW]: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
[minimap2]: https://github.com/lh3/minimap2
[blat]: https://kentinformatics.com
[blat-license]: https://wonder.cdc.gov/amd/flu/irma/blat-license-noncommercial.html
[blat-permission]: https://github.com/CDCgov/irma/tree/master/IRMA_RES/third_party/copyright_and_licenses/blat-permission-redistribution.pdf
[pigz]: https://github.com/madler/pigz
[pigz-license]: https://raw.githubusercontent.com/madler/pigz/refs/heads/master/README
[samtools]: https://github.com/samtools/samtools
[samtools-license]: https://raw.githubusercontent.com/samtools/samtools/refs/heads/develop/LICENSE
[manuscript]: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6