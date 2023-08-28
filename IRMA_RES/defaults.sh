### PERFORMANCE ###
GRID_ON=0            # grid computation on [1,0] for on or off
GRID_PATH=""         # grid path, defaults to the IRMA_RES path if left empty string, do not include quotes for tilde prefix
LIMIT_BLAT=60000     # threshold before grid
LIMIT_SSW=80000      # threshold before grid
LIMIT_SAM=500        # threshold before grid
SINGLE_LOCAL_PROC=16 # local maximum processes
DOUBLE_LOCAL_PROC=8  # local maximum processes (double this number)
ALLOW_TMP=1          # if GRID_ON=0, try to use /tmp for working directory
TMP=/tmp             # the scratch/tmpfs for working on the assemblies
IRMA_QUEUE=""        # queue for qsub, e.g.: QUEUE="-q all.q"  Leave unset normally.
ALLOW_DISK_CHECK=1   # check for disk free space before running, relies on Perl & unix df

### REFERENCE ###
MIN_FA=1         # no alternative reference [0..1]
MIN_CA=20        # minimum count for alternative finished assembly
SKIP_E=1         # skip reference elongation
REF_SET=$DEF_SET # Same as the "consensus.fasta" in the reference folder for the module.
ASSEM_REF=0      # start with same reference set for the final assembly, sort if custom file
NO_SORT_REFS=0   # do not sort refs for custom references

### READ GATHERING ###
MAX_ROUNDS=5      # round of read gathering
USE_MEDIAN=1      # use the median quality or the average [1,0]
QUAL_THRESHOLD=30 # minimum read statistic
MIN_LEN=125       # minimum read length
INCL_CHIM=0       # includes chimera or not [0,1]
# transposase adapter, clips 5' of the adapter on the forward strand and 3' on the reverse strand
# applicable to NextTera pair-end reads
ADAPTER="AGATGTGTATAAGAGACAG"
FUZZY_ADAPTER=1            # Allow 1 adapter mismatch
ENFORCE_CLIPPED_LENGTH=1   # Reads are filtered for minimum length post adapter trimming.
MERGE_SECONDARY=0          # Merge secondary data after the first round to the unmatched pool. Useful if co-infections are not expected.
RESIDUAL_ASSEMBLY_FACTOR=0 # The ratio of match read patterns to altmatch is measured, and if it's less than factor, residual assembly is performed. Set 0 for off.

## MATCH STEP
MATCH_PROC=20     # grid maximum processes for the MATCH
MATCH_PROG="BLAT" # match (all or any match) program [BLAT]
MIN_RP=15         # minimum read pattern count to continue
MIN_RC=15         # minimum read count to continue
MIN_BLAT_MATCH=0  # minimum blat match, default settings within the program practically limit to 30 bp, only useful if set higher.

## SORT STEP
SORT_PROG="BLAT"            # [LABEL,BLAT]
SORT_PROC=80                # currently not used
SORT_GROUPS=""              # if blank, default is not to group any genes into primary & secondary
BAN_GROUPS="UNRECOGNIZABLE" # patterns not allowed, useful for LABEL
NONSEGMENTED=0              # segmented! [0,1]
# LABEL
LIMIT_LABEL=1000 # limit on LABEL grid activation (self-adjusts)
GENE_GROUP=""    # item1,item2:otherwise
SECONDARY_SORT=0 # sort via BLAT then the SORT program
LABEL_MODULE=""  # primary LABEL sort

## ALIGN STEP ##
ALIGN_PROG="SAM" # rough assembly / alignment to working reference [SAM,BLAT]
ALIGN_PROC=20    # grid maximum processes for the rough align'
DEL_TYPE=""      # how to handle deletions in the rough alignment: NNN, REF, or DEL (blank = OLD DEFAULTS)

### FINISHING ASSEMBLY ###
MAX_ITER_ASSEM=5         # max assembly iteration [5]
NO_MERGE=0               # do not merge read pairs [0]
ASSEM_PROG="SSW"         # assembly program [SSW]
ASSEM_PROC=20            # grid maximum processes for assembly
PHASE_PROC=80            # grid processes for variant phasing
LIMIT_PHASE=200          # number of variants before using grid
INS_T=0.25               # minimum frquenncy threshold for insertion refinement
DEL_T=0.60               # minimum frequency threshold for deletion refinement
INS_T_DEPTH=1            # minimum coverage depth for insertion refinement
DEL_T_DEPTH=1            # minimum coverage depth for deletion refinement
SILENCE_COMPLEX_INDELS=0 # silences reads with complex indels (having 4 indels with at least one greater than 9)
MIN_AMBIG=0.25           # minimum called SNV frequency for mixed base in amended consensus folder
MIN_CONS_SUPPORT=""      # consensus allele minimum count, used to amend final consensus as 'N' if support fails. Blank is off.
MIN_CONS_QUALITY=""      # consensus allele minimum average quality, used to amend final consensus as 'N' if threshold fails. Blank is off.
SSW_M=2                  # smith-waterman match score
SSW_X=5                  # smith-waterman mismatch penalty
SSW_O=10                 # smith-waterman gap open penalty
SSW_E=1                  # smith-waterman gap extension penalty
SEG_NUMBERS=""           # amended consensus suffix for segments (if applicable). list separated by commas and key-value separated by colon.
ALIGN_AMENDED=0          # do global alignment of the amended consensus to the HMM profile
PADDED_CONSENSUS=0       # attempt to pad amended_consensus with Ns for amplicaton dropout: requires ALIGN_AMENDED=1 and ASSEM_REF=1

### VARIANT CALLING ###
# HEURISTICS
AUTO_F=1          # auto-adjust frequency threshold [1,0]
MIN_FI=0.005      # minimum insertion variant frequency
MIN_FD=0.005      # minimum deletion variant frequency
MIN_F=0.008       # minimum frequency for single nucleotide variants
MIN_C=2           # minimum count for variants
MIN_AQ=24         # minimum average variant quality, does not apply to deletions
MIN_TCC=100       # minimum non-ambiguous column coverage
MIN_CONF=0.80     # minimum confidence not machine error
PHASE_DISTANCE="" # threshold used for EXPENRD to assign phases, leave blank for default

# CONFIDENCE INTERVALS
SIG_LEVEL=0.999 # significance test level for variant calling (.90,.95,.99,.999).
