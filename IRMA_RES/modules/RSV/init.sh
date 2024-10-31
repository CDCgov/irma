### PERFORMANCE ###
GRID_ON=0   # grid computation on [1,0] for on or off
ALLOW_TMP=1 # if GRID_ON=0, try to use /tmp for working directory
TMP=/tmp    # the scratch/tmpfs for working on the assemblies

### REFERENCE ###
MIN_FA=1         # no alternative reference [0..1]
MIN_CA=20        # minimum count for alternative finished assembly
SKIP_E=1         # skip reference elongation
REF_SET=$DEF_SET # Same as the "consensus.fasta" in the reference folder for the module.

### READ GATHERING ###
MAX_ROUNDS=2      # round of read gathering
USE_MEDIAN=1      # use the median quality or the average [1,0]
QUAL_THRESHOLD=19 # minimum read statistic
MIN_LEN=100       # minimum read length
MERGE_SECONDARY=1

## MATCH STEP
MATCH_PROC=64     # grid maximum processes for the MATCH
MATCH_PROG="BLAT" # match (all or any match) program [BLAT]
MIN_RP=15         # minimum read pattern count to continue
MIN_RC=15         # minimum read count to continue

## SORT STEP
SORT_GROUPS="__ALL__"
SORT_PROG="BLAT" # [LABEL,BLAT]
SORT_PROC=64     # currently not used
NONSEGMENTED=0   # segmented! [0,1]
# LABEL

## ALIGN STEP ##
ALIGN_PROG="BLAT" # rough assembly / alignment to working reference [SAM,BLAT]
ALIGN_PROC=64     # grid maximum processes for the rough align

### FINISHING ASSEMBLY ###
ASSEM_PROG="SSW" # assembly program [SSW]
ASSEM_PROC=64    # grid maximum processes for assembly
INS_T=0.25       # minimum frquenncy threshold for insertion refinement
DEL_T=0.60       # minimum frequency threshold for deletion refinement
MIN_AMBIG=0.25   # minimum called SNV frequency for mixed base in amended consensus folder
ALIGN_AMENDED=1
PADDED_CONSENSUS=1
ASSEM_REF=1

### VARIANT CALLING ###
# HEURISTICS
AUTO_F=1      # auto-adjust frequency threshold [1,0]
MIN_FI=0.01   # minimum insertion variant frequency
MIN_FD=0.01   # minimum deletion variant frequency
MIN_F=0.016   # minimum frequency for single nucleotide variants
MIN_C=5       # minimum count for variants
MIN_AQ=48     # minimum average variant quality, does not apply to deletions
MIN_TCC=100   # minimum non-ambiguous column coverage
MIN_CONF=0.80 # minimum confidence not machine error

# CONFIDENCE INTERVALS
SIG_LEVEL=0.999 # significance test level for variant calling (.90,.95,.99,.999).

PACKAGED_FASTQ=0
