### PERFORMANCE ###
GRID_ON=0		# grid computation on [1,0] for on or off
SINGLE_LOCAL_PROC=${NSLOTS:-32}	# local maximum processes
DOUBLE_LOCAL_PROC=$((SINGLE_LOCAL_PROC / 2))	# local maximum processes (double this number)
ALLOW_TMP=1		# if GRID_ON=0, try to use /tmp for working directory
TMP=/tmp		# the scratch/tmpfs for working on the assemblies

### REFERENCE ###
MIN_FA=1		# no alternative reference [0..1]
MIN_CA=20		# minimum count for alternative finished assembly
SKIP_E=1		# skip reference elongation
REF_SET=$DEF_SET	# Same as the "consensus.fasta" in the reference folder for the module.

### READ GATHERING ###
MAX_ROUNDS=5		# round of read gathering
USE_MEDIAN=1		# use the median quality or the average [1,0]
QUAL_THRESHOLD=27	# minimum read statistic. May wish to set lower for 2x300 reads.
MIN_LEN=80		# minimum read length
MERGE_SECONDARY=1

## MATCH STEP
MATCH_PROC=64		# grid maximum processes for the MATCH
MATCH_PROG="BLAT"	# match (all or any match) program [BLAT]
MIN_RP=1		# minimum read pattern count to continue
MIN_RC=15		# minimum read count to continue
MIN_BLAT_MATCH=65	# minimum blat match, default settings within the program practically limit to 30 bp, only useful if set higher.

## SORT STEP
SORT_GROUPS="__ALL__"
SORT_PROG="BLAT"	# [LABEL,BLAT]
SORT_PROC=64		# currently not used
NONSEGMENTED=0		# segmented! [0,1]
# LABEL
SECONDARY_SORT=0	# LABEL sorting fast-mode

## ALIGN STEP ##
ALIGN_PROG="BLAT"	# rough assembly / alignment to working reference [SAM,BLAT]
ALIGN_PROC=64		# grid maximum processes for the rough align
DEL_TYPE="REF"		# how to handle deletions in the rough alignment: NNN, REF, or DEL (blank = OLD DEFAULTS)

### FINISHING ASSEMBLY ###
ASSEM_PROG="SSW"	# assembly program [SSW]
ASSEM_PROC=64		# grid maximum processes for assembly
INS_T=0.75		# minimum frquenncy threshold for insertion refinement
DEL_T=0.75		# minimum frequency threshold for deletion refinement 
INS_T_DEPTH=15		# minimum coverage depth for insertion refinement
DEL_T_DEPTH=10		# minimum coverage depth for deletion refinement (in addition to plurality and frequency)
MIN_AMBIG=0.20		# minimum called SNV frequency for mixed base in amended consensus folder
ASSEM_REF=1		# use the same reference(s) for the final assembly as read gathering
			# NOTE: uses the default reference seed, any mention of "user supplied reference sorting" is merely perfunctory
ALIGN_AMENDED=1		# align the amended consensus to the HMM profile
PADDED_CONSENSUS=1	# attempt to pad amended_consensus with Ns for amplicaton dropout: requires ALIGN_AMENDED=1 and ASSEM_REF=1
MIN_CONS_SUPPORT=15	# consensus allele minimum count
MIN_CONS_QUALITY=15	# consensus allele minimum average quality

### VARIANT CALLING ###
# HEURISTICS
AUTO_F=1		# auto-adjust frequency threshold [1,0]
MIN_FI=0.005		# minimum insertion variant frequency
MIN_FD=0.005		# minimum deletion variant frequency
MIN_F=0.008		# minimum frequency for single nucleotide variants
MIN_C=2			# minimum count for variants
MIN_AQ=24		# minimum average variant quality, does not apply to deletions
MIN_TCC=100		# minimum non-ambiguous column coverage
MIN_CONF=0.80		# minimum confidence not machine error

# CONFIDENCE INTERVALS
SIG_LEVEL=0.999		# significance test level for variant calling (.90,.95,.99,.999). 

