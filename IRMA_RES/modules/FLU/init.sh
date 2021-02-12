### PERFORMANCE ###
GRID_ON=0		# grid computation on [1,0] for on or off
GRID_PATH=""		# grid path, defaults to the IRMA_RES path if left empty string, do not include quotes for tilde prefix
SINGLE_LOCAL_PROC=16	# local maximum processes
DOUBLE_LOCAL_PROC=8	# local maximum processes (double this number)
ALLOW_TMP=1		# if GRID_ON=0, try to use /tmp for working directory
TMP=/tmp		# the scratch/tmpfs for working on the assemblies

### REFERENCE ###
MIN_FA=1		# no alternative reference [0..1]
MIN_CA=20		# minimum count for alternative finished assembly
SKIP_E=1		# skip reference elongation
REF_SET=$DEF_SET	# Same as the "consensus.fasta" in the reference folder for the module.

### READ GATHERING ###
MAX_ROUNDS=5			# round of read gathering
USE_MEDIAN=1			# use the median quality or the average [1,0]
QUAL_THRESHOLD=30		# minimum read statistic
MIN_LEN=125			# minimum read length
ENFORCE_CLIPPED_LENGTH=0	# Off. Reads are filtered for minimum length post adapter trimming.	

## MATCH STEP
MATCH_PROC=20		# grid maximum processes for the MATCH
MATCH_PROG="BLAT"	# match (all or any match) program [BLAT]
MIN_RP=15		# minimum read pattern count to continue
MIN_RC=15		# minimum read count to continue

## SORT STEP 
SORT_PROG="BLAT"	# [LABEL,BLAT]
SORT_PROC=80		# currently not used
NONSEGMENTED=0
# Pattern list to group gene segment lineages into gene for primary/secondary sorting.
SORT_GROUPS="PB2,PB1,PA,HA,NP,NA,MP,NS"
SEG_NUMBERS="B_PB1:1,B_PB2:2,A_PB2:1,A_PB1:2,PA:3,HA:4,NP:5,NA:6,M:7,NS:8"

# LABEL
SECONDARY_SORT=1						# LABEL sorting fast-mode
LABEL_MODULE="irma-FLU"						# if LABEL SECONDARY SORT is 0, use LABEL_MODULE
SECONDARY_LABEL_MODULES="irma-FLU-HA,irma-FLU-NA:irma-FLU-OG"	# otherwise, search for primary classification from BLAT and use the modules accordingly
GENE_GROUP="HA,NA:OG"						# specify primary sorting gene groups for BLAT

## ALIGN STEP ##
ALIGN_PROG="SAM"	# rough assembly / alignment to working reference [SAM,BLAT]
ALIGN_PROC=20		# grid maximum processes for the rough align

### FINISHING ASSEMBLY ###
ASSEM_PROG="SSW"	# assembly program [SSW]
ASSEM_PROC=20		# grid maximum processes for assembly
INS_T=0.25		# minimum frquenncy threshold for insertion refinement
DEL_T=0.60		# minimum frequency threshold for deletion refinement 
MIN_AMBIG=0.20		# minimum called SNV frequency for mixed base in amended consensus folder


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
