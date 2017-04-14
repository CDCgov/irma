# MODULE SPECIFIC
NONSEGMENTED=1				# non-segmented virus
SORT_GROUPS="__ALL__"			# Pattern list for sorting

# STAGES
MATCH_PROG="BLAT"
SORT_PROG="BLAT"
ALIGN_PROG="SAM"
ASSEM_PROG="SSW"

# PERFORMANCE
MATCH_PROC=140				# grid maximum processes for the MATCH
ALIGN_PROC=140				# grid maximum processes for the rough align
ASSEM_PROC=140				# grid maximum processes for assembly
SINGLE_LOCAL_PROC=16			# local maximum processes
DOUBLE_LOCAL_PROC=8			# local half maximum processes for doubled work
GRID_ON=1				# grid computation on
ALLOW_TMP=1				# use TMP for ppath
TMP=/tmp				# temperary space

# VARIANT CALLING HEURISTICS & STATS
MIN_FI=0.0045				# minimum insertion variant frequency
MIN_FD=0.0045				# minimum deletion variant frequency
MIN_F=0.015				# minimum frequency for variants
MIN_C=2					# minimum count for variants
MIN_AQ=24				# minimum average variant quality, does not apply to deletions
MIN_TCC=100				# minimum non-ambiguous column coverage
MIN_CONF=0.80				# minimum confidence not machine error
SIG_LEVEL=0.999				# significance test level for variant calling (.90,.95,.99,.999). 

# CONSENSUS REFINEMENT & READ SELECTION
MAX_ROUNDS=5
QUAL_THRESHOLD=30			# average or median threshold for QUALITY reads
MIN_LEN=125				# minimum read length for QUALITY reads
INS_T=0.15				# threshold for insertion refinement
DEL_T=0.50				# threshold for deletion refinement
SKIP_E=1				# skip reference elongation
INCL_CHIM=0				# whether or not to get rid of chimera
MIN_RP=15				# minimum read pattern count to continue
MIN_RC=15				# minimum read count to continue
MIN_AMBIG=0.25				# min SNV freq for ambig nts in final amended consensus
MERGE_SECONDARY=1

# ASSEMBLY
MAX_ITER_SSW=5				# max num of SSW iterations to perform, 3 should be sufficient w/4 to prove
SSW_M=2					# smith-waterman match score
SSW_X=5					# smith-waterman mismatch penalty
SSW_O=10				# smith-waterman gap open penalty
SSW_E=1					# smith-waterman gap extension penalty
