# HEADER
PARAM_FILE_NAME="FLU"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2015-01-29"

# PERFORMANCE
GRID_ON=0				# grid computation on
MATCH_PROC=20				# grid maximum processes for the MATCH
ALIGN_PROC=20				# grid maximum processes for the rough align
ASSEMBLE_PROC=20			# grid maximum processes for assembly
SINGLE_LOCAL_PROC=12			# local maximum processes
DOUBLE_LOCAL_PROC=6			# local half maximum processes for doubled work

# VARIANT CALLING HEURISTICS & STATS
MIN_FI=0.005				# minimum insertion variant frequency
MIN_FD=0.005				# minimum deletion variant frequency
MIN_F=0.008				# minimum frequency for variants
MIN_C=2					# minimum count for variants
MIN_AQ=24				# minimum average variant quality, does not apply to deletions
MIN_TCC=100				# minimum non-ambiguous column coverage
MIN_CONF=0.80				# minimum confidence not machine error
SIG_LEVEL=0.999				# significance test level for variant calling (.90,.95,.99,.999). 

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=30			# average or median threshold for QUALITY reads
MIN_LEN=125				# minimum read length for QUALITY reads
INS_T=0.15				# threshold for insertion refinement
DEL_T=0.75				# threshold for deletion refinement
SKIP_E=1				# skip reference elongation
INCL_CHIM=0				# whether or not to get rid of chimera
MIN_RP=15				# minimum read pattern count to continue
MIN_RC=15				# minimum read count to continue
MIN_AMBIG=0.25				# min SNV freq for ambig nts in final amended consensus

# ASSEMBLY
MAX_ITER_SSW=5				# max num of SSW iterations to perform, 3 should be sufficient w/4 to prove
SSW_M=2					# smith-waterman match score
SSW_X=5					# smith-waterman mismatch penalty
SSW_O=10				# smith-waterman gap open penalty
SSW_E=1					# smith-waterman gap extension penalty

# DO NOT ALTER
NONSEGMENTED=0				# segmented!
LFASTM=1				# LABEL sorting fast-mode

# STAGES
MATCH_PROG="BLAT"
SORT_PROG="BLAT"
ALIGN_PROG="BLAT"
ASSEM_PROG="SSW"
