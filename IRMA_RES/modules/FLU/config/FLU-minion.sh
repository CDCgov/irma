# HEADER
PARAM_FILE_NAME="FLU MinION"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2016-08-23"

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=0			# average or median threshold for QUALITY reads
MIN_LEN=150				# minimum read length for QUALITY reads
INS_T=0.75				# threshold for insertion refinement
DEL_T=0.75				# threshold for deletion refinement
MIN_RP=3				# minimum read pattern count to continue
MIN_RC=3				# minimum read count to continue
MIN_CONS_SUPPORT=50
MIN_CONS_QUALITY=10

# VARIANT CALLING HEURISTICS & STATS
MIN_AQ=8			# minimum average variant quality, does not apply to deletions

SORT_PROG="BLAT"
ALIGN_PROG="SAM BLAT"

SSW_M=2			# smith-waterman match score
SSW_X=3			# smith-waterman mismatch penalty
SSW_O=6			# smith-waterman gap open penalty
SSW_E=1			# smith-waterman gap extension penalty
