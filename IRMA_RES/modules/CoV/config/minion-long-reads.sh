# HEADER
PARAM_FILE_NAME="CoV MinION Long Reads"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2020-04"

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=0			# average or median threshold for QUALITY reads
MIN_LEN=150				# minimum read length for QUALITY reads
INS_T=0.75				# threshold for insertion refinement
DEL_T=1.00				# threshold for deletion refinement : 1 => turn OFF deletion editing
MIN_RP=3				# minimum read pattern count to continue
MIN_RC=3				# minimum read count to continue

# VARIANT CALLING HEURISTICS & STATS
MIN_AQ=8			# minimum average variant quality, does not apply to deletions

DEL_TYPE="REF"		# rough alignment keeps reference
ALIGN_PROG="SAM"	# rough alignment with HMM
ASSEM_PROG="SSW"	# final assembly with MINIMAP2

MM2_A=2			# minimap2 match score
MM2_B=3			# minimap2 mismatch penalty
MM2_O=6			# minimap2 gap open penalty
MM2_E=1			# minimap2 gap extension penalty
