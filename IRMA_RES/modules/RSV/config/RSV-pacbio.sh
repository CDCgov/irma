# HEADER
# Use the read Of Inserts file or CCS reads
PARAM_FILE_NAME="RSV-PacBio"
PARAM_FILE_AUTHOR="X. Qiu"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2016-08-12"

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=20			# average or median threshold for QUALITY reads
MIN_LEN=125				# minimum read length for QUALITY reads
INS_T=0.25				# threshold for insertion refinement
DEL_T=0.50				# threshold for deletion refinement
MIN_RP=3				# minimum read pattern count to continue
MIN_RC=3				# minimum read count to continue

# VARIANT CALLING HEURISTICS & STATS
MIN_AQ=37				# minimum average variant quality, does not apply to deletions
