# HEADER
PARAM_FILE_NAME="RSV-pgm"
PARAM_FILE_AUTHOR="X. Qiu"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2016-08-12"

GRID_ON=0

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=19			# average or median threshold for QUALITY reads
MIN_LEN=100				# minimum read length for QUALITY reads

MIN_FI=0.02				# minimum insertion variant frequency
MIN_FD=0.03				# minimum deletion variant frequency

# STAGES
MATCH_PROG="BLAT"
SORT_PROG="BLAT"
ALIGN_PROG="SAM BLAT"
ASSEM_PROG="SSW"
