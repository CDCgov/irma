# HEADER
PARAM_FILE_NAME="CoV MinION Long Reads"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.2"
PARAM_FILE_DATE="2023-09"

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=20 # average or median threshold for QUALITY reads
USE_MEDIAN=1
MIN_LEN=125 # minimum read length for QUALITY reads
INS_T=0.5   # threshold for insertion refinement
DEL_T=0.75  # threshold for deletion refinement : 1 => turn OFF deletion editing
MIN_RP=1    # minimum read pattern count to continue
MIN_RC=3    # minimum read count to continue
MIN_CONS_SUPPORT=20
MIN_CONS_QUALITY=10

# VARIANT CALLING HEURISTICS & STATS
MIN_AQ=8    # minimum average variant quality, does not apply to deletions
MIN_FI=0.01 # minimum insertion variant frequency
MIN_FD=0.02 # minimum deletion variant frequency

ALIGN_PROG="BLAT" # rough alignment with BLAT
ASSEM_PROG="SSW"  # final assembly, if reads and template are greater than 28k each, consider MINIMAP2.
