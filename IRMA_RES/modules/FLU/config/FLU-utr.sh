# HEADER
PARAM_FILE_NAME="FLU"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2015-01-29"

# CONSENSUS REFINEMENT & READ SELECTION
INS_T=0.25				# threshold for insertion refinement
DEL_T=0.60				# threshold for deletion refinement
SKIP_E=0				# skip reference elongation
MIN_LEN=125				# relaxed threshold for UTR elongation with SAM

# STAGES
MATCH_PROG="BLAT"
SORT_PROG="BLAT"
ALIGN_PROG="SAM BLAT"
ASSEM_PROG="SSW"
