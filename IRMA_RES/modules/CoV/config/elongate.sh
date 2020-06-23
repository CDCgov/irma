# HEADER
PARAM_FILE_NAME="elongate"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2020-04"

# Performs reference elongation

# CONSENSUS REFINEMENT & READ SELECTION
INS_T=0.25				# threshold for insertion refinement
DEL_T=1.00				# threshold for deletion refinement
SKIP_E=0				# skip reference elongation
MIN_LEN=75				# relaxed threshold for UTR elongation with SAM
MAX_ROUNDS=5				# more rounds than default needed for elongation

# STAGES
DEL_TYPE="REF NNN"
MATCH_PROG="BLAT"
SORT_PROG="BLAT"
ALIGN_PROG="SAM BLAT"
ASSEM_PROG="SSW"
