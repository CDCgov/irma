# HEADER
PARAM_FILE_NAME="sensitive"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2020-04"

# CONSENSUS REFINEMENT & READ SELECTION
INS_T=0.25				# threshold for insertion refinement
DEL_T=0.75				# threshold for deletion refinement
SKIP_E=1				# skip reference elongation
MIN_LEN=70				# relaxed threshold for UTR elongation with SAM
MAX_ROUNDS=5				# diminishing returns fo later rounds, but leave no read behind

# DEFAULT settings but with SGE/OGE/UGE execution turned on.
GRID_ON=1		# needs "qsub" and a NFS install of IRMA

# STAGES
DEL_TYPE="REF"
MATCH_PROG="BLAT"
SORT_PROG="BLAT"
ALIGN_PROG="SAM"
ASSEM_PROG="SSW"
