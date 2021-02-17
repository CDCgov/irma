# HEADER
PARAM_FILE_NAME="sensitive"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2020-04"

# CONSENSUS REFINEMENT & READ SELECTION
SKIP_E=1		# skip reference elongation
MIN_LEN=70		# relaxed threshold for UTR elongation with SAM
ASSEM_REF=0

# DEFAULT settings but with SGE/OGE/UGE execution turned on.
GRID_ON=1		# needs "qsub" and a NFS install of IRMA

# STAGES
MATCH_PROG="BLAT"
SORT_PROG="BLAT"
ALIGN_PROG="SAM"
ASSEM_PROG="SSW"
