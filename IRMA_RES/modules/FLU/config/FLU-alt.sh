# HEADER
PARAM_FILE_NAME="FLU-alt"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2018-02"

# VARIANT CALLING HEURISTICS & STATS
MIN_FA=0.050				# minimum frequency for alternative reference
MIN_CA=20				# minimum count for alternative reference

ALIGN_PROG="SAM"
SORT_PROG="LABEL"
MATCH_PROG="BLAT"

# CONSENSUS REFINEMENT & READ SELECTION
INS_T=0.50				# threshold for insertion refinement
DEL_T=0.60				# threshold for deletion refinement

GRID_ON=0
