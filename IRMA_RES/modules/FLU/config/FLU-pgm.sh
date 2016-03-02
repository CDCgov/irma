# HEADER
PARAM_FILE_NAME="FLU-pgm-flc"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2015-01-29"

# PERFORMANCE
MATCH_PROC=16				# grid maximum processes for the MATCH
ALIGN_PROC=16				# grid maximum processes for the rough align
ASSEM_PROC=16			# grid maximum processes for assembly

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=19			# average or median threshold for QUALITY reads
MIN_LEN=50				# minimum read length for QUALITY reads
