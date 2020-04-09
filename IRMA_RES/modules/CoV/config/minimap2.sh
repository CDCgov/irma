# HEADER
PARAM_FILE_NAME="minimap2"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2020-04"

# Default settings but with minimap2 used for final assembly.
# Experimental. SSW still recommended for short reads.

# STAGES
DEL_TYPE="REF"			# don't delete in rough alignment
ALIGN_PROG="SAM"		# use HMMs for rough alignment
ASSEM_PROG="MINIMAP2"		# final assembly in minimap2

# match, mismatch, gap open, gap extend
MM2_A=2
MM2_B=5
MM2_O=10
MM2_E=1
