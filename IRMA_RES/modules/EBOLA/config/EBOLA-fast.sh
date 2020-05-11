# HEADER
PARAM_FILE_NAME="EBOLA-fast"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2015-04-17"

# STAGES
MATCH_PROG="BLAT"
SORT_PROG="BLAT"
ALIGN_PROG="BLAT"
ASSEM_PROG="SSW"

# READ GATHERING
MAX_ROUNDS=2

# ASSEMBLY
MAX_ITER_ASSEM=3			# max num of SSW iterations to perform, 3 should be sufficient
SSW_M=2					# smith-waterman match score
SSW_X=5					# smith-waterman mismatch penalty
SSW_O=10				# smith-waterman gap open penalty
SSW_E=1					# smith-waterman gap extension penalty
