# HEADER
PARAM_FILE_NAME="Permissive mode"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.1"
PARAM_FILE_DATE="2021-02"

INS_T_DEPTH=1		# minimum coverage depth for insertion refinement
DEL_T_DEPTH=1		# minimum coverage depth for deletion refinement 

DEL_TYPE="REF NNN"
ALIGN_PROG="BLAT"
ASSEM_PROG="SSW"
MAX_ROUNDS=5

MIN_BLAT_MATCH=0
QUAL_THRESHOLD=10	# average or median threshold for QUALITY reads
MIN_LEN=50		# minimum read length for QUALITY reads
MIN_RP=1		# minimum read pattern count to continue
MIN_RC=1		# minimum read count to continue

MIN_CONS_SUPPORT=1	# consensus allele minimum count
MIN_CONS_QUALITY=1	# consensus allele minimum average quality
