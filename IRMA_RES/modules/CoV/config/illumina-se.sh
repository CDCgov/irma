# HEADER
PARAM_FILE_NAME="Illumina short Single-end configuration"
PARAM_FILE_AUTHOR="K. Lacek"
PARAM_FILE_VERSION="1.1"
PARAM_FILE_DATE="2021-02"

SINGLE_LOCAL_PROC=48    # local maximum processes
DOUBLE_LOCAL_PROC=24    # local maximum processes (double this number)
MATCH_PROC=48           # grid maximum processes for the MATCH
SORT_PROC=48            # currently not used
ALIGN_PROC=48           # grid maximum processes for the rough align
ASSEM_PROC=48           # grid maximum processes for assembly

MIN_LEN=15
MIN_BLAT_MATCH=0
DEL_T_DEPTH=3
INS_T_DEPTH=3

ALIGN_PROG="BLAT"
ASSEM_PROG="SSW"
MAX_ROUNDS=5
MAX_ITER_ASSEM=2
