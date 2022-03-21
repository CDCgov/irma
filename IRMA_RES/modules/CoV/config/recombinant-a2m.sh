# HEADER
PARAM_FILE_NAME="CoV Recombinants"
PARAM_FILE_AUTHOR="K. Lacek"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2021-06"

SINGLE_LOCAL_PROC=16    # local maximum processes
DOUBLE_LOCAL_PROC=8    # local maximum processes (double this number)
MATCH_PROC=80           # grid maximum processes for the MATCH
SORT_PROC=80            # currently not used
ALIGN_PROC=80           # grid maximum processes for the rough align
ASSEM_PROC=80           # grid maximum processes for assembly

TMP=/tmp
ALIGN_AMENDED=1
PADDED_CONSENSUS=0
#MIN_DROPOUT_EDGE_DEPTH=3
ASSEM_REF=1 
GRID_ON=0
CUSTOM_REF_FILE="SARS-CoV-2_recombinant.fasta"   # custom ref file
DEF_SET=$(dirname $DEF_SET)/$CUSTOM_REF_FILE
REF_SET=$DEF_SET
MERGE_SECONDARY=1
SORT_GROUPS="SARS-CoV-2"
