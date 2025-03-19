# HEADER
PARAM_FILE_NAME="CoV Recombinants"
PARAM_FILE_AUTHOR="K. Lacek"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2021-06"

MATCH_PROC=80 # grid maximum processes for the MATCH
SORT_PROC=80  # currently not used
ALIGN_PROC=80 # grid maximum processes for the rough align
ASSEM_PROC=80 # grid maximum processes for assembly

TMP=/tmp
ALIGN_AMENDED=1
PADDED_CONSENSUS=0
ASSEM_REF=1
GRID_ON=0
CUSTOM_REF_FILE="SARS-CoV-2_recombinant.fasta" # custom ref file
DEF_SET=$(dirname $DEF_SET)/$CUSTOM_REF_FILE
REF_SET=$DEF_SET
MERGE_SECONDARY=1
SORT_GROUPS="SARS-CoV-2"
