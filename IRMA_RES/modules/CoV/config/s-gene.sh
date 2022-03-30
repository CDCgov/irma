# HEADER
PARAM_FILE_NAME="CoV S-gene"
PARAM_FILE_AUTHOR="K. Lacek"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2021-02"
#Config developed with Nanopore Data

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=0        # average or median threshold for QUALITY reads
MIN_LEN=125             # minimum read length for QUALITY reads
INS_T=0.5              # threshold for insertion refinement
DEL_T=0.5              # threshold for deletion refinement : 1 => turn OFF deletion editing
MIN_RP=1                # minimum read pattern count to continue
MIN_RC=3                # minimum read count to continue
MIN_CONS_SUPPORT=3
#MIN_DROPOUT_EDGE_DEPTH=3
ASSEM_REF=1

# VARIANT CALLING HEURISTICS & STATS
MIN_AQ=8                # minimum average variant quality, does not apply to deletions
MIN_FI=0.01             # minimum insertion variant frequency
MIN_FD=0.02             # minimum deletion variant frequency

ALIGN_PROG="BLAT"        # rough alignment with HMM
ASSEM_PROG="SSW"        # final assembly, if reads and template are greater than 28k each, consider MINIMAP2.

PADDED_CONSENSUS=1
CUSTOM_REF_FILE="SARS-CoV-2_spike.fasta"   # custom ref file
REF_SET=$(dirname $REF_SET)/$CUSTOM_REF_FILE
