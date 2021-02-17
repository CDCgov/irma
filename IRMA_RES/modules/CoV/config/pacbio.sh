# HEADER
# Use the read Of Inserts file or CCS reads
PARAM_FILE_NAME="PacBio"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="2.0"
PARAM_FILE_DATE="2021-01"

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=20			# average or median threshold for QUALITY reads
MIN_LEN=125				# minimum read length for QUALITY reads
INS_T=0.75				# threshold for insertion refinement
DEL_T=0.75				# threshold for deletion refinement
MIN_RP=1				# minimum read pattern count to continue
MIN_RC=3				# minimum read count to continue
MIN_CONS_SUPPORT=3      		# minimum consensus alleles to not ambiguate

#ALIGNMENT STEP
ALIGN_PROG="SAM"

# FINISHING ASSEMBLY
ASSEM_PROG="MINIMAP2"

# match, mismatch, gap open, gap extend
MM2_A=2
MM2_B=5
MM2_O=10
MM2_E=1

# VARIANT CALLING HEURISTICS & STATS
MIN_AQ=37		# minimum average variant quality, does not apply to deletions
MIN_FI=0.02		# minimum insertion variant frequency
MIN_FD=0.01		# minimum deletion variant frequency
