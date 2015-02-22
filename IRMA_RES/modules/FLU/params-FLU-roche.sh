# Module specific OVERRIDE parameters here
# FLU-fast module for IRMA
# 11.20.2014, Sam Shepard

MATCH_PROC=20				# maximum processes for the MATCH per genetic unit
ALIGN_PROC=20				# maximum processes for the rough align per genetic unit
ASSEMBLE_PROC=20			# maximum processes for assembly per genetic unit
SINGLE_PREPRO_PROC=12			# maximum processes for QC on local node
DOUBLE_PREPRO_PROC=6			# maximum processes for the left or right pairs
MIN_FI=0.10				# minimum insertion variant frequency
MIN_FD=0.10				# minimum deletion variant frequency
MIN_F=0.10				# minimum frequency for variants
INCL_CHIM=0				# whether or not to include chimera
QUAL_THRESHOLD=10			# average or median threshold
MIN_LEN=50				# minimum read length for QUALITY reads
MAX_ROUNDS=2				# maximum number of iterations to BLAT
SKIP_E=0				# perform elongation
AUTO_F=0				# auto-adjust frequency threshold
MIN_RP=1				# minimum read pattern count
MIN_RC=1				# minimum read count
INS_T=0.50				# threshold for insertion refinement of references
DEL_T=0.50				# threshold for deletion refinement of references
MIN_C=2					# minimum count for variants
MIN_AQ=20				# minimum average variant quality, does not apply to deletions
MIN_TCC=30				# minimum non-ambiguous column coverage
SIG_LEVEL=0.95				# significance test level for variant calling (.90,.95,.99,.999). 


# DO NOT ALTER
BLAT_SORT=0				# sort using BLAT
NONSEGMENTED=0				# non-segmented virus
LFASTM=1					# LABEL sorting fast-mode
SSW_M=2					# smith-waterman match score
SSW_X=5					# smith-waterman mismatch penalty
SSW_O=10				# smith-waterman gap open penalty
SSW_E=1					# smith-waterman gap extension penalty
#Lmpath=$mpath				# HMM module path
