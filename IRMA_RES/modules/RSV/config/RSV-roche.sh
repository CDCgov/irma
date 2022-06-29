# HEADER
PARAM_FILE_NAME="RSV-roche"
PARAM_FILE_AUTHOR="X. Qiu"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2016-08-12"

# VARIANT CALLING HEURISTICS & STATS
AUTO_F=0				# auto-adjust frequency threshold
MIN_FI=0.10				# minimum insertion variant frequency
MIN_FD=0.10				# minimum deletion variant frequency
MIN_F=0.10				# minimum frequency for variants
MIN_AQ=20				# minimum average variant quality, does not apply to deletions
MIN_TCC=30				# minimum non-ambiguous column coverage
SIG_LEVEL=0.95				# significance test level for variant calling (.90,.95,.99,.999). 

# CONSENSUS REFINEMENT & READ SELECTION
QUAL_THRESHOLD=10			# average or median threshold
MIN_LEN=50				# minimum read length for QUALITY reads
MIN_RP=1				# minimum read pattern count
MIN_RC=1				# minimum read count
INS_T=0.50				# threshold for insertion refinement of references
DEL_T=0.50				# threshold for deletion refinement of references
