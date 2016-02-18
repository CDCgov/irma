

cat <<EOF > $ppath/logs/FLU-$RUN.sh
### BACKGROUND INFO ###
# $PROGRAM, v$VERSION, $DATE
# $AUTHOR ($AFFIL), $EMAIL
# Last commit: $LAST_COMMIT

### START CONFIG ###
PARAM_FILE_NAME=$PARAM_FILE_NAME
PARAM_FILE_AUTHOR=$PARAM_FILE_AUTHOR
PARAM_FILE_VERSION=$PARAM_FILE_VERSION
PARAM_FILE_DATE=$PARAM_FILE_DATE

### PERFORMANCE ###
GRID_ON=$GRID_ON	# grid computation on [1,0] for on or off
LIMIT_BLAT=$LIMIT_BLAT	# threshold before grid
LIMIT_SSW=$LIMIT_SSW		# threshold before grid
LIMIT_SAM=$LIMIT_SAM		# threshold before grid
SINGLE_LOCAL_PROC=$SINGLE_LOCAL_PROC	# local maximum processes
DOUBLE_LOCAL_PROC=$DOUBLE_LOCAL_PROC	# local maximum processes (double this number)
ALLOW_TMP=$ALLOW_TMP		# if GRID_ON=0, try to use /tmp for working directory
TMP=$TMP		# the scratch/tmpfs for working on the assemblies

### REFERENCE ###
MIN_FA=$MIN_FA		# no alternative reference [0..1]
MIN_CA=$MIN_CA		# minimum count for alternative finished assembly
SKIP_E=$SKIP_E		# skip reference elongation
REF_SET=$REF1_SET	# Starting reference, usually default for \$DEF_SET
ASSEM_REF=$ASSEM_REF

### READ GATHERING ###
FASTA=$FASTA			# accept fasta format
MAX_ROUNDS=$MAX_ROUNDS		# round of read gathering
USE_MEDIAN=$USE_MEDIAN		# use the median quality or the average [1,0]
QUAL_THRESHOLD=$QUAL_THRESHOLD	# minimum read statistic
MIN_LEN=$MIN_LEN		# minimum read length
INCL_CHIM=$INCL_CHIM		# includes chimera or not [0,1]

## MATCH STEP
MATCH_PROC=$MATCH_PROC		# grid maximum processes for the MATCH
MATCH_PROG="${MATCH_PROGS[@]}"	# match (all or any match) program [BLAT]
MIN_RP=$MIN_RP		# minimum read pattern count to continue
MIN_RC=$MIN_RC		# minimum read count to continue

## SORT STEP 
SORT_PROG="${SORT_PROGS[@]}"	# [LABEL,BLAT]
SORT_PROC=$SORT_PROC		# currently not used
NONSEGMENTED=$NONSEGMENTED		# segmented! [0,1]
# LABEL
LFASTM=$LFASTM		# LABEL sorting fast-mode

## ALIGN STEP ##
ALIGN_PROG="${ALIGN_PROGS[@]}"	# rough assembly / alignment to working reference [SAM,BLAT]
ALIGN_PROC=$ALIGN_PROC		# grid maximum processes for the rough align

### FINISHING ASSEMBLY ###
MAX_ITER_ASSEM=$MAX_ITER_ASSEM	# max assembly iteration [5]
NO_MERGE=$NO_MERGE		# do not merge read pairs [0]
ASSEM_PROG="$ASSEM_PROG"	# assembly program [SSW]
ASSEM_PROC=$ASSEM_PROC		# grid maximum processes for assembly
INS_T=$INS_T		# minimum frquenncy threshold for insertion refinement
DEL_T=$DEL_T		# minimum frequency threshold for deletion refinement 
MIN_AMBIG=$MIN_AMBIG		# minimum called SNV frequency for mixed base in amended consensus folder
SSW_M=$SSW_M			# smith-waterman match score
SSW_X=$SSW_X			# smith-waterman mismatch penalty
SSW_O=$SSW_O		# smith-waterman gap open penalty
SSW_E=$SSW_E			# smith-waterman gap extension penalty

### VARIANT CALLING ###
# HEURISTICS
AUTO_F=$AUTO_F		# auto-adjust frequency threshold [1,0]
MIN_FI=$MIN_FI		# minimum insertion variant frequency
MIN_FD=$MIN_FD		# minimum deletion variant frequency
MIN_F=$MIN_F		# minimum frequency for single nucleotide variants
MIN_C=$MIN_C			# minimum count for variants
MIN_AQ=$MIN_AQ		# minimum average variant quality, does not apply to deletions
MIN_TCC=$MIN_TCC		# minimum non-ambiguous column coverage
MIN_CONF=$MIN_CONF		# minimum confidence not machine error

# CONFIDENCE INTERVALS
SIG_LEVEL=$SIG_LEVEL		# significance test level for variant calling (.90,.95,.99,.999). 
EOF
