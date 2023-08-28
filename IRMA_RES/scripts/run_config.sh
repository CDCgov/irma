cat << EOF > "$ppath"/logs/${MODULE}-$RUN.sh
### BACKGROUND INFO ###
# $PROGRAM, v$VERSION, $DATE
# $AUTHOR ($AFFIL), $EMAIL
# Last commit: $LAST_COMMIT

### START CONFIG ###
PARAM_FILE_NAME="$PARAM_FILE_NAME"
PARAM_FILE_AUTHOR="$PARAM_FILE_AUTHOR"
PARAM_FILE_VERSION="$PARAM_FILE_VERSION"
PARAM_FILE_DATE="$PARAM_FILE_DATE"

### PERFORMANCE ###
GRID_ON=$GRID_ON	# grid computation on [1,0] for on or off
GRID_PATH="${GRID_PATH}"	# optional grid path for working directory, leave "" for no option
LIMIT_BLAT=$LIMIT_BLAT	# threshold before grid
LIMIT_SSW=$LIMIT_SSW		# threshold before grid
LIMIT_SAM=$LIMIT_SAM		# threshold before grid
SINGLE_LOCAL_PROC=$SINGLE_LOCAL_PROC	# local maximum processes
DOUBLE_LOCAL_PROC=$DOUBLE_LOCAL_PROC	# local maximum processes (double this number)
ALLOW_TMP=$ALLOW_TMP		# if GRID_ON=0, try to use /tmp for working directory
TMP=$TMP		# the scratch/tmpfs for working on the assemblies
PACKAGED_FASTQ=${PACKAGED_FASTQ:-1}		# Create *.tar.gz for final fastq, otherwise use *.gz for each file
USE_IRMA_CORE=${USE_IRMA_CORE:-0}	# Use IRMA CORE where available for improved performance

### REFERENCE ###
MIN_FA=$MIN_FA		# no alternative reference [0..1]
MIN_CA=$MIN_CA		# minimum count for alternative finished assembly
SKIP_E=$SKIP_E		# skip reference elongation
REF_SET=$REF1_SET	# Starting reference, usually default for \$DEF_SET
ASSEM_REF=$ASSEM_REF

### READ GATHERING ###
MAX_ROUNDS=$MAX_ROUNDS		# round of read gathering
USE_MEDIAN=$USE_MEDIAN		# use the median quality or the average [1,0]
QUAL_THRESHOLD=$QUAL_THRESHOLD	# minimum read statistic
MIN_LEN=$MIN_LEN		# minimum read length
INCL_CHIM=$INCL_CHIM		# includes chimera or not [0,1]
			# transposase adapter, clips 5' of the adapter on the forward strand and 3' on the reverse strand
			# applicable to nexttera pair-end reads 
ADAPTER="$ADAPTER"
FUZZY_ADAPTER=${FUZZY_ADAPTER:-1}			# Allow 1 adapter mismatch
ENFORCE_CLIPPED_LENGTH=${ENFORCE_CLIPPED_LENGTH:-1}	# Reads are filtered for minimum length post adapter trimming.	
MERGE_SECONDARY=$MERGE_SECONDARY	# merge secondary data after the first round, good if no co-infections
BLAT_IDENTITY=${BLAT_IDENTITY:-80}	# blat identity for read gathering, usually 80
DO_SECONDARY=${DO_SECONDARY:-0}		# do secondary assembly
RESIDUAL_ASSEMBLY_FACTOR=${RESIDUAL_ASSEMBLY_FACTOR:-400}	# the primary match to alternative match (secondary) ratio limit. If less than this number, do the residual.

## MATCH STEP
MATCH_PROC=$MATCH_PROC		# grid maximum processes for the MATCH
MATCH_PROG="${MATCH_PROGS[@]}"	# match (all or any match) program [BLAT]
MIN_RP=$MIN_RP		# minimum read pattern count to continue doing primary assembly
MIN_RC=$MIN_RC		# minimum read count to continue doing primary assembly
MIN_RP_RESIDUAL=${MIN_RP_RESIDUAL:-150}		# minimum read pattern count to continue doing residual assembly
MIN_RC_RESIDUAL=${MIN_RC_RESIDUAL:-150}		# minimum read count to continue doing residual assembly
MIN_BLAT_MATCH=${MIN_BLAT_MATCH:-0}	# minimum blat match, default settings within the program practically limit to 30 bp, only useful if set higher.

## SORT STEP 
SORT_PROG="${SORT_PROGS[@]}"	# [LABEL,BLAT]
SORT_PROC=$SORT_PROC		# currently not used
NONSEGMENTED=$NONSEGMENTED		# segmented! [0,1]
# LABEL
LFASTM=$LFASTM		# LABEL sorting fast-mode
GENE_GROUP="$GENE_GROUP"

## ALIGN STEP ##
ALIGN_PROG="${ALIGN_PROGS[@]}"	# rough assembly / alignment to working reference [SAM,BLAT]
ALIGN_PROC=$ALIGN_PROC		# grid maximum processes for the rough align
DEL_TYPE="${DEL_TYPES[@]}"	# rough assembly deletion handling: NNN, REF, or DEL

### FINISHING ASSEMBLY ###
MAX_ITER_ASSEM=$MAX_ITER_ASSEM	# max assembly iteration [5]
NO_MERGE=$NO_MERGE		# do not merge read pairs [0]
ASSEM_PROG="$ASSEM_PROG"	# assembly program [SSW]
ASSEM_PROC=$ASSEM_PROC		# grid maximum processes for assembly
INS_T=$INS_T		# minimum frquenncy threshold for insertion refinement
DEL_T=$DEL_T		# minimum frequency threshold for deletion refinement 
INS_T_DEPTH=$INS_T_DEPTH		# minimum coverage depth for insertion refinement
DEL_T_DEPTH=$DEL_T_DEPTH		# minimum coverage depth for deletion refinement 
SILENCE_COMPLEX_INDELS=${SILENCE_COMPLEX_INDELS:-0}	# silences reads with complex indels (having 4 indels with at least one greater than 9)
MIN_AMBIG=$MIN_AMBIG		# minimum called SNV frequency for mixed base in amended consensus folder
SSW_M=$SSW_M			# smith-waterman match score
SSW_X=$SSW_X			# smith-waterman mismatch penalty
SSW_O=$SSW_O			# smith-waterman gap open penalty
SSW_E=$SSW_E			# smith-waterman gap extension penalty
MM2_A=${MM2_A:-$SSW_M}		# minimap2 match
MM2_B=${MM2_B:-$SSW_X}		# minimap2 mismatch
MM2_O=${MM2_O:-$SSW_O}		# minimap2 gap open
MM2_E=${MM2_E:-$SSW_E}		# minimap2 gap extend
MIN_CONS_SUPPORT=${MIN_CONS_SUPPORT:-1}	# minimum consensus support for amended consensus
MIN_CONS_QUALITY=${MIN_CONS_QUALITY:-0}	# minimum consensus average quality for amended consensus
ALIGN_AMENDED=${ALIGN_AMENDED:-0}		# do global alignment of the amended consensus to the HMM profile
PADDED_CONSENSUS=${PADDED_CONSENSUS:-1}	# attempt to pad amended_consensus with Ns for amplicaton dropout: requires ALIGN_AMENDED=1 and ASSEM_REF=1
MIN_DROPOUT_EDGE_DEPTH=${MIN_DROPOUT_EDGE_DEPTH:-0} # minimum threshold before dropouts (91+) with flanking regions (6+) are masked during final iterative refinement

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
