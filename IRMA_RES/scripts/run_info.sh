[ "$USE_MEDIAN" == "-M" ] && USE_MEDIAN=1 || USE_MEDIAN=0
if [ "$MATCH_PROG" == "BLAT" ];then
	MATCH_OPTS=" ($BLAT_OPTS)"
else
	MATCH_OPTS=""	
fi

if  [ "$SORT_PROG" == "LABEL" ];then
	SORT_OPTS=" (LFASTM=$LFASTM)"
elif [ "$SORT_PROG" == "BLAT" ];then
	SORT_OPTS=""
else
	SORT_OPTS=""
fi

if [ "$ALIGN_PROG" == "SAM" ];then
	ALIGN_OPTS=""
elif [ "$ALIGN_PROG" == "BLAT" ];then
	ALIGN_OPTS=""
else
	ALIGN_OPTS=""
fi

if [ "$ASSEM_PROG" == "SSW" ];then
	ASSEM_OPTS=" ($SSW_OPTS)"
else
	ASSEM_OPTS=""
fi

if [ "$PAIRED" -eq "1" ];then
	INPUT_READ_DATA="LEFT;RIGHT	$LEFT;$RIGHT"
else
	INPUT_READ_DATA="LEFT	$LEFT"
fi

cat <<EOF > $ppath/logs/run_info.txt
program_name	PROGRAM	$PROGRAM
program_version	VERSION	$VERSION
last_git_commit_hash	LAST_COMMIT	$LAST_COMMIT
sample	RUN	$RUN
paired_end_reads	PAIRED	$PAIRED
data	$INPUT_READ_DATA
module_name	MODULE	$MODULE
module_param_call	MODULE_CONFIG $MODULE_CONFIG
parameter_file_name	PARAM_FILE_NAME	$PARAM_FILE_NAME
parameter_file_version	PARAM_FILE_VERSION	$PARAM_FILE_VERSION
match_last_program	MATCH_PROG	$MATCH_PROG$MATCH_OPTS
align_last_program	ALIGN_PROG	$ALIGN_PROG$ALIGN_OPTS
sort_last_program	SORT_PROG	$SORT_PROG$SORT_OPTS
assemble_last_program	ASSEM_PROG	$ASSEM_PROG$ASSEM_OPTS
match_programs	MATCH_PROGS	${MATCH_PROGS[@]}
sort_programs	SORT_PROGS	${SORT_PROGS[@]}
align_programs	ALIGN_PROGS	${ALIGN_PROGS[@]}
profiles	phmms	$phmms
do_not_merge_read_pairs	NO_MERGE	$NO_MERGE
starting_reference	REF_SET	$REF1_SET
grid_acceleration_on	GRID_ON	$GRID_ON
working_directory	ppath	$ppath
alternative_frequency	MIN_FA	$MIN_FA
alternative_count	MIN_CA	$MIN_CA
skip_reference_elongation	SKIP_E	$SKIP_E
read_statistic	USE_MEDIAN	$USE_MEDIAN
minimum_statistic_quality	QUAL_THRESHOLD	$QUAL_THRESHOLD
minimum_read_length	MIN_LEN	$MIN_LEN
interval_significance_level	SIG_LEVEL	$SIG_LEVEL
maximum_read_gathering_rounds	MAX_ROUNDS	$MAX_ROUNDS
read_fasta_input	FASTA	$FASTA
minimum_read_patterns	MIN_RP	$MIN_RP
minimum_reads	MIN_RC	$MIN_RC
minimum_insertion_edit_threshold	INS_T	$INS_T
minimum_deletion_edit_threshold	DEL_T	$DEL_T
minimum_frequency_for_mixed_consensus_call	MIN_AMBIG	$MIN_AMBIG
auto_adjust_minimum_variant_frequency	AUTO_F	$AUTO_F
minimum_deletion_variant_frequency	MIN_FD	$MIN_FD
minimum_insertion_variant_frequency	MIN_FI	$MIN_FI
minimum_variant_frequency	MIN_F	$MIN_F
minimum_variant_count	MIN_C	$MIN_C
minimum_variant_average_quality	MIN_AQ	$MIN_AQ
minimum_variant_depth	MIN_TCC	$MIN_TCC
minimum_confidence_not_sequencer_error	MIN_CONF	$MIN_CONF
EOF
