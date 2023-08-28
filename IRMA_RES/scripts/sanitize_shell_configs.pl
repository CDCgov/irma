#!/usr/bin/env perl
# Filename:         sanitize_shell_configs
#
# Description:      Sanitizes shell variable configuration input by checking for
#                   valid formats. Valid formats should be safe to source.
#
# Date dedicated:   2023-08-21
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
#
# Citation:         Shepard SS, Meno S, Bahl J, Wilson MM, Barnes J, Neuhaus E.
#                   Viral deep sequencing needs an adaptive approach: IRMA, the
#                   iterative refinement meta-assembler. BMC Genomics.
#                   2016;17(1). doi:10.1186/s12864-016-3030-6
#
# =============================================================================
#
#                            PUBLIC DOMAIN NOTICE
#
#  This source code file or script constitutes a work of the United States
#  Government and is not subject to domestic copyright protection under 17 USC §
#  105. This file is in the public domain within the United States, and
#  copyright and related rights in the work worldwide are waived through the CC0
#  1.0 Universal public domain dedication:
#  https://creativecommons.org/publicdomain/zero/1.0/
#
#  The material embodied in this software is provided to you "as-is" and without
#  warranty of any kind, express, implied or otherwise, including without
#  limitation, any warranty of fitness for a particular purpose. In no event
#  shall the Centers for Disease Control and Prevention (CDC) or the United
#  States (U.S.) government be liable to you or anyone else for any direct,
#  special, incidental, indirect or consequential damages of any kind, or any
#  damages whatsoever, including without limitation, loss of profit, loss of
#  use, savings or revenue, or the claims of third parties, whether or not CDC
#  or the U.S. government has been advised of the possibility of such loss,
#  however caused and on any theory of liability, arising out of or in
#  connection with the possession, use or performance of this software.
#
#  Please provide appropriate attribution in any work or product based on this
#  material.

use strict;
use warnings;
use Carp qw(croak);
use English qw(-no_match_vars);

my $bool      = '[01]';
my $integer   = '\d+';
my $float     = '\d|0.\d+';
my $path      = '[\sA-Za-z0-9/_.-]+';
my $date      = '[0-9./-]+';
my $free_text = '[\sa-zA-Z0-9_.-]+';

my %valid_config = (
                     ADAPTER                  => '[AGCTN]+',
                     ALIGN_AMENDED            => $bool,
                     ALIGN_PROC               => $integer,
                     ALIGN_PROG               => '(SAM|BLAT)(\s(SAM|BLAT))*',
                     ALLOW_DISK_CHECK         => $bool,
                     ALLOW_TMP                => $bool,
                     ASSEM_PROC               => $integer,
                     ASSEM_PROG               => 'SSW|MINIMAP2',
                     ASSEM_REF                => $bool,
                     AUTO_F                   => $bool,
                     BAN_GROUPS               => '\w+',
                     BLAT_IDENTITY            => $integer,
                     CUSTOM_REF_FILE          => $path,
                     DEF_SET                  => $path,
                     DEL_T                    => $float,
                     DEL_T_DEPTH              => $integer,
                     DEL_TYPE                 => '(REF|NNN|DEL)(\s(REF|NNN|DEL))*',
                     DO_SECONDARY             => $bool,
                     DOUBLE_LOCAL_PROC        => $integer,
                     ENFORCE_CLIPPED_LENGTH   => $bool,
                     FUZZY_ADAPTER            => $bool,
                     GENE_GROUP               => '[a-zA-Z0-9_,:-]+',
                     GRID_ON                  => $bool,
                     GRID_PATH                => $path,
                     INCL_CHIM                => $bool,
                     INS_T                    => $float,
                     INS_T_DEPTH              => $integer,
                     IRMA_QUEUE               => '[A-Za-z.]+',
                     LABEL_MODULE             => '[a-zA-Z0-9_-]+',
                     LIMIT_BLAT               => $integer,
                     LIMIT_LABEL              => $integer,
                     LIMIT_PHASE              => $integer,
                     LIMIT_SAM                => $integer,
                     LIMIT_SSW                => $integer,
                     MATCH_PROC               => $integer,
                     MATCH_PROG               => 'BLAT',
                     MAX_ITER_ASSEM           => $integer,
                     MAX_ITER_SSW             => $integer,
                     MAX_ROUNDS               => $integer,
                     MERGE_SECONDARY          => $bool,
                     MIN_AMBIG                => $float,
                     MIN_AQ                   => $integer,
                     MIN_BLAT_MATCH           => $integer,
                     MIN_C                    => $integer,
                     MIN_CA                   => $integer,
                     MIN_CONF                 => $float,
                     MIN_CONS_QUALITY         => $integer,
                     MIN_CONS_SUPPORT         => $integer,
                     MIN_DROPOUT_EDGE_DEPTH   => $integer,
                     MIN_F                    => $float,
                     MIN_FA                   => $float,
                     MIN_FD                   => $float,
                     MIN_FI                   => $float,
                     MIN_LEN                  => $integer,
                     MIN_RC                   => $integer,
                     MIN_RC_RESIDUAL          => $integer,
                     MIN_RP                   => $integer,
                     MIN_RP_RESIDUAL          => $integer,
                     MIN_TCC                  => $integer,
                     MM2_A                    => $integer,
                     MM2_B                    => $integer,
                     MM2_E                    => $integer,
                     MM2_O                    => $integer,
                     NO_MERGE                 => $bool,
                     NONSEGMENTED             => $bool,
                     NO_SORT_REFS             => $bool,
                     PACKAGED_FASTQ           => $bool,
                     PADDED_CONSENSUS         => $bool,
                     PARAM_FILE_AUTHOR        => $free_text,
                     PARAM_FILE_DATE          => $date,
                     PARAM_FILE_NAME          => $free_text,
                     PARAM_FILE_VERSION       => $free_text,
                     PHASE_DISTANCE           => $float,
                     PHASE_PROC               => $integer,
                     QUAL_THRESHOLD           => $integer,
                     REF_SET                  => $path . '|\$DEF_SET',
                     RESIDUAL_ASSEMBLY_FACTOR => $integer,
                     SECONDARY_LABEL_MODULES  => '[a-zA-Z0-9_.,:-]+',
                     SECONDARY_SORT           => $bool,
                     SEG_NUMBERS              => '[a-zA-Z0-9_.,:-]+',
                     SIG_LEVEL                => '0?(.90|.95|.99|.999)',
                     SILENCE_COMPLEX_INDELS   => $bool,
                     SINGLE_LOCAL_PROC        => $integer,
                     SKIP_E                   => $bool,
                     SORT_GROUPS              => '[a-zA-Z0-9/_,:-]+',
                     SORT_PROC                => $integer,
                     SORT_PROG                => '(LABEL|BLAT)(\s(LABEL|BLAT))*',
                     SSW_E                    => $integer,
                     SSW_M                    => $integer,
                     SSW_O                    => $integer,
                     SSW_X                    => $integer,
                     TMP                      => $path,
                     USE_MEDIAN               => $bool,
                     USE_IRMA_CORE            => $bool
);

open( my $CONFIG, '<', $ARGV[0] ) or die "Error: cannot open file '$ARGV[0]' for reading. See: $OS_ERROR\n";
local $RS = "\n";
my @sanitized_configs = ();
while ( my $line = <$CONFIG> ) {
    chomp($line);
    my ( $key, $value ) = map { trim($_) } ( split /[#=]/smx, $line );

    if ( !defined $key || $key eq q{} ) {
        next;
    }

    if ( !defined $valid_config{$key} ) {
        die "Error: configuration '$key' is not a valid IRMA configuration value! Please check file: '$ARGV[0]'\n";
    } else {
        my $re = $valid_config{$key};
        if ( $value eq q{} || $value =~ /^$re|"$re"|'$re'|""|''$/smx ) {
            push( @sanitized_configs, "$key=$value" );
        } else {
            die "Error: configuration '$key' expects pattern /$re/ but found: $value\n";
        }
    }
}
close $CONFIG or croak("Error: failed to close file. See: $OS_ERROR\n");

if ( scalar @sanitized_configs > 0 ) {
    print STDOUT join( ";", @sanitized_configs );
}

sub trim {
    my $string = shift // q{};
    if ( $string =~ /^\s*(.*?)\s*$/smx ) {
        return $1;
    } else {
        return $string;
    }
}

