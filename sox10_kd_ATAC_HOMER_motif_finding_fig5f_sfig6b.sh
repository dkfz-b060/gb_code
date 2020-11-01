#!/bin/bash

########################################
#
# sox10_kd_ATAC_HOMER_motif_finding_fig4f_sfig6b.sh
#
# Mike Fletcher
# 20201101
#
# (original name: sox10_kd_ATAC_HOMER_motif_finding.sh; code adapted from ipynb 20191108)
#
########################################
#
# WHAT THIS DOES
#
########################################
#
# takes the differential ATAC peaks for the two cell lines (LN229, ZH487)
# between the two conditions (shSOX10, NT/control)
#
# does motif finding (HOMER) on each peakset, using the other condition as background
#
########################################
# INPUTS
########################################
#
# need, for both the LN229 and ZH487 cell lines:
#	
# - the diffBind up/down shSOX10/NT peaks, output by sox10_kd_ATAC_ChIPseq_analysis_fig4.Rscript
#	example filenames (note datestamped!):
#		2019-05-07_LN229_ATAC_diff_bound_peaks_NT_up_forHOMER.bed
#		2019-05-07_LN229_ATAC_diff_bound_peaks_shSOX10_up_forHOMER.bed
#		
#
########################################
# OUTPUTS
########################################
#
# HOMER motif finding (matching against known motifs - no de novo!) results for the up/down peaks used as input
#
########################################
# RESOURCE USAGE
########################################
#
# this does a lot of motif finding using HOMER; so to speed it up parallelise over CPUs
#		(it's quite lightweight mem-wise so I use 1g per core)
#
########################################

########################################
#
# START CODE
#
########################################
# define HOMER dir: use Chris'
HOMER_DIR="/icgc/dkfzlsdf/analysis/B060/aichmueller_pan_cancer/software/HOMER/bin/"
# add to PATH
PATH=${HOMER_DIR}:${PATH}
    
# go to analysis dir
cd /icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_sox10_kd/

###############
#
# do for ATAC diff bound peaksets, per cell line
# 
# so here: we want to use one condition as the foreground and the other as the background
# which will give fold changes interpretable as "up in condition A or B"
#
#
###############
# do for LN229
#
# generate output dir - use date, cell line, analysis name.
OUT_DIR="20191108_HOMER_results_LN229_ATAC_diff_bound_for_Volcano_plot_fg_NT"
# define input fg/bg files
fg="2019-05-07_LN229_ATAC_diff_bound_peaks_NT_up_forHOMER.bed"
bg="2019-05-07_LN229_ATAC_diff_bound_peaks_shSOX10_up_forHOMER.bed"
# run HOMER
# use as input the fg bed (ie NT up)
# genome = hg19
# output dir as generated above
# parallelise with 16 threads each
# use as background regions: bg bed (ie shSOX10 up)
# do not do de novo finding
#
# use default peak size, which will be 200bp around the peak summits given in the beds
${HOMER_DIR}/findMotifsGenome.pl ${fg} hg19 ./${OUT_DIR} -p 16 -bg ${bg} -nomotif
#
# do reverse: 
OUT_DIR="20191108_HOMER_results_LN229_ATAC_diff_bound_for_Volcano_plot_fg_shSOX10"
${HOMER_DIR}/findMotifsGenome.pl ${bg} hg19 ./${OUT_DIR} -p 16 -bg ${fg} -nomotif

# do for ZH487
# as above:
OUT_DIR="20191108_HOMER_results_ZH487_ATAC_diff_bound_for_Volcano_plot_fg_NT"
fg="2019-05-07_ZH487_ATAC_diff_bound_peaks_NT_up_forHOMER.bed"
bg="2019-05-07_ZH487_ATAC_diff_bound_peaks_shSOX10_up_forHOMER.bed"
# run HOMER
# use as input the fg bed (ie NT up)
# genome = hg19
# output dir as generated above
# parallelise with 16 threads each
# use as background regions: bg bed (ie shSOX10 up)
# do not do de novo finding
#
# use default peak size, which will be 200bp around the peak summits given in the beds
${HOMER_DIR}/findMotifsGenome.pl ${fg} hg19 ./${OUT_DIR} -p 16 -bg ${bg} -nomotif
#
# do reverse: 
OUT_DIR="20191108_HOMER_results_ZH487_ATAC_diff_bound_for_Volcano_plot_fg_shSOX10"
${HOMER_DIR}/findMotifsGenome.pl ${bg} hg19 ./${OUT_DIR} -p 16 -bg ${fg} -nomotif
