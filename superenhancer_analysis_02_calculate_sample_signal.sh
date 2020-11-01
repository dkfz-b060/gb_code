#!/bin/bash

########################################
#
# superenhancer_analysis_02_calculate_sample_signal.sh
#
# Mike Fletcher
# 20201025
#
# (original name: superenhancers_subtype_02_calculate_sample_signal.sh)
#
########################################
#
# WHAT THIS DOES
#
########################################
#
# based on Carl's script define_enhancers.sh
#
#	takes the 12.5kbp-stitched 'enhancer union' from each GBM subtype (bed file, output of script 01)
# 	uses bigWigAverageOverBed to calculate signal (for **all** samples, not just that subtype!) per-sample
#
########################################
# INPUTS
########################################
#
# for each subtype analysed: the stitched H3K27ac peak .bed file from superenhancer_analysis_01_prepare_data.Rscript
#
# for each sample: the H3K27ac signal intensities, in genome-wide .bw/.bigWig format
#		available from GEO, see https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121719
#
########################################
# OUTPUTS
########################################
#
# per-sample .tab files giving the signal per region
#
########################################
# RESOURCE USAGE
########################################
#
# bigWigAverageOverBed is very lightweight and runs quickly (<1min/sample)
#
########################################

##############################
##############################
#
# DEFINITIONS
#
##############################
# histone mark to use
mark="H3K27ac"

# data dir: view-by-pid dir for the ChIPseq data
datadir="/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_newAlignments/wgs_bowtie/"

# output dir: my analysis dir
outdir="/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_superenhancers/"

##############################
#
# START CODE
#
##############################
# work in my output dir:
cd ${outdir}

# make new output dir for the per-sample signal files, if doesn't already exist
if [ ! -d ${outdir}/sample_signal/ ]; then
    mkdir ${outdir}/sample_signal/
fi

echo "Scoring with bigwig ..."
for st in MES IDH RTK_I RTK_II; do 
    echo "=== ${st} ==="
	
	# get the 12.5kbp stitched per-sample enhancer union; take first 3 columns (chrom/start/end) and 8th (peak ID), then remove header line
    cut -f1-3,8 ${outdir}/${st}_subtype_enhancer_union_12500_stitched.bed | tail -n +2 > lite
	
	# now: for all sample IDs (found by looking in view-by-pid dir and taking anything matching AK*)
	# calculate the per-sample signal for the subtype
    for i in `ls ${datadir} | grep ^AK`; do
	echo "       $i"
	if [ -f ${datadir}/${i}/${mark}/bw/${i}_${mark}_SES_subtract.bw ]; then
	    /ibios/tbi_cluster/11.4/x86_64/bin/bigWigAverageOverBed ${datadir}/${i}/${mark}/bw/${i}_${mark}_SES_subtract.bw lite ${outdir}/sample_signal/${i}_${mark}_in_${st}_enhancer_union.tab
	fi
    done
	
	# clean up tmp files
	rm lite
done