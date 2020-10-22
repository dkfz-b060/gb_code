#!/bin/R-3.4.3

# gbm_rtn_analysis_common.R
#
# 20170124
#
# Common code for the RTN analysis on GBM data:
# - definitions for the analysis (e.g. paths)
# - functions that are reused throughout the code

####################################################################################################
#
# FIRST: load common HIPO016 analysis code...
#
source("/home/fletcher/git_repos/gbm-master-regulators/gbm_hipo016_common.R")
#
####################################################################################################

####################################################################################################
#
# DEFINITIONS
# paths, options, etc., for use throughout the analysis
#
####################################################################################################
message( "\nLoading common code for RTN analysis...\n")
###############################################
#
# PATHS to resources, etc.
#
###############################################
# path to working directory for output
# must be already created...
working.dir <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/rtn/results_gbm_rerun_AK088/"
# set working dir to this
setwd(working.dir)

###############
# input array GX data (raw Affy array .cel files) for TNI
###############
# path to data directory containing input .CEL files, in view-by-sample id format - as in samplesheet
data.dir.test <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/microarray_tcga_HT_HG-U133A/view_by_sampleid/"
# path to data directory containing input .CEL files, in view-by-sample id format - as in samplesheet
data.dir.validation <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/microarray_assorted_HG-U133A-Plus2/cel_files/"

###############
# samplesheets for input array GX data for TNI
###############
# path to TCGA samplesheet = MAGE-TAB table
# it's up one dir from the actual view-by-sampleID dir
samplesheet.path.test <- paste( data.dir.test, "../tcga_gbm_HT_HG-U133A_samplesheet.csv" , sep="")
# path to samplesheet = MAGE-TAB table
# it's up one dir from the actual view-by-sampleID dir
samplesheet.path.validation <- paste( data.dir.validation, "../validation_GBMs_samplesheet_final.csv" , sep="")

###############
# TF lists for TNI 
###############
# list of TFs from Vaquerizas et al
vaquerizas.tfs.path <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/rtn/vaquerizas_tfs_nrg2538-s3.txt"

###############
# TCGA GX signature analysis
###############
# path to the d-statistics from the TCGA paper, in .csv
# see above for details
tcga.d.stats.path <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/rtn/tcga_gx_SAM_pairwise.csv"

# path to dir containing TCGA GX sigs - used in GSEA analysis
# these are header-less, single column files of gene symbols, named:
# TCGA_GBM_GX_signature_CL.txt   TCGA_GBM_GX_signature_NL.txt
# TCGA_GBM_GX_signature_MES.txt  TCGA_GBM_GX_signature_PN.txt
tcga.sigs.path <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/gsea/"

###############################################
#
# PARALLEL PROCESSING RESOURCES
#
###############################################
# #cores to use in parallel analyses - used in tni permutation
parallel.cores.tni <- 32
# #cores to use in parallel analyses - used in tna analyses
# 20190805 both 32 and 16 cores causes crashes (ffs)
parallel.cores.tna <- 12

###############################################
#
# GENE EXPRESSION DATA
# load and process:
# 	- limma Renv for subtype genes
#	- Zuguang's GX matrices for GX
#
###############################################
#
# limma Renv
#	loads many R objects from Renv
#
########################
# load limma Renv object
load( limma.renv )

########################
#
# RNAseq matrices from Zuguang
#	final R object name: 'expr'
#	with rows = genes (ENSG IDs), columns = sample AK IDs
#
########################
# load RNAseq matrix:
load(path.rnaseq.matrix)
# get rpkm
expr <- rpkm

####################################################################################################
#
#
# FUNCTIONS
#
#
####################################################################################################

####################################################################################################
#
# process.tni()
#	produces a processed TNI object (using parallel processing)
#
# options:
#		gx.matrix = filtered gene expression matrix, rows = genes X columns = samples
#		array.anno = annotation array containing probeIDs and gene IDs
#		tfs = list of TFs to use in the analysis
#		parallel.cores.tni = #cores to use in the analysis (defined at start of the script)
#
####################################################################################################
process.tni <- function( gx.matrix, array.anno, tfs, parallel.cores.tni )
{
	####################
	# first: create TNI from GX data
	####################
	# create new TNI and pre-process!:
	# for TFs, use those probeIDs which have gene symbols matching those present in the TF list
	# coerce to character and take unique hits only
	rtni <- new( "TNI", gexp=gx.matrix, transcriptionFactors=as.character( unique( array.anno$probe_id[ array.anno$symbol %in% tfs ] ) ) )
	# then preprocess
	# use as gexpIDs the constructed array annotation df
	rtni <- tni.preprocess( rtni, gexpIDs=array.anno )

	#########################
	# run permutation analysis, infer MI
	#########################
	# use SNOW to parallelise
	suppressMessages( library(snow) )
	# start cluster, run permutation, then stop when finished
	# use #cores as defined above
	options(cluster=makeCluster(parallel.cores.tni, "SOCK"))
	rtni <- tni.permutation(rtni, parChunks=50)
	stopCluster(getOption("cluster"))

	# run bootstrap analysis
	rtni <- tni.bootstrap(rtni)

	# apply DPI filter
	rtni <- tni.dpi.filter(rtni)
	
	# return processed TNI obj
	return(rtni)
}
####################################################################################################
 
####################################################################################################
#
# analyse.tna()	
#	function to run tna analysis + GSEA + output GSEA plots for regulons
# 
# args:
#  rtni = pre-processed rtni object
#  cohort = character string for labelling output - test or validation!
#  subtype = character string of subtype (used in output file names)
#  phenotype = named, numeric vector of phenotypes used in GSEA
#     e.g. the calculated logFC of a particular GXsig
#     use full "universe" of genes here
#  hits = character vector of geneIDs associated with a particular phenotype
#     e.g. DEG from a particular analytical contrast
#     use only the "significant" genes here
#  pval = number: adjusted pvalue cutoff for GSEA output
#     default = 0.01
#  gsea.nperm = number: #permutations to use to calculate pvals in GSEA
#     default = 10000 (will take some time to run!!)
#  parallel.cores.tna = number: #cores to parallelise over
#	  default = as defined in variables at start of script
#  gsea.stepfilter = logical T/F : whether to filter GSEA input MRs based on MRA results
#	  default = TRUE
#  gsea.tfs = NULL or character vector with gene IDs (=symbols): perform 1-/2-GSEA on the specified TFs
#	  default = NULL
#  gsea.minRegulonSize = number: #genes in regulon required before testing in GSEA
#	  default = 15
#
# outputs:
#    the processed TNA object, with MRA, 1- and 2-tail GSEA, overlap, synergy and shadow analyses all completed.
#    also, in working directory: PDFs of 1- and 2-tailed GSEA results.
#
analyse.tna <- function( rtni, cohort, subtype, phenotype, hits, pval=0.01, gsea.nperm=10000, parallel.cores.tna, mra.pval=0.05, gsea.tfs=NULL, gsea.stepfilter=TRUE, gsea.minRegulonSize=15 )
{
    # clean up stuff, for mem usage
	gc()
	
	message( paste("Running TNA analysis for subtype:", subtype, "in cohort:", cohort, sep=" ") )
	#
	# define and create per-subtype output directory (for many output plots)
	subtype.dir <- paste( "tna_results", subtype, sep="_" )
	# make new working dir for today
	system2("mkdir", args=subtype.dir)
	# pre-process tni -> tna
    #
    # options:
    # object = tni obj
    # phenotype = numeric vector of GX values (with names), ranked in the phenotype of interest in MRA
    #    i.e. the universe of genes measured in RTK_I ranked by some statistic - logFC
    #       used in gsea (quantitative) methods
    #
    # hits = char vector of gene IDs considered as 'hits'
    #    i.e. 
    #       used in MRA (hypergeometric/overlap-based) methods.
    #
    # n.b. hits and phenotype are different representations of the same phenotype!
    #
    # phenoIDs = 2-column df with col1 = probe IDs, col2 = gene IDs.
    rtna <- tni2tna.preprocess(object=rtni, phenotype=phenotype, hits=hits )
    #####
    # run MRA
    rtna <- tna.mra(rtna, pValueCutoff=mra.pval)
	#####
	# start SNOW cluster to speed up following analyses:
	# tna.gsea1, tna.gsea2, tna.synergy, tna.shadow
	message( "Start SNOW cluster with ", parallel.cores.tna, " cores" )
	options(cluster=makeCluster(parallel.cores.tna, "SOCK"))
	#####
	###############  
    # run GSEA analysis
    #
    # run both one- and two-tailed analyses
    # set higher pvalcutoff
    # increase #permutations in attempt to calculate exact pval
	#####
	# run 1-tail analysis
	#####
    rtna <- tna.gsea1( rtna, pValueCutoff=pval, minRegulonSize=gsea.minRegulonSize, nPermutations=gsea.nperm, tfs=gsea.tfs, stepFilter=gsea.stepfilter )
    # plot GSEA results
    outputname <- paste( subtype.dir, "/", Sys.Date(), "_RTN_GBM_", subtype, "_", cohort, "_GSEA_1-tail", sep="" )
    # options:
    # labPheno = phenotype label
    labpheno <- paste( subtype, " - ", cohort, " - logFC", sep="")
    # order regulons by adj. pval
    #
    # plot 1-tail
    tna.plot.gsea1(rtna, labPheno=labpheno, file=outputname, 
            regulon.order="adj.pvalue", width=6, height=15, heightPanels=c(1,4,4), ylimPanels=c(0,3.5,0,1))
			
    # clean up stuff, for mem usage
	gc()		

    #####
	# do for 2-tail
	#####
    rtna <- tna.gsea2( rtna, pValueCutoff=pval, minRegulonSize=gsea.minRegulonSize, nPermutations=gsea.nperm, tfs=gsea.tfs, stepFilter=gsea.stepfilter )
    #
    # plot 2-tail
    outputname <- paste( subtype.dir, "/", Sys.Date(), "_RTN_GBM_", subtype, "_", cohort, "_GSEA_2-tail", sep="" )
    tna.plot.gsea2(rtna, labPheno=labpheno, file=outputname, 
            regulon.order="adj.pvalue", width=6, height=6, heightPanels=c(1.5,1,4), ylimPanels=c(1,1.5,0,1))
			
	# clean up stuff, for mem usage
	gc()
	
	#####
    # run overlap, synergy and shadow analyses
    #
    # use pval cutoff of 0.001 (as per our paper)
    # use min regulon size of 20 (as per our paper)
    #
    # n.b. default option for synergy/shadow is to use stepFilter=T - removes non-sig regulons based on GSEA analysis
    # default nPermutations = 1000
    #rtna <- tna.overlap( rtna, pValueCutoff = 0.001, minRegulonSize = 20 )
    #rtna <- tna.synergy( rtna, pValueCutoff = 0.001, minRegulonSize = 20 )
    #rtna <- tna.shadow( rtna, pValueCutoff = 0.001, minRegulonSize = 20 )
	#####
	# stop cluster
	message( "Stop SNOW cluster" )
	stopCluster(getOption("cluster"))
	
    # clean up stuff, for mem usage
	gc()
	
    #####
    # also return processed tna obj
    return(rtna)
	#####
}
####################################################################################################

####################################################################################################
#
# function to do TNA analysis using the TCGA d-statistics for each SUBTYPEvsALL comparison
#
# tweaks to the analyse.tna above are:
#
# 1. don't have pvals, so have no 'hits'; so don't run tna.mra
# 2. will generate the 'phenotype' vector for each subtype within function
# 3. as don't have pvals, need to set stepFilter=FALSE in gsea funs...
# 4. adjust pval and minregulon thresholds for GSEA functions as we have too many MRs otherwise...
#
# 
analyse.tna.tcga <- function( rtni, cohort, subtype, pval1tail=0.01, pval2tail=0.01, gsea.nperm=10000, minreg=20, parallel.cores.tna, gsea.tfs=NULL, gsea.stepfilter=TRUE, gsea.minRegulonSize=15 )
{
    message( paste("Running TNA analyses for TCGA subtype: ", subtype) )

	# define and create per-subtype output directory (for many output plots)
	subtype.dir <- paste( "tna_results_TCGA", subtype, sep="_" )
	# make new working dir for today
	system2("mkdir", args=subtype.dir)

	# first: generate the 'phenotype' vector for GSEA
	# 	
	# get the name of the column name needed 
	# this is in 'SUBTYPE'vsAll format (no spaces), with SUBTYPE in caps!
	colname <- paste( toupper(subtype), "vsAll", sep="")
	# extract this column from the d-stats dataframe
	phenotype <- tcga.d.stats[,colname]
	# add names! needed for tna analysis
	names(phenotype) <- tcga.d.stats$Gene
	# order d-statistics
	phenotype <- phenotype[ order(phenotype, decreasing=T) ]
	#
	# get the 'hits' (= vector of gene symbols) from the loaded TCGA sigs
	# use the 1st column
	hits.name <- paste( "tcga.sig.", toupper(subtype), sep="" )
	hits <- get(hits.name)[,1]
	#
    # pre-process tni -> tna
    #
    # options:
    # object = tni obj
    # phenotype = numeric vector of GX values (with names), ranked in the phenotype of interest in MRA
    #    i.e. the universe of genes measured in RTK_I ranked by some statistic - logFC
    #       used in gsea (quantitative) methods
	# hits = vector of gene symbols for use in MRA
    #
    # n.b. hits and phenotype are different representations of the same phenotype!
    #
    # phenoIDs = 2-column df with col1 = probe IDs, col2 = gene IDs.
    rtna <- tni2tna.preprocess( object=rtni, phenotype=phenotype, hits=hits )
    #####
    # run MRA
    rtna <- tna.mra(rtna)
	#####
	# start SNOW cluster to speed up following analyses:
	# tna.gsea1, tna.gsea2, tna.synergy, tna.shadow
	message( "Start SNOW cluster with ", parallel.cores.tna, " cores" )
	options(cluster=makeCluster(parallel.cores.tna, "SOCK"))
	#####
	###############  
    # run GSEA analysis
    #
    # run both one- and two-tailed analyses
    # set higher pvalcutoff
    # increase #permutations in attempt to calculate exact pval
	#####
	# run 1-tail analysis
	#####
    rtna <- tna.gsea1( rtna, pValueCutoff=pval1tail, nPermutations=gsea.nperm, minRegulonSize=gsea.minRegulonSize, tfs=gsea.tfs, stepFilter=gsea.stepfilter )
    # plot GSEA results
    outputname <- paste( subtype.dir, "/", Sys.Date(), "_RTN_GBM_TCGAsig_", subtype, "_", cohort, "_GSEA_1-tail", sep="" )
    # options:
    # labPheno = phenotype label
    labpheno <- paste( subtype.dir, "/", subtype, " (TCGA) - (d-stat)", sep="")
    # order regulons by adj. pval
    #
    # plot 1-tail
    tna.plot.gsea1(rtna, labPheno=labpheno, file=outputname, 
            regulon.order="adj.pvalue", width=6, height=15, heightPanels=c(1,4,4), ylimPanels=c(0,3.5,0,1))

    #####
	# do for 2-tail
	#####
    rtna <- tna.gsea2( rtna, pValueCutoff=pval2tail, nPermutations=gsea.nperm, minRegulonSize=gsea.minRegulonSize, tfs=gsea.tfs, stepFilter=gsea.stepfilter )
    #
    # plot 2-tail
    outputname <- paste( subtype.dir, "/", Sys.Date(), "_RTN_GBM_TCGAsig_", subtype, "_", cohort, "_GSEA_2-tail", sep="" )
    tna.plot.gsea2(rtna, labPheno=labpheno, file=outputname, 
            regulon.order="adj.pvalue", width=6, height=8, heightPanels=c(1,4,4), ylimPanels=c(0,3.5,0,1))

	#####
	# stop cluster
	message( "Stop SNOW cluster" )
	stopCluster(getOption("cluster"))
    #####
    # also return processed tna obj
    return(rtna)
	#####
}
####################################################################################################

####################################################################################################
#
# output.tna.maps()
#	 get list of TFs to plot from GSEA 1-tail results
#	 then use this to get a graph of subtype-sig MRs in either
#   	(i) amap = association map = MRs only
#    	(ii) rmap = regulon map = MRs + their targets
#
# options: 
# rtna = rtna object for input
# subtype = character, subtype for output filenames
# cohort = cohort label - test or validation
output.tna.maps <- function( rtna, subtype, cohort )
{
	# define and per-subtype output directory (as per tna.analyse)
	subtype.dir <- paste( "tna_results", subtype, sep="_" )
	
	# get list of significant MRs from the GSEA1 results
    res.tfs <- rtna@results$GSEA2.results$differential$Regulon

    # get regulon map (gtype="rmap") using dpi-filt network (tnet="dpi")
    rtna.g.regulon <- tna.graph( rtna, tnet="dpi", gtype="rmap", tfs=res.tfs )
    # get association map (gtype="amap"), dpi-filt, filtered by quantile (amapFilter="quantile" - default), default value 0.75
    rtna.g.assoc <- tna.graph( rtna, tnet="dpi", gtype="amap", amapFilter="quantile", tfs=res.tfs )
    # output graphs as Rdata objs
    outputname <- paste( subtype.dir, "/", Sys.Date(), "_RTN_GBM_", subtype, "_TNA_", cohort, "_rmap_GSEA2.Rdata", sep="" )
    save( rtna.g.regulon, file=outputname )
    outputname <- paste( subtype.dir, "/", Sys.Date(), "_RTN_GBM_", subtype, "_TNA_", cohort, "_amap_GSEA2.Rdata", sep="" )
    save( rtna.g.assoc, file=outputname )
}
####################################################################################################

####################################################################################################
# 
# get.consensus.mrs()
#	for any pair of TNA objects, find consensus 2-tail GSEA MRs 
#		(i.e. significant in both TNA analyses, in same dES direction) 
#
# options:
# tna1, tna2 = two processed TNA objects with 2-tail GSEA analyses
# adjpval = adjusted p-value (2-tail GSEA) threshold to filter on; default = 0.01 
#	(to match TNA analysis above)
#
####################################################################################################
get.consensus.mrs <- function( tna1, tna2, adjpval=0.01 )
{
	# get significant regulon names for each cohort
	# from @results$GSEA2.results$differential$Regulon
	#
	# get only for the pval stated, default 0.01:
	out <- intersect( tna1@results$GSEA2.results$differential$Regulon[tna1@results$GSEA2.results$differential$Adjusted.Pvalue < adjpval], 
		tna2@results$GSEA2.results$differential$Regulon[tna2@results$GSEA2.results$differential$Adjusted.Pvalue < adjpval] )
		
	# next: look up dES for each sig regulon
	out.1.des <- tna1@results$GSEA2.results$differential$Observed.Score[ match(out, tna1@results$GSEA2.results$differential$Regulon) ]
	out.2.des <- tna2@results$GSEA2.results$differential$Observed.Score[ match(out, tna2@results$GSEA2.results$differential$Regulon) ]
	
	# check if dES have same sign; if they do, keep the regulon
	out <- out[ sign(out.1.des) == sign(out.2.des) ]
	
	# return regulons
	return(out)
}
####################################################################################################

####################################################################################################
#
# function to calculate Jaccard Coefficients on transcriptional network matrices ()
# stolen from code in the 'supplements.R' of Fletcher2013b 
#
##-----------------------------------------------------------------------------
#compute jaccard coefficient on tnets (among regulons)
getJC <- function(tnet)
{
  bin<-tnet
  bin[bin!=0]=1
  jc<-function(x,xmat){
    c<-x+xmat
    a<-colSums(c==2)
    b<-colSums(c>0)
    b[b==0]=1
    a/b
  }
  jcmat<-apply(bin,2,jc,xmat=bin)
  colnames(jcmat)<-colnames(bin)
  diag(jcmat)=0
  jcmat
}
####################################################################################################

####################################################################################################
# 
# heatmap.consensus.mrs()
#	generate heatmaps showing, for a set of MRs,
#		1) 2-tail GSEA dES
#		2) 2-tail GSEA significance
#		3) GX in HIPO016 RNAseq data
#
# options:
#	gsea.res = numeric matrix of differential enrichment scores, rows=MRs, columns=subtype/cohort
#	gsea.res.cellnotes = character matrix of "significant"/"non-significant", rows=MRs, columns=subtype/cohort
#	mr.gx.summary = numeric matrix of GX values, rows=MRs, columns=subtype/cohort
#	outputprefix = filename for output (will append .pdf and .png)
#
####################################################################################################
heatmap.consensus.mrs <- function( gsea.res, gsea.res.cellnotes, mr.gx.summary, outputprefix )
{
	# colours for diffES heatmap
	cols <- colorRamp2( breaks=c(-3,0,3), colors=c("blue", "white", "red") )
	# colours for significance heatmap
	sig.cols <- structure( c("white", "red"), names = c("non-significant", "significant"))
	
	# generate Heatmap object	
	p <- Heatmap( gsea.res, col=cols, name="2-tail GSEA differential ES", width=unit(6, "cm"), cluster_columns=FALSE, cluster_rows=TRUE,
		column_title="2-tail GSEA\ndifferential enrichment score", show_row_names=FALSE, column_title_gp = gpar(fontsize = 10) ) +
	# significance heatmap
	Heatmap( gsea.res.cellnotes, col=sig.cols, name="MR significance in subtype/cohort", width=unit(6, "cm"), 
		column_title="MR significance\nin subtype/cohort", show_row_names=FALSE, column_title_gp = gpar(fontsize = 10) ) +
	# MR GX heatmap
	Heatmap( as.matrix(mr.gx.summary) , col=c("white", "orange"), name="MR GX in HIPO016\nlog2(rpkm+1)", width=unit(2.5, "cm"),
		column_title="MR GX\n(HIPO016)", column_title_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 6), cluster_columns=FALSE )	
	# output as PDF: 
	# generate outputname
	outputname <- paste( outputprefix, "pdf", sep="." )
	# open plotting device
	pdf(file=outputname, width=10, height=12)
	print(p)
	dev.off()
	# output as PNG: 
	# set #ppi
	ppi <- 300
	# generate outputname
	outputname <- paste( outputprefix, "png", sep="." )
	# open plotting device
	png(file=outputname, width=10*ppi, height=12*ppi, res=ppi)
	print(p)
	dev.off()
}
####################################################################################################