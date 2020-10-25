#!/bin/R

########################################
#
# gbm_hipo016_common.R
#
# Mike Fletcher
# 20201025
#
# R version used: both 3.4.3 and 3.5.1
#
# (original name: gbm_hipo016_common.R)
#
########################################
#
# WHAT THIS DOES
#
########################################
#
# Common R code for the HIPO016 GBM analysis:
# - definitions for the analysis (e.g. paths)
# - functions that are reused throughout the code
#
########################################

####################################################################################################
#
# DEFINITIONS
# paths, options, etc., for use throughout the analysis
#
####################################################################################################
message( "\nLoading common code for HIPO016 analysis...\n")

# set library paths: use personal lib first, then R-bundle, then R
.libPaths(new=.libPaths()[c(2,1,3)])

################
# data:
################
# location of metadata df, in RDS format:
path.metadata.rds <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/resources_mike/HIPO016_GBM_subtyping_20190804.rds"

# location of Zuguang's processed RNAseq matrices
# contains Gencode v19 genes as rows, samples as columns
path.rnaseq.matrix <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/expression/hipo16_rnaseq_count_rpkm.RData"

# location of limma analysis dir (non-datestamped, created by subtype_genes_limma.Rscript)
path.limma.analysis <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/RNAseq_subtype_genes/"

# define path to the by-sample-id basedir for the ChIPseq data:
path.chipseq <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_newAlignments/wgs_bowtie/" 

# path to limma Renv
limma.renv <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/RNAseq_subtype_genes/R_objects/2019-08-04_pipeline_end_Renv.Rdata"

# path to (minCov2) subtype enhancer dir with per-subtype .bed files:
path.enh.dir <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/cluster_enh/subtypeEnhancers/output/"

# path to Gencode v19 transcriptome in GR, saved as Rdata file
path.gencode <- "/icgc/dkfzlsdf/analysis/hipo/hipo_016/resources_mike/transcriptome-refs/gencode.v19.annotation.Rdata"

################
# plotting options:
################
# set #ppi for output PNGs
ppi <- 300

# define colours for subtypes, as named vector:
subtype.colours <- c(IDH="#e41a1c", MES="#377eb8", RTK_I="#4daf4a", RTK_II="#984ea3")

####################################################################################################
#
#
# FUNCTIONS
#
#
####################################################################################################

####################################################################################################
#
# my.mds.plot()
# 	SINGLE FUN TO PLOT MDS FROM SOME SORT OF EXPRESSION MATRIX
# 	plot MDS of data - from limma/voom analysis pipeline
#
# options:
# 		expr = some sort of expression matrix with rows=genes, cols=samples
# 		samp_id = vector of sample ids in same order as columns of expr
# 		samp_type = vector of sample subtypes in same order as columns of expr
# 		outname = name of output plot file
# 		title = title for output MDS plot
#
####################################################################################################
my.mds.plot <- function(expr, samp_id, samp_type, outname, title) 
{

    # calculate and plot Euclidean sample distances
    # calculate Euclidean distances
    sampleDists <- dist( t( expr ) )

    # coerce distances to matrix
    sampleDistMatrix <- as.matrix( sampleDists )
    # define rownames
    rownames(sampleDistMatrix) <- paste( samp_id, samp_type, sep="-" )

    # generate mds data frame:
    mds <- data.frame(cmdscale(sampleDistMatrix))
    mds$class <- samp_type
    
    # plot as png
    png(outname)
    print( qplot(X1,X2,color=class, data=mds, main=title) + scale_colour_manual(values=subtype.colours) + theme_bw() )
    dev.off()
}
####################################################################################################

####################################################################################################
#
# my.labelled.mds.plot()
# 	SINGLE FUN TO PLOT MDS FROM SOME SORT OF EXPRESSION MATRIX
# 	plot MDS of data - from limma/voom analysis pipeline
#	colour by subtype
#
# options:
# 		expr = some sort of expression matrix with rows=genes, cols=samples
# 		samp_id = vector of sample ids in same order as columns of expr
# 		samp_type = vector of sample subtypes in same order as columns of expr
# 		outname = name of output plot file
# 		title = title for output MDS plot
#		other_lab = vector of [some other info] of same length as #samples
#
####################################################################################################
my.labelled.mds.plot <- function(expr, samp_id, samp_type, outname, title, other_lab)
{

    # calculate and plot Euclidean sample distances
    # calculate Euclidean distances
    sampleDists <- dist( t( expr ) )

    # coerce distances to matrix
    sampleDistMatrix <- as.matrix( sampleDists )
    # define rownames
    rownames(sampleDistMatrix) <- paste( samp_id, samp_type, sep="-" )

    # generate mds data frame:
    mds <- data.frame(cmdscale(sampleDistMatrix))
    mds$class <- samp_type
    
    # plot as png
    # add label with row# 
	# use palette "Set1" from RColorBrewer
    png(outname, width=1024, height=1024)
    print( qplot(X1,X2, color=class, data=mds, main=title) + geom_text( aes(label=other_lab), size=8, hjust=1, vjust=1) + 		scale_colour_manual(values=subtype.colours) + theme_bw() ) 
    dev.off()
}
####################################################################################################

####################################################################################################
#
# convertIDs()
# 	function to interconvert gene IDs
#	need to load Bioc annotation libraries: AnnotationDbi, org.Hs.eg.db
#	from the RNAseq Bioconductor workflow tutorial
#
# options:
# 		ids = vector of gene IDs to convert
# 		from = keytype (see Bioc annotation docs) to convert from
# 		to = keytype to convert to
# 		db = annotation db to use (must contain keytypes in from/to)
# 		ifMultiple = how to deal with multi-mapping IDs: either NA, or use first match
#
####################################################################################################
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst"))
{
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}
####################################################################################################

####################################################################################################
#
# read.bed()
# 	generic function to read in headerless bed files
#
# args:
# 		bed.path = path to bed file to load
# 		enh.colnames = character vector of column headers for that bed file
#		header = T/F whether the bed file has a header line
#
####################################################################################################
read.bed <- function( bed.path, colnames, header=F )
{
    # read in the bed as a headerless table
    l <- read.table(file=bed.path, header=header, col.names=colnames)
    # then set up GR with chrom/start/end and metadata cols:
    # remove leading metadata columns (chrom/start/end) already in GR
    l <- GRanges( seqnames=Rle(l[,1]), ranges=IRanges(start=l[,2], end=l[,3]),
              mcols=data.frame( l[,-(1:3)] ) )
    # get rid of leading "mcols." in colnames
    colnames(mcols(l)) <- unlist( lapply( strsplit(colnames(mcols(l)), split="[.]"), "[", 2 ) )
    # return the GR obj
    return(l)
}
####################################################################################################