library(GetoptLong)
library(circlize)
library(ComplexHeatmap)
library(epik)
load_epik_config("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/script_publicate/hipo16_config_epik.R")

enforce_rerun = FALSE

methylation_hooks$set_chr("chr21")
sample_id = methylation_hooks$sample_id

######################## Figure 1C: the circular plot ######################
library(EnrichedHeatmap)
library(GenomicRanges)
window_size = 1e6
if(file.exists(qq("@{PROJECT_DIR}/rds/methylation_different_to_normal_w@{window_size}.RData")) && !enforce_rerun) {
	load(qq("@{PROJECT_DIR}/rds/methylation_different_to_normal_w@{window_size}.RData"))
} else {
	chr_len = read.chromInfo(species = "hg19")$chr.len[CHROMOSOME]
	gr_window = GRanges()
	for (chr in CHROMOSOME) {
	    methylation_hooks$set_chr(chr, verbose = FALSE)
	    meth = methylation_hooks$meth[, sample_id]
	    gr = methylation_hooks$gr
	    chr_gr = GRanges(seqnames = chr, ranges = IRanges(1, chr_len[chr]))
	    chr_window = makeWindows(chr_gr, w = window_size)
	    mtch = as.matrix(findOverlaps(chr_window, gr))
	    gr2 = chr_window[unique(mtch[, 1])]
	    meth = tapply(mtch[, 2], mtch[, 1], function(i) colMeans(meth[i, , drop = FALSE]))
	    meth = do.call("rbind", meth)
	    mcols(gr2) = meth
	    gr_window = c(gr_window, gr2)
	}
	names(gr_window) = NULL
	df_window = as.data.frame(gr_window)
	
	ind_normal = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "normal"])
	ind_IDH = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "IDH"])
	ind_MES = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "MES"])
	ind_RTK_I = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "RTK_I"])
	ind_RTK_II = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "RTK_II"])

	normal_mean = rowMeans(df_window[, ind_normal])
	m = as.matrix(df_window[, ind_IDH] - normal_mean); od = hclust(dist(t(m)))$order; df_window[, ind_IDH] = m[, od]
	m = as.matrix(df_window[, ind_MES] - normal_mean); od = hclust(dist(t(m)))$order; df_window[, ind_MES] = m[, od]
	m = as.matrix(df_window[, ind_RTK_I] - normal_mean); od = hclust(dist(t(m)))$order; df_window[, ind_RTK_I] = m[, od]
	m = as.matrix(df_window[, ind_RTK_II] - normal_mean); od = hclust(dist(t(m)))$order; df_window[, ind_RTK_II] = m[, od]

	save(df_window, file = qq("@{PROJECT_DIR}/rds/methylation_different_to_normal_w@{window_size}.RData"))
}

ind_normal = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "normal"])
ind_IDH = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "IDH"])
ind_MES = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "MES"])
ind_RTK_I = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "RTK_I"])
ind_RTK_II = which(colnames(df_window) %in% rownames(SAMPLE)[SAMPLE$subgroup == "RTK_II"])


pdf(qq("@{PROJECT_DIR}/image/circos_general_methylation.pdf"), width = 7, height = 7)
figure1_A = function() {
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
diff_col_fun = colorRamp2(c(-0.15, 0, 0.15), c("#3794bf", "#FFFFFF", "#df8640"))

circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(convert_height(0.5, "mm"), convert_height(0.5, "mm")),
	gap.after = c(rep(1, 21), 20), start.degree = 80)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), plotType = c("axis"), axis.labels.cex = 0.5)
circos.track(track.index = 1, panel.fun = function(x, y) {
	circos.text(CELL_META$xcenter, CELL_META$ylim[2]*2, gsub("chr", "", CELL_META$sector.index), 
		niceFacing = TRUE, facing = "inside", cex = 1)
}, bg.border = NA)
th = 0.008
#### IDH
circos.genomicTrack(df_window, track.height = th*length(ind_IDH), 
	stack = TRUE, numeric.column = ind_IDH, 
	panel.fun = function(region, value, ...) {
		i = getI(...)
		qqcat("IDH @{i}th sample, @{CELL_META$sector.index}\n")
		circos.genomicRect(region, value, col = diff_col_fun(value[[1]]), border = NA, ...)
})
draw.sector(start.degree = 99, end.degree = 81, rou1 = CELL_META$cell.top.radius, 
	rou2 = CELL_META$cell.bottom.radius, col = COLOR$subgroup["IDH"], border = NA)
text(0, mean(CELL_META$yplot), "IDH", cex = 0.8)

#### MES
circos.genomicTrack(df_window, track.height = th*length(ind_MES), 
	stack = TRUE, numeric.column = ind_MES, 
	panel.fun = function(region, value, ...) {
		i = getI(...)
		qqcat("MES @{i}th sample, @{CELL_META$sector.index}\n")
		circos.genomicRect(region, value, col = diff_col_fun(value[[1]]), border = NA, ...)
})
draw.sector(start.degree = 99, end.degree = 81, rou1 = CELL_META$cell.top.radius, 
	rou2 = CELL_META$cell.bottom.radius, col = COLOR$subgroup["MES"], border = NA)
text(0, mean(CELL_META$yplot), "MES",  cex = 0.8)


### RTK_I
circos.genomicTrack(df_window, track.height = th*length(ind_RTK_I), 
	stack = TRUE, numeric.column = ind_RTK_I, 
	panel.fun = function(region, value, ...) {
		i = getI(...)
		qqcat("RTK_I @{i}th sample, @{CELL_META$sector.index}\n")
		circos.genomicRect(region, value, col = diff_col_fun(value[[1]]), border = NA, ...)
})
draw.sector(start.degree = 99, end.degree = 81, rou1 = CELL_META$cell.top.radius, 
	rou2 = CELL_META$cell.bottom.radius, col = COLOR$subgroup["RTK_I"], border = NA)
text(0, mean(CELL_META$yplot), "RTK_I",  cex = 0.8)

#### RTK_II
circos.genomicTrack(df_window, track.height = th*length(ind_RTK_II), 
	stack = TRUE, numeric.column = ind_RTK_II, 
	panel.fun = function(region, value, ...) {
		i = getI(...)
		qqcat("RTK_II @{i}th sample, @{CELL_META$sector.index}\n")
		circos.genomicRect(region, value, col = diff_col_fun(value[[1]]), border = NA, ...)
})
draw.sector(start.degree = 99, end.degree = 81, rou1 = CELL_META$cell.top.radius, 
	rou2 = CELL_META$cell.bottom.radius, col = COLOR$subgroup["RTK_II"], border = NA)
text(0, mean(CELL_META$yplot), "RTK_II", cex = 0.8)

lgd = Legend(col_fun = diff_col_fun, at = c(-0.2, -0.1, 0, 0.1, 0.2), direction = "horizontal",
	title = "Methylation difference\nto normal", title_position = "topcenter")
grid.draw(lgd)
}
dev.off()

################## Figure 1D ##############

#### figure extra, mean methylation in genomic features ####
gf_list = get_mean_methylation_in_genomic_features(sample_id, GENOMIC_FEATURE_LIST)
gf_mean_meth = sapply(gf_list, function(gf) {
	m = mcols(gf)
	m = m[, grep("mean_meth", colnames(m))]
	m = as.matrix(m)
	colnames(m) = gsub("mean_meth_", "", colnames(m))
	tapply(seq_len(ncol(m)), SAMPLE[sample_id, "subgroup"], function(ind) mean(m[, ind], na.rm = TRUE))
})
gf_mean_meth = gf_mean_meth[, order(gf_mean_meth[1, ])]
pdf(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/image/mean_meth_in_gf.pdf"))
plot(NULL, xlim = c(0, ncol(gf_mean_meth)), ylim = c(0, 1), axes = FALSE,
	ylab = "mean methylation", xlab = "")
for(i in seq_len(ncol(gf_mean_meth))) {
	x = gf_mean_meth[, i]
	points(rep(i-0.5, 5) + runif(5, max = 0.05, min = -0.05), x, col = COLOR$subgroup[names(x)], pch = 16)
	text(i - 0.5, max(x) + 0.03, colnames(gf_mean_meth)[i], srt = 90, adj = c(0, 0.5))
}
axis(side = 2)
title("Mean methylation in genomic features")
box()
legend("bottomright", pch = 16, col = COLOR$subgroup[1:5], legend = names(COLOR$subgroup)[1:5])
dev.off()



####### Figure 1G #######
library(GetoptLong)

load_namespace("/desktop-home/guz/project/development/epik")
load_epik_config("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/script_publicate/hipo16_config_epik.R")
	


pushTrackList = function(trackList, track) {
    if(!is.null(track)) {
        trackList[[length(trackList) + 1]] = track
    }
    return(trackList)
}

constructAnnotationTrack = function(gr, chr, name = NULL, genome = "hg19", start = 0, end = Inf, fill = "blue", ...) {
    if(length(fill) == 1) fill = rep(fill, length(gr))
    l1 = as.vector(seqnames(gr) == chr)
    gr2 = gr[l1]
    fill = fill[l1]
    l2 = end(gr2) > start & start(gr2) < end
    gr2 = gr2[l2]
    fill = fill[l2]

    if(length(gr2)) {

        AnnotationTrack(name = name,
                        start = start(gr2),
                        end = end(gr2),
                        chromosome = seqnames(gr2),
                        genome = genome, 
                        stacking = "dense",
                        showTitle = TRUE, 
                        col = "transparent",
                        fill = fill,
                        height = unit(5, "mm"),
                        ...)
    } else {
        NULL
    }
}



chr = "chr20"
gene_start = 62769128
gene_end = 62880190

sample_id = rownames(SAMPLE)

	methylation_hooks$set_chr(chr, verbose = FALSE)

	gr = methylation_hooks$gr
	meth_mat = methylation_hooks$meth
	l = start(gr) >= gene_start & end(gr) <= gene_end
	gr = gr[l]
	meth_mat = meth_mat[l, , drop = FALSE]
	

chrom_list = get_chromHMM_list(c(rownames(SAMPLE)[1:60], paste0("E0", 67:74)))
chrom_list = lapply(chrom_list, function(x) x[seqnames(x) == chr])
## convert chrom_list to a matrix
window = makeWindows(GRanges(seqnames = chr, ranges = IRanges(62769100, 62880200)), w = 200)
chrom_mat = matrix(18, nrow = length(window), ncol = length(chrom_list))
colnames(chrom_mat) = names(chrom_list)
for(i in seq_along(chrom_list)) {
	qqcat("@{i}/@{length(chrom_list)}...\n")
	mtch = as.matrix(findOverlaps(window, chrom_list[[i]]))
	chrom_mat[mtch[, 1], i] = as.integer(gsub("E", "", chrom_list[[i]][ mtch[, 2] ]$states))
}


	message("add transcripts to gviz tracks...")
    options(Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/")
    trackList = list()
    grtrack = GeneRegionTrack(TXDB, chromosome = chr, start = gene_start, end = gene_end, 
        name="Gene\nmodel", showId = TRUE, rotate.title = TRUE, col = NA, showTitle = FALSE,
        size = 0.5)
    
    trackList = pushTrackList(trackList, grtrack)

    message("add methylation to gviz tracks...")
    subgroup = SAMPLE$subgroup

    for(t in unique(subgroup)) {
        mat = meth_mat[, subgroup == t]
        co = as.hclust(reorder(as.dendrogram(hclust(dist(t(mat[apply(mat, 1, function(x) all(!is.na(x))), ])))), wts = colMeans(mat, na.rm = TRUE)))$order
        trackList = pushTrackList(trackList, DataTrack(name = t,
                                    start = start(gr),
                                    end = end(gr),
                                    chromosome = seqnames(gr),
                                    genome = "hg19",
                                    data = t(mat[, co]),
                                    type = "heatmap",
                                    ylim = c(0, 1),
                                    showSampleNames = FALSE,
                                    gradient = c("blue", "white", "red"),
                                    size = 0.15*ncol(mat),
                                    col = NA,
                                    panelFun = local({t=t;function() {grid.text(t, x = unit(-2, "mm"), just = "right")}})
                                ))
    }

    for(t in unique(subgroup)) {
    	if(t == "normal") {
    		mat = chrom_mat[, intersect(colnames(chrom_mat), paste0("E0", 67:74))]
    	} else {
       		mat = chrom_mat[, intersect(colnames(chrom_mat), sample_id[subgroup == t])]
        }
        trackList = pushTrackList(trackList, DataTrack(name = t,
                                    start = start(window),
                                    end = end(window),
                                    chromosome = seqnames(window),
                                    genome = "hg19",
                                    data = t(mat),
                                    type = "heatmap",
                                    ylim = c(1, 18),
                                    showSampleNames = FALSE,
                                    gradient = state_col,
                                    size = 0.15*ncol(mat),
                                    col = NA,
                                    panelFun = local({t=t;function() {grid.text(t, x = unit(-2, "mm"), just = "right")}})
                                ))
    }

 

pdf(qq("@{PROJECT_DIR}/image/MYT1_chr20_62769128-62880190.pdf"), width = 12, height = 6)
# convert fontsize to cex
    grid.newpage()
    pushViewport(viewport(width = 0.8, height =0.8))
    plotTracks(trackList, from = gene_start, to = gene_end, chromosome = chr, panel.only = TRUE,
        main = title, cex.main = 1)
    popViewport()
    pushViewport(viewport(y = 0.9, just = "bottom", height = 0.1))
    popViewport()

dev.off()

    
