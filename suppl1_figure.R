
suppressPackageStartupMessages(library(GetoptLong))

library(epik)
load_epik_config("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/script_publicate/hipo16_config_epik.R")

library(EnrichedHeatmap)
library(matrixStats)
library(ComplexHeatmap)
library(cola)
library(circlize)

data_type = data.frame(WGBS = rep(1, nrow(SAMPLE)),
	RNASeq = rep(1, nrow(SAMPLE)),
	"H3K27ac" = 0,
	"H3K27me3" = 0,
	"H3K36me3" = 0,
	"H3K4me1" = 0,
	"H3K4me3" = 0,
	"H3K9me3" = 0)
rownames(data_type) = rownames(SAMPLE)
for(mk in MARKS) {
	sample_id = chipseq_hooks$sample_id(mk)
	data_type[sample_id, mk] = 1
}

SAMPLE$subgroup = factor(SAMPLE$subgroup, levels = c("IDH", "MES", "RTK_I", "RTK_II", "normal"), ordered = TRUE)
pdf(qq("@{PROJECT_DIR}/image/data_matrix.pdf"), width = 5, height = 8)
rect_gp = gpar(col = "white", type = "none")
cell_fun = function(j, i, x, y, w, h, fill) {
	# grid.roundrect(x, y, w, h, r = unit(0.8, "mm"), gp = gpar(fill = fill, col = "white"))
	grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "white"))
}
ht = Heatmap(SAMPLE[, 1], name = "Subtype", show_row_names = FALSE, col = COLOR$subgroup,
	combined_name_fun = NULL) + 
Heatmap(data_type[, 1, drop = F], col = c("0" = "grey", "1" = "#7570B3"), 
	show_row_names = FALSE, rect_gp = rect_gp, cell_fun = cell_fun,
	show_heatmap_legend = FALSE) +
Heatmap(data_type[, 2, drop = F], col = c("0" = "grey", "1" = "#7570B3"), 
	show_row_names = FALSE, rect_gp = rect_gp, cell_fun = cell_fun,
	show_heatmap_legend = FALSE) +
Heatmap(data_type[, -(1:2)], show_row_names = FALSE, cluster_columns = FALSE,
	col = c("0" = "grey", "1" = "#7570B3"), rect_gp = rect_gp, cell_fun = cell_fun,
	show_heatmap_legend = FALSE)
figureS1_A = grid.grabExpr(draw(ht, split = SAMPLE[, 1], row_order = order(SAMPLE[, 1], data_type[, 3], data_type[, 5]),
	heatmap_legend_list = list(Legend(at = 1:2, labels = c("available", "not available"), 
		legend_gp = gpar(fill = c("#7570B3", "grey"), col = "white"), title = "Status")),
	column_title = "Data matrix")
)
dev.off()


subgroup = c("IDH", "MES", "RTK_I", "RTK_II")
SAMPLE = SAMPLE[SAMPLE$subgroup %in% subgroup, , drop = FALSE]

sample_id = rownames(SAMPLE)

#### generate a list of genomic features ###
gene = GENOMIC_FEATURE_LIST[["gene"]]
intergenic = GENOMIC_FEATURE_LIST[["intergenic"]]
strand(gene) = "*"
mcols(gene) = NULL
gr = reduce(c(gene, intergenic))

GENOMIC_FEATURE_LIST$non_cgi = setdiff(gr, c(GENOMIC_FEATURE_LIST$cgi))
GENOMIC_FEATURE_LIST$non_cgi_and_shore = setdiff(gr, c(GENOMIC_FEATURE_LIST$cgi, GENOMIC_FEATURE_LIST$cgi_shore))

# GENOMIC_FEATURE_LIST = GENOMIC_FEATURE_LIST[c("cgi", "cgi_shore", "non_cgi", "non_cgi_and_shore")]
GENOMIC_FEATURE_LIST = GENOMIC_FEATURE_LIST[c("cgi", "non_cgi")]

###### split each genomic feature by 1kb window ##########
windowsize = 1000
library(EnrichedHeatmap)
gr_list = lapply(GENOMIC_FEATURE_LIST, function(gr) {
	gr = makeWindows(gr, w = windowsize, short.keep = TRUE)
	gr[width(gr) >= windowsize/2]
})

##### we don't need that many regions
gr_list$non_cgi = gr_list$non_cgi[sample(c(TRUE, FALSE), length(gr_list$non_cgi), p = c(0.1, 0.9), replace = TRUE)]
# gr_list$non_cgi_and_shore = gr_list$non_cgi_and_shore[sample(c(TRUE, FALSE), length(gr_list$non_cgi_and_shore), p = c(0.1, 0.9), replace = TRUE)]

##### generate the methylation matrix #########
gr_list = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, genomic_features = gr_list)
mat_list = lapply(gr_list, function(gr) {
	m = as.matrix(mcols(gr))
	m = m[apply(m, 1, function(x) all(!is.na(x))), grep("mean_meth", colnames(m))]
	m
})

#### perform consensus clustering ###########
library(cola)

set.seed(123)
res_list = lapply(mat_list, function(mat) {
	hierarchical_partition(mat, top_value_method = "sd", partition_method = "kmeans", scale_rows = FALSE,
		anno = data.frame(subtype = SAMPLE$subgroup), anno_col = list(subtype = COLOR$subgroup))
})

set.seed(123)
mat_expr = adjust_matrix(EXPR[, sample_id])
res_list$expr = hierarchical_partition(mat_expr, 
	anno = data.frame(subtype = SAMPLE$subgroup), anno_col = list(subtype = COLOR$subgroup))

df = as.data.frame(do.call("cbind", lapply(res_list, get_classes)))
df = cbind(df, subtype = SAMPLE$subgroup)

saveRDS(res_list, qq("@{PROJECT_DIR}/rds/unsupervised_clustering_methylation_cgi_and_non_cgi.rds"))
library(ComplexHeatmap)
library(circlize)


########### predict subtypes based on publicated signature genes ###############
library(ComplexHeatmap)
library(cola)

gene_name = extract_field_from_gencode(GTF_FILE)

### read in the signature gens as well as the expression #####
set.seed(12345)
gbm_2017 = read.table("/home/guz/project/analysis/hipo16_new/tcga_signature/gbm_2017.txt", header = TRUE, row.names = 1)
colnames(gbm_2017) = gsub("avg_expr_", "", colnames(gbm_2017))
tcga_2010 = read.table("/home/guz/project/analysis/hipo16_new/tcga_signature/tcga_microarray.txt", header = TRUE, row.names = 1, sep = "\t")[, -1]
gbm_2017 = as.matrix(gbm_2017)
tcga_2010 = as.matrix(tcga_2010)

# note there can be duplicated gene symbol!
mat = EXPR[, sample_id]
rownames(mat) = gene_name[rownames(mat)]
mat = cola:::adjust_matrix(mat)
mat = t(scale(t(mat)))

gbm_2017 = as.data.frame(gbm_2017)
gbm_2017$class = colnames(gbm_2017)[1:3][epik:::rowWhichMax(as.matrix(gbm_2017[, 1:3]))]

cn = intersect(rownames(mat), rownames(gbm_2017))
dist_to_signatures = as.matrix(pdist::pdist(t(mat[cn, ]), t(gbm_2017[cn, 1:3])))
colnames(dist_to_signatures) = colnames(gbm_2017)[1:3]
diff_ratio = apply(dist_to_signatures, 1, function(x) { 
	y = diff(sort(x))
	y[1]/y[2]
})
diff_ratio[diff_ratio > 3] = 3
predicted_class = colnames(gbm_2017)[epik:::rowWhichMax(-dist_to_signatures)]
predicted_class[diff_ratio < 1] = NA
pdf(qq("@{PROJECT_DIR}/image/subgroup_by_gbm_2017.pdf"), width = 10, height = 10)
ht_list = Heatmap(mat[cn, ], top_annotation = HeatmapAnnotation(dist = dist_to_signatures,
	diff_ratio = anno_barplot(diff_ratio, axis = TRUE),
	predicted_class = predicted_class,
	subgroup = SAMPLE$subgroup,
	col = COLOR,
	show_annotation_name = TRUE, annotation_height = unit(c(2, 3, 0.5, 0.5), "cm")),
	show_column_names = FALSE, show_row_names = FALSE
) +
Heatmap(gbm_2017[cn, 1:3], cluster_columns = FALSE, show_row_names = FALSE,
	width = unit(4, "cm"))
draw(ht_list, split = gbm_2017[cn, ]$class, column_title = "GBM 2017 classification")
decorate_annotation("diff_ratio", {
	grid.lines(unit(c(0, 1), "npc"), unit(c(1,1), "native"), gp = gpar(lty = 2))
})
dev.off()

predicted_by_gbm_2017 = predicted_class

## only use signiture genes which are highly expressed
set.seed(123)
km = kmeans(tcga_2010, centers = 2)$cluster
tm = tapply(seq_len(nrow(tcga_2010)), km, function(x) mean(tcga_2010[x, ]))
tm = tm[order(names(tm))]
tcga_2010 = tcga_2010[km == which.max(tm), ]

tcga_2010 = as.data.frame(tcga_2010)
tcga_2010$class = colnames(tcga_2010)[1:4][epik:::rowWhichMax(as.matrix(tcga_2010[, 1:4]))]

cn = intersect(rownames(mat), rownames(tcga_2010))
dist_to_signatures = as.matrix(pdist::pdist(t(mat[cn, ]), t(tcga_2010[cn, 1:4])))
colnames(dist_to_signatures) = colnames(tcga_2010)[1:4]
diff_ratio = apply(dist_to_signatures, 1, function(x) { 
	y = diff(sort(x))
	y[1]/y[2]
})
diff_ratio[diff_ratio > 3] = 3
predicted_class = colnames(tcga_2010)[epik:::rowWhichMax(-dist_to_signatures)]
predicted_class[diff_ratio < 1] = NA
pdf(qq("@{PROJECT_DIR}/image/subgroup_by_tcga_2010.pdf"), width = 10, height = 10)
ht_list = Heatmap(mat[cn, ], top_annotation = HeatmapAnnotation(dist = dist_to_signatures,
	diff_ratio = anno_barplot(diff_ratio, axis = TRUE),
	predicted_class = predicted_class,
	subgroup = SAMPLE$subgroup,
	col = COLOR,
	show_annotation_name = TRUE, annotation_height = unit(c(2, 3, 0.5, 0.5), "cm")),
	show_column_names = FALSE, show_row_names = FALSE
) +
Heatmap(tcga_2010[cn, 1:4], cluster_columns = FALSE, show_row_names = FALSE, width = unit(4, "cm"))
draw(ht_list, split = tcga_2010[cn, ]$class, column_title = "tcga_2010 classification")
decorate_annotation("diff_ratio", {
	grid.lines(unit(c(0, 1), "npc"), unit(c(1,1), "native"), gp = gpar(lty = 2, col = "grey"))
})
dev.off()
predicted_by_tcga_2010 = predicted_class

### Ceccarelli	M,	Cell	2016;	164:550-63;	PMID:	26824661	
# find differential genes between subgroups as signature genes
mm = read.table("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/illuminahiseq_rnaseqv2-RSEM_genes_normalized.tab", header = TRUE, sep = "\t", check.names = FALSE)
gi = gsub("^.*\\|", "", mm[, 1])
library(org.Hs.eg.db)
x <- org.Hs.egENSEMBL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
xx = xx[sapply(xx, length) == 1]
xx = unlist(xx)
l = gi %in% names(xx)
mm = mm[l, -1]
gi = gi[l]
l = !duplicated(xx[gi])
mm = mm[l, ]
rownames(mm) = xx[gi[l]]
colnames(mm) = gsub("^(\\w+-\\w+-\\w+)?-.*$", "\\1", colnames(mm))
anno = read.table("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/tcga_lgg_subtypes.tab", row.names = 1, header = TRUE, sep = "\t")
cn = intersect(colnames(mm), rownames(anno))
mm = mm[, cn]
anno = anno[cn, "Pan.Glioma.RNA.Expression.Cluster"]
l = !is.na(anno)
mm = mm[, l]
anno = anno[l]
l = anno != "unclassified"
mm = mm[, l]
anno = anno[l]
mm = as.matrix(mm)
sig_df = do.call("rbind", apply(mm, 1, function(x) {
	gm = tapply(x, anno, mean)
	gn = names(gm)
	gmax = gn[which.max(gm)]
	p = NULL
	for(i in setdiff(names(gm), gmax)) {
		p[[i]] = t.test(x[anno == gmax], x[anno == i])$p.value
	}
	data.frame(group = gmax, max_p = max(p, na.rm = TRUE), stringsAsFactors = FALSE)
}))

sig = do.call(rbind, tapply(seq_len(nrow(sig_df)), sig_df$group, function(ind) {
	ind = ind[order(sig_df$max_p[ind])[1:100]]
	data.frame(group = sig_df$group[ind], index = ind)
}))

Heatmap(t(scale(t(mm[sig$ind, ]))), show_row_names = FALSE, show_column_names = FALSE,
	column_split = anno, row_split = sig$group,
	top_annotation = HeatmapAnnotation(subtype = anno)) + rowAnnotation(group = sig$group)

m_scale = t(scale(t(mm[sig$ind, ])))
sig_expr = do.call("cbind", tapply(seq_len(ncol(m_scale)), anno, function(ind) rowMeans(m_scale[, ind])))

mat = EXPR[, sample_id]
rownames(mat) = gsub("\\.\\d+$", "", rownames(mat))
mat = cola:::adjust_matrix(mat)
mat = t(scale(t(mat)))

tcga_2016 = as.data.frame(sig_expr)
tcga_2016$class = colnames(tcga_2016)[1:4][epik:::rowWhichMax(as.matrix(tcga_2016[, 1:4]))]

cn = intersect(rownames(mat), rownames(tcga_2016))
dist_to_signatures = as.matrix(pdist::pdist(t(mat[cn, ]), t(tcga_2016[cn, 1:4])))
colnames(dist_to_signatures) = colnames(tcga_2016)[1:4]
diff_ratio = apply(dist_to_signatures, 1, function(x) { 
	y = diff(sort(x))
	y[1]/y[2]
})
diff_ratio[diff_ratio > 3] = 3
predicted_class = colnames(tcga_2016)[epik:::rowWhichMax(-dist_to_signatures)]
predicted_class[diff_ratio < 1] = NA
pdf(qq("@{PROJECT_DIR}/image/subgroup_by_tcga_2016.pdf"), width = 10, height = 10)
ht_list = Heatmap(mat[cn, ], top_annotation = HeatmapAnnotation(dist = dist_to_signatures,
	diff_ratio = anno_barplot(diff_ratio, axis = TRUE),
	predicted_class = predicted_class,
	subgroup = SAMPLE$subgroup,
	col = COLOR,
	show_annotation_name = TRUE, annotation_height = unit(c(2, 3, 0.5, 0.5), "cm")),
	show_column_names = FALSE, show_row_names = FALSE
) +
Heatmap(tcga_2016[cn, 1:4], cluster_columns = FALSE, show_row_names = FALSE,
	width = unit(4, "cm"))
draw(ht_list, split = tcga_2016[cn, ]$class, column_title = "TCGA 2016 classification")
decorate_annotation("diff_ratio", {
	grid.lines(unit(c(0, 1), "npc"), unit(c(1,1), "native"), gp = gpar(lty = 2))
})
dev.off()

predicted_by_tcga_2016 = predicted_class





res_list = readRDS(qq("@{PROJECT_DIR}/rds/unsupervised_clustering_methylation_cgi_and_non_cgi.rds"))
df = as.data.frame(do.call("cbind", lapply(res_list, get_classes)))
df$gbm_2017 = predicted_by_gbm_2017
df$tcga_2010 = predicted_by_tcga_2010
df = cbind(df, subtype = SAMPLE$subgroup)

compare_partition = function(x1, x2, method = "jaccard") {
	get_stat = function(x1, x2, method) {
		cl1 = as.cl_hard_partition(x1)
		cl2 = as.cl_hard_partition(x2)
		cl = cl_ensemble(cl1, cl2)
		cl_agreement(cl, method = method)[[1]]
	}

	s = get_stat(x1, x2, method)
	sr = numeric(1000)
	for(i in 1:1000) {
		sr[i] = get_stat(x1, sample(x2, length(x2)), method)
	}

	sum(sr >= s)/1000
}

### add hierarchy as legends
library(data.tree)
library(grid)
library(dendextend)
subgroup_dend = function(object, hierarchy = object@hierarchy) {

	lt = list()
	lt[["0"]] = Node$new("all samples")
	cn = colnames(object@list[["0"]]@.env$data)
	max_depth = max(nchar(hierarchy))
	lt[["0"]]$node_height = max_depth + 1
	for(i in seq_len(nrow(hierarchy))) {
		lt[[ hierarchy[i, 2] ]] = lt[[ hierarchy[i, 1] ]]$AddChildNode({
			node = Node$new(hierarchy[i, 2])
			node$node_height = max_depth - nchar(hierarchy[i, 2]) + 1
			node
		})
		l = hierarchy[, 1] == hierarchy[i, 2]
	}
	dend = as.dendrogram(lt[["0"]], heightAttribute = "node_height")

	od = structure(1:nleaves(dend), names = labels(dend))
	dend_env = new.env()
	dend_env$dend = dend

	update_midpoint = function(index = NULL) {
		if(is.null(index)) {
			if(is.leaf(dend_env$dend)) {
				pos = od[attr(dend_env$dend, "label")]
				midpoint = 0
			} else {
				x = NULL
				for(i in seq_len(length(dend_env$dend))) {
					if(is.null(attr(dend_env$dend[[i]], "x"))) {
						update_midpoint(i)
					}
					x[i] = attr(dend_env$dend[[i]], "x")
				}
				pos = (max(x) + min(x))/2
				midpoint = (max(x) - min(x))/2
			}
		} else {
			if(is.leaf(dend_env$dend[[index]])) {
				pos = od[attr(dend_env$dend[[index]], "label")]
				midpoint = 0
			} else {
				x = NULL
				for(i in seq_len(length(dend_env$dend[[index]]))) {
					if(is.null(attr(dend_env$dend[[c(index, i)]], "x"))) {
						update_midpoint(c(index, i))
					}
					x[i] = attr(dend_env$dend[[c(index, i)]], "x")
				}
				pos = (max(x) + min(x))/2
				midpoint = (max(x) - min(x))/2
			}
		}
		if(is.null(index)) {
			attr(dend_env$dend, "x") = pos
		} else {
			attr(dend_env$dend[[index]], "x") = pos
			attr(dend_env$dend[[index]], "midpoint") = midpoint
		}
	}
	update_midpoint()

	# reorder(dend_env$dend, wts = order(labels(dend_env$dend)))
	dend_env$dend
}

legend_with_dend = function(dend, title, col = NULL, labels = get_leaves_attr(dend, "label")) {
	n = length(labels)
	if(is.null(col)) {
		col = structure(rand_color(n), names = labels)
	}
	col = col[labels]

	obj = grid.grabExpr({
		
		pushViewport(viewport(width = unit(1.05, "cm") + unit(4, "mm") + max_text_width(labels, gp = gpar(fontsize = 10)),
			height = max_text_height(title) + unit(1.2, "mm") + unit(4, "mm")*n))
		pushViewport(viewport(y = 1, height = max_text_height(title) + unit(1.2, "mm"), just = "top"))
		grid.text(title, x = 0, just = "left", gp = gpar(fontsize = 10, fontface = "bold"))
		popViewport()

		pushViewport(viewport(x = 0, y = 0, width = unit(4, "mm") + max_text_width(labels), height = unit(4, "mm")*n, just = c("left", "bottom")))
		grid.rect(x = 0, y = (1:n-0.5)/n, width = unit(4, "mm"), height = unit(4, "mm"), just = "left", gp = gpar(fill = col, col = "white"))
		grid.text(labels, x = unit(5, "mm"), y = (1:n-0.5)/n, just = "left", gp = gpar(fontsize = 10))
		popViewport()
		
		pushViewport(viewport(x = 1, y = 0, width = unit(1, "cm"), height = unit(4, "mm")*n, just = c("right", "bottom")))
		grid.dendrogram(dend, facing = "left")
		popViewport()

		popViewport()

	})

	gf = frameGrob(layout = grid.layout(nrow = 1, ncol = 1,
            widths = unit(1.05, "cm") + unit(4, "mm") + max_text_width(labels, gp = gpar(fontsize = 10)), 
            heights = max_text_height(title) + unit(1.2, "mm") + unit(4, "mm")*n))
	gf = placeGrob(gf, obj, 1, 1)
	gf
}

df2 = df[order(df$subtype, df$expr, df$gbm_2017, df$tcga_2010, df$cgi), ]
ht_list = Heatmap(df2[, "subtype", drop = FALSE], name = "subtype",
	col = COLOR$subgroup, show_row_names = FALSE, show_heatmap_legend = FALSE)
ht_list = ht_list + Heatmap(df2[, "cgi"], name = "cgi", show_row_names = FALSE,
	col = c("01" = COLOR$subgroup[["IDH"]],
		    "001" = COLOR$subgroup[["MES"]],
		    "000" = COLOR$subgroup[["RTK_II"]]),
	show_heatmap_legend = FALSE)
ht_list = ht_list + Heatmap(df2[, "non_cgi"], name = "non_cgi", show_row_names = FALSE,
	col = c("002" = COLOR$subgroup[["IDH"]],
		    "0001" = COLOR$subgroup[["MES"]],
		    "02" = COLOR$subgroup[["RTK_I"]],
		    "0000" = COLOR$subgroup[["RTK_II"]]),
	show_heatmap_legend = FALSE)
ht_list = ht_list + Heatmap(df2[, "expr"], name = "expr", show_row_names = FALSE,
	col = c("000" = COLOR$subgroup[["IDH"]],
		    "031" = COLOR$subgroup[["MES"]],
		    "001" = COLOR$subgroup[["RTK_II"]],
		    "030" = "#A65628"),
	show_heatmap_legend = FALSE)
ht_list = ht_list + Heatmap(df2[, "tcga_2010"], name = "tcga_2010", show_row_names = FALSE,
	col = c("Proneural" = COLOR$subgroup[["IDH"]],
		    "Classical" = COLOR$subgroup[["RTK_II"]],
		    "Mesenchymal" = COLOR$subgroup[["MES"]],
		    "Neural" = "#A65628"),
	show_heatmap_legend = FALSE)
ht_list = ht_list + Heatmap(df2[, "gbm_2017"], name = "gbm_2017", show_row_names = FALSE,
	col = c("Proneural" = COLOR$subgroup[["IDH"]],
		    "Classical" = COLOR$subgroup[["RTK_II"]],
		    "Mesenchymal" = COLOR$subgroup[["MES"]]),
	show_heatmap_legend = FALSE)
lgd0 = color_mapping_legend(ht_list@ht_list$subtype@matrix_color_mapping, plot = FALSE)
dend = subgroup_dend(res_list$cgi)
attr(dend[[1]], "height") = 1
lgd1 = legend_with_dend(dend, title = "CGI",
	col = c("01" = COLOR$subgroup[["IDH"]],
		    "001" = COLOR$subgroup[["MES"]],
		    "000" = COLOR$subgroup[["RTK_II"]]))
dend = subgroup_dend(res_list$non_cgi)
attr(dend[[2]], "height") = 3
lgd2 = legend_with_dend(dend, title = "non_cgi",
	col = c("002" = COLOR$subgroup[["IDH"]],
		    "0001" = COLOR$subgroup[["MES"]],
		    "02" = COLOR$subgroup[["RTK_I"]],
		    "0000" = COLOR$subgroup[["RTK_II"]]))
dend = subgroup_dend(res_list$expr)
attr(dend[[2]], "height") = 3
lgd3 = legend_with_dend(dend, title = "expr",
	col = c("000" = COLOR$subgroup[["IDH"]],
		    "001" = COLOR$subgroup[["RTK_II"]],
		    "031" = COLOR$subgroup[["MES"]],
		    "030" = "#A65628"))
lgd4 = color_mapping_legend(ht_list@ht_list$tcga_2010@matrix_color_mapping, plot = FALSE)
lgd5 = color_mapping_legend(ht_list@ht_list$gbm_2017@matrix_color_mapping, plot = FALSE)

pdf(qq("@{PROJECT_DIR}/image/compare_subtypes.pdf"), width = 4, height = 8)
figureS1_C = grid.grabExpr(
draw(ht_list, heatmap_legend_list = list(lgd0, lgd1, lgd2, lgd3, lgd4, lgd5),
	column_title = "compare subtypes")
)
dev.off()


### circular comparison ####
figureS1_C = function() {
tb = table(as.character(df$subtype))
par(mfrow = c(3, 2))
col_list = list(
	structure(brewer.pal(4, "Pastel1")[-3], names = c("01", "001", "000")),
	structure(brewer.pal(4, "Pastel1"), names = c("002", "0001", "02", "0000")),
	structure(brewer.pal(5, "Pastel1")[-3], names = c("000", "031", "001", "030")),
	structure(c(brewer.pal(4, "Pastel1")[-3], "white"), names = c("Proneural", "Mesenchymal", "Classical", "NA")),
	structure(c(brewer.pal(5, "Pastel1")[-3], "white"), names = c("Proneural", "Mesenchymal", "Classical", "Neural", "NA"))
	structure(c(brewer.pal(5, "Pastel1")[-3], "white"), names = c("LGr2", "LGr4", "LGr1", "LGr3", "NA"))
)
for(i in 1:5) {
	circos.par(cell.padding = c(0, 0), track.margin = c(0, 0), gap.degree = 10)
	circos.initialize(names(tb), xlim = cbind(rep(0, length(tb)), tb))
	circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
		x = df[df$subtype == CELL_META$sector.index, i]
		tb2 = sort(table(x, useNA = "ifany"), decreasing = TRUE)
		tb2 = c(tb2[!is.na(names(tb2))], tb2[is.na(names(tb2))])
		names(tb2)[is.na(names(tb2))] = "NA"
		n = length(tb2)
		cm = cumsum(tb2)
		circos.rect(c(0, cm[-n]), rep(0, n), cm, rep(1, n), col = col_list[[i]][names(tb2)], border = "white")
	}, bg.border = NA)
	circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
		n = CELL_META$xlim[2]
		circos.rect(0:(n-1), rep(0, n), 1:n, rep(1, n), col = COLOR$subgroup[CELL_META$sector.index], border = "white")
		circos.text(CELL_META$xcenter, -0.5, CELL_META$sector.index, col = "black", niceFacing = TRUE, facing = "inside")
	}, bg.border = NA)
	text(0, 0, colnames(df)[i])
	circos.clear()
}
}

################### subgroup by chromatin states ##################

state_name = c("TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "Tx", "TxWk", "EnhG1", "EnhG2", "EnhA1", "EnhA2", "EnhWk", "ZNF/Rpts", "Het",
    "TssBiv", "EnhBiv", "ReprPC", "ReprPCWk", "Quies")
names(state_name) = paste0("E", ifelse(nchar(1:18) == 1, "0", ""), 1:18)

state_col = c("#FF0000", "#FF4500", "#FF4500", "#FF4500", "#008000", "#006400",
    "#C2E105", "#C2E105", "#FFC34D", "#FFC34D", "#FFFF00", "#66CDAA", "#8A91D0",
    "#CD5C5C", "#BDB76B", "#808080", "#C0C0C0", "#000000")
names(state_col) = names(state_name)

chrom_list = get_chromHMM_list(rownames(SAMPLE[SAMPLE$subgroup != "normal", , drop =FALSE]))
chrom_list = get_chromHMM_list(rownames(SAMPLE))

pdf(qq("@{PROJECT_DIR}/image/chrom_mds.pdf"))
sample_id = names(chrom_list)
for(st in names(state_name)) {	
	cat(st, "\n")
	n = length(chrom_list)
	mat = matrix(1, nrow = n, ncol = n)
	for(i in 1:(n-1)) {
	    for(j in (i+1):n) {
	        mat[i, j] = genomic_corr_jaccard(chrom_list[[i]][chrom_list[[i]]$states == st], chrom_list[[j]][chrom_list[[j]]$states == st])
	        mat[j, i] = mat[i, j]
	    }
	}
	loc = cmdscale(as.dist(1 - mat), k = 2)
	plot(loc, pch = 16, col = COLOR$subgroup[SAMPLE[sample_id, ]], main = qq("@{st}, @{state_name[st]}, dissimilarity = 1 - jaccard"))
}
dev.off()


	st = "E01"
	n = length(chrom_list)
	mat = matrix(1, nrow = n, ncol = n)
	for(i in 1:(n-1)) {
	    for(j in (i+1):n) {
	        mat[i, j] = genomic_corr_jaccard(chrom_list[[i]][chrom_list[[i]]$states == st], chrom_list[[j]][chrom_list[[j]]$states == st])
	        mat[j, i] = mat[i, j]
	    }
	}
	loc = cmdscale(as.dist(1 - mat), k = 2)
figureS1_D = local({
	loc = loc
	sample_id = names(chrom_list)
	function() {
		par(mar = c(4.1, 4.1, 4.1, 1))
		plot(loc, pch = 16, cex = 1.5, col = COLOR$subgroup[SAMPLE[sample_id, ]], main = qq("MDS of E01, active TSS\ndissimilarity = 1 - jaccard"),
		xlab = "Dimension 1", ylab = "Dimension 2", axes = FALSE)
		axis(side = 1, at = c(-0.2, -0.1, 0, 0.1, 0.2))
		axis(side = 2, at = c(-0.2, -0.1, 0, 0.1))
		box()
	}
})

library(gridGraphics)

obj_names = load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/methylation_array_subtyping/2018-02-26_methylation_clustering_end_Renv.Rdata")

title <- paste("Hierarchical clustering of n=", dim(dat.filt)[2], " GBM samples (methylation array data)\n", 
	dim(dat.filt)[1], " probes from Sturm et al. (2012) common to 450k and EPIC platforms; Euclidean distance, hclust method: ", hclust.method, sep="")
# for muts map as described
ha.bottom <- ComplexHeatmap::HeatmapAnnotation( 
	Age = anno_points(cnas.muts.gx$age, pch = 16, ylim = c(0, max(cnas.muts.gx[, 1])), axis = TRUE, axis_side = "left"),
    Gender = cnas.muts.gx$gender,
    "IDH1/2" = cnas.muts.gx$mutIDH1_2,
    TP53 = cnas.muts.gx$mutTP53,
    "Chr7 gain" = cnas.muts.gx$chr7gain,
    "Chr10 loss" = cnas.muts.gx$chr10loss,
    "Chr19 gain" = cnas.muts.gx$chr19gain,
    "Chr20 gain" = cnas.muts.gx$chr20gain,
    "EGFR amp" = cnas.muts.gx$ampEGFR,
    "PDGFRA amp" = cnas.muts.gx$ampPDGFRA,
    "CDKN2A/B del" = cnas.muts.gx$delCDKN2A_B,
    "CDK4 amp" = cnas.muts.gx$ampCDK4,
    "RB1 del" = cnas.muts.gx$delRB1,
    "MET amp" = cnas.muts.gx$ampMET,
    show_legend = FALSE,
    col=list(
		Gender=c("m"="#1B9E77", "f"="#E2798A"),
		"Chr7 gain"=c("normal"="white", "gain"="#E41A1C"),
		"Chr10 loss"=c("normal"="white", "loss"="#377EB8"),
		"Chr19 gain"=c("normal"="white", "gain"="#E41A1C"),	
		"Chr20 gain"=c("normal"="white", "gain"="#E41A1C"),	
		"IDH1/2"=c("wt"="white", "mut"="black", "R132H"="#377EB8"),
		TP53=c("wt"="white", "mut"="black"),
		"EGFR amp"=c("normal"="white", "amp"="#E41A1C"),		
		"PDGFRA amp"=c("normal"="white", "amp"="#E41A1C"),		
		"CDKN2A/B del"=c("normal"="white", "del"="#377EB8"),		
		"CDK4 amp"=c("normal"="white", "amp"="#E41A1C"),		
		"RB1 del"=c("normal"="white", "del"="#377EB8"),	
		"MET amp"=c("normal"="white", "amp"="#E41A1C") 
	),
    na_col="grey", show_annotation_name= c(FALSE, rep(TRUE, 13)), annotation_name_side = "left",
    annotation_height = unit(c(2, rep(0.5, 13)), "cm")
 )
ha <- HeatmapAnnotation( 
    df=data.frame(Subtype=predicted.subtypes),
    col=list(Subtype=COLOR$subgroup),
    na_col="black", show_annotation_name=TRUE, show_legend=FALSE,
    annotation_name_side = "left"
 )
 
hm <- Heatmap( matrix=as.matrix(dat.filt), col= colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
  column_title=title, name="Methylation", row_title=paste( dim(dat.filt)[1], "probes" ),
  cluster_rows=row.clust, cluster_columns=dat.clust, show_row_names=FALSE, show_column_names=FALSE,
  show_row_dend = FALSE,
  use_raster = TRUE, raster_quality = 2,
  top_annotation=ha, bottom_annotation=ha.bottom )
lgd_list = list(
	Legend(at = 1:4, labels = names(COLOR$subgroup)[-5], title = "Subtype", legend_gp = gpar(fill = COLOR$subgroup[-5])),
	Legend(at = 1:2, labels = c("Male", "Female"), title = "Gender", legend_gp = gpar(fill = c("#1B9E77", "#E2798A"))),
	Legend(at = 1:3, labels = c("Mutation", "R132H", "Not available"), title = "Gene mutation", legend_gp = gpar(fill = c("black", "#377EB8", "grey"))),
	Legend(at = 1:2, labels = c("Gain", "Loss"), title = "Chromosome abbreviation", legend_gp = gpar(fill = c("#E41A1C", "#377EB8"))),
	Legend(at = 1:2, labels = c("Amplification", "Deletion"), title = "Gene alteration", legend_gp = gpar(fill = c("#E41A1C", "#377EB8")))
)
figureS1_B = grid.grabExpr({
draw(hm, heatmap_legend_list = lgd_list, annotation_legend_side = "right", heatmap_legend_side = "right",
	padding = unit(c(0.2, 3, 0.2, 0.2), "cm"))
decorate_annotation("Age", {grid.text("Age", unit(-10, "mm"), just = "right")})
for(nm in names(ha.bottom@anno_list)) {
	decorate_annotation(nm, { grid.rect(gp = gpar(fill = NA, col = "black")) })
}
decorate_annotation("EGFR amp", {grid.lines(unit(c(-30, 0), "mm"), unit(c(1, 1), "npc"))})
decorate_annotation("IDH1/2", {grid.lines(unit(c(-30, 0), "mm"), unit(c(1, 1), "npc"))})
decorate_annotation("Chr7 gain", {grid.lines(unit(c(-30, 0), "mm"), unit(c(1, 1), "npc"))})
}, width = 12, height = 12)




pdf(qq("@{PROJECT_DIR}/image/figureS1.pdf"), width = 18, height = 16)
grid.newpage()
pushViewport(viewport(x = unit(1, "cm"), 
	width = unit(1, "npc") - unit(1, "cm"), height = 0.95,
	just = "left"))
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, width = c(2, 1.2))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, width = 0.95))
grid.draw(figureS1_B)
grid.text("A", x = unit(0, "mm"), y = 1, just = c("right", "top"), gp = gpar(fontsize = 30))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
pushViewport(viewport(y = 0.4, height = 0.6, just = "bottom"))
grid.echo(figureS1_C, newpage = FALSE)
grid.text("B", x = unit(0, "mm"), y = 1, just = c("right", "top"), gp = gpar(fontsize = 30))
popViewport()
pushViewport(viewport(x = unit(-1, "cm"), y = 0, height = 0.4, just = c("left", "bottom")))
grid.echo(figureS1_D, newpage = FALSE)
grid.text("C", x = unit(5, "mm"), y = 1, just = c("right", "top"), gp = gpar(fontsize = 30))
popViewport(4)
dev.off()

############### new summary of the cohort plot ###########
pdf(qq("~/figure1A.pdf"), width = 14, height = 8)

df = read.csv("/icgc/dkfzlsdf/analysis/hipo/hipo_016/manuscript/supplementary_tables/Supplementary Table 1 - Sample Metadata.csv", check.names = FALSE)
df = df[order(factor(df[, "subtype_final"], levels = c("IDH", "MES", "RTK_I", "RTK_II", "normal"))), ]

subtype_col = RColorBrewer::brewer.pal(5, "Set1")
names(subtype_col) = c("IDH", "MES", "RTK_I", "RTK_II", "normal")
gender_col = c("m"="#1B9E77", "f"="#E2798A")
seq_type_col = c("TRUE" = "#BF5B17", "FALSE" = NA)
mut_col = c("wt"="white", "mut"="black", "R132H"="#7570B3")
gain_loss_col = c("normal"="white", "gain"="#F0027F", "loss"="#386CB0")
amp_del_col = c("normal"="white", "amp"="#F0027F", "del"="#386CB0") 

ha <- ComplexHeatmap::HeatmapAnnotation( 
	Subtype = df$subtype_final,
	Age = anno_points(df$age, pch = 16, ylim = c(0, max(df[, "age"], na.rm = T)*1.1), axis = TRUE, axis_side = "left"),
    Gender = df$gender,

    "Epic/450K" = df$data_meth_array,
    "WGBS" = df$data_wgbs,
    "RNA-seq" = df$data_rnaseq,
    "ChIP-seq" = df$data_chipseq_H3K27ac | df$data_chipseq_H3K27me3 | df$data_chipseq_H3K36me3 | df$data_chipseq_H3K4me1 | df$data_chipseq_H3K4me3 | df$data_chipseq_H3K9me3,
    "WGS" = df$data_chipseq_H3K27ac | df$data_chipseq_H3K27me3 | df$data_chipseq_H3K36me3 | df$data_chipseq_H3K4me1 | df$data_chipseq_H3K4me3 | df$data_chipseq_H3K9me3,
    "IDH1/2" = df$mutIDH1_2,
    "Chr7 gain" = df$chr7gain,
    "Chr10 loss" = df$chr10loss,
    "Chr19 gain" = df$chr19gain,
    "Chr20 gain" = df$chr20gain,
    "EGFR amp" = df$ampEGFR,
    "PDGFRA amp" = df$ampPDGFRA,
    "CDKN2A/B del" = df$delCDKN2A_B,
    "CDK4 amp" = df$ampCDK4,
    "RB1 del" = df$delRB1,
    "MET amp" = df$ampMET,
    "MDM2 amp" = df$MDM2,
    "MDM4 amp" = df$MDM4,
    "PTEN del" = gsub(" ", "", df$delPTEN),
    show_legend = FALSE,
    col=list(
    	Subtype = subtype_col,
		Gender= gender_col,
		"Epic/450K" = seq_type_col,
		"WGBS" = seq_type_col,
		"RNA-seq" = seq_type_col,
		"ChIP-seq" = seq_type_col,
		"WGS" = seq_type_col,
		"Chr7 gain"= gain_loss_col,
		"Chr10 loss"=gain_loss_col,
		"Chr19 gain"=gain_loss_col,	
		"Chr20 gain"=gain_loss_col,	
		"IDH1/2"= mut_col,
		"EGFR amp"= amp_del_col,		
		"PDGFRA amp"= amp_del_col,		
		"CDKN2A/B del"= amp_del_col,		
		"CDK4 amp"= amp_del_col,		
		"RB1 del"= amp_del_col,	
		"MET amp"= amp_del_col, 
		"MDM2 amp"= amp_del_col,
		"MDM4 amp"= amp_del_col,
		"PTEN del"= amp_del_col 
	),
    na_col="grey", show_annotation_name= c(TRUE, FALSE, rep(TRUE, 24)), 
    annotation_name_side = "left",
    annotation_height = unit(c(0.5, 2, rep(0.5, 24)), "cm")
 )

 
hm <- Heatmap(matrix(nrow = 0, ncol = nrow(df)), top_annotation = ha)

lgd_list = list(
	Legend(at = 1:5, labels = names(subtype_col), title = "Subtype", legend_gp = gpar(fill = subtype_col)),
	Legend(at = 1:2, labels = c("Male", "Female"), title = "Gender", legend_gp = gpar(fill = gender_col)),
	Legend(at = 1, labels = "Available", title = "Sequence data", legend_gp = gpar(fill = seq_type_col[1])),
	Legend(at = 1:3, labels = c("WT", "Mutation", "R132H"), title = "Gene mutation", legend_gp = gpar(fill = mut_col)),
	Legend(at = 1:2, labels = c("Gain", "Loss"), title = "Chromosome abbreviation", legend_gp = gpar(fill = gain_loss_col[-1])),
	Legend(at = 1:2, labels = c("Amplification", "Deletion"), title = "Gene alteration", legend_gp = gpar(fill = amp_del_col[-1])),
	Legend(at = 1, labels = "Infomation not available", legend_gp = gpar(fill = "grey"))
)

grid.newpage()
pushViewport(viewport(width = 0.8))
draw(hm, heatmap_legend_list = lgd_list, annotation_legend_side = "right", heatmap_legend_side = "right", newpage = FALSE)
upViewport()

decorate_annotation("Age", {grid.text("Age", unit(-10, "mm"), just = "right")})
for(nm in names(ha@anno_list)) {
	decorate_annotation(nm, { grid.rect(gp = gpar(fill = NA, col = "black")) })
}
decorate_annotation("EGFR amp", {grid.lines(unit(c(-30, 0), "mm"), unit(c(1, 1), "npc"))})
decorate_annotation("IDH1/2", {grid.lines(unit(c(-30, 0), "mm"), unit(c(1, 1), "npc"))})
decorate_annotation("Chr7 gain", {grid.lines(unit(c(-30, 0), "mm"), unit(c(1, 1), "npc"))})
decorate_annotation("Epic/450K", {grid.lines(unit(c(-30, 0), "mm"), unit(c(1, 1), "npc"))})
dev.off()



######## compare top 1kb windows and probes
library(circlize)
library(EnrichedHeatmap)
chrominfo = read.chromInfo(chromosome.index = CHROMOSOME)$df
chrom_gr = GRanges(chrominfo[, 1], ranges = IRanges(chrominfo[, 2]+1, chrominfo[, 3]))
chrom_window = makeWindows(chrom_gr, w = 1000)
chrom_window = get_mean_methylation_in_genomic_features(rownames(SAMPLE)[SAMPLE[, 1] != "normal"], chrom_window)
chrom_window = chrom_window[chrom_window$ncpg > 5]

chrom_window$sd = rowSds(as.matrix(as.data.frame(mcols(chrom_window)[, 3:62])))

probes = scan("/icgc/dkfzlsdf/analysis/hipo/hipo_016/cluster_450k/20170118_volker_hccpgs_order.txt", what = "character")

suppressMessages( library(IlluminaHumanMethylation450kanno.ilmn12.hg19) )

data(Locations)

probes_df = Locations[probes, ]

probes_gr = GRanges(seqnames = probes_df[, 1], ranges = IRanges(probes_df[, 2], probes_df[, 2] + 1))

mtch = as.matrix(findOverlaps(probes_gr, chrom_window))
chrom_window$has_probe = 0
chrom_window[mtch[, 2]]$has_probe = 1


png(qq("@{PROJECT_DIR}/image/probes_WGBS_window_overlap.png"), height = 1000)
x = chrom_window$sd
od = order(x)
x2 = x[od]
l2 = chrom_window$has_probe[od] == 1
par(mfrow = c(2, 1))
plot(seq_along(x2), x2, type = "l", xlim = c(1, length(x2)),
	xlab = "Index", ylab = "SD", main = "SD of 1kb-window methylation from WGBS")
tb = table(round(which(l2), -4))
plot(as.numeric(names(tb)), tb, type = "h", xlim = c(1, length(x2)),
	xlab = "Index", ylab = "Number of windows", main = "1kb-windows that overlap to top 8k most variable probes")
axis(side = 2)
dev.off()


