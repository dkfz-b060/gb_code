#### Figure 2D ######
library(GetoptLong)


load_namespace("/desktop-home/guz/project/development/epik")
load_epik_config("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/script_publicate/hipo16_config_epik.R")


# lmr_list, dmv_list, pmd_cluster
load(qq("@{PROJECT_DIR}/rds/fastseg_lmr_dmv_pmd_segmentations_new_settings.rds"))

pmd_cluster = pmd_cluster[rownames(SAMPLE)]
dmv_list = dmv_list[rownames(SAMPLE)]
lmr_list = lmr_list[rownames(SAMPLE)]

# load(qq("@{PROJECT_DIR}/rds/block.RData"))
pmd_cluster = list()
for(sid in rownames(SAMPLE)) {
    pmd_cluster[[sid]] = GRanges()
    for(chr in CHROMOSOME) {
        qqcat("@{sid}, @{chr}...\n")
        pmd_cluster[[sid]] = c(pmd_cluster[[sid]], readRDS(qq("@{PROJECT_DIR}/rds/pmd_methylseekr_@{chr}_@{sid}.rds")))
    }
    pmd_cluster[[sid]] = pmd_cluster[[sid]][ pmd_cluster[[sid]]$type == "PMD" ]
    # pmd_cluster[[sid]] = intersect(pmd_cluster[[sid]], bb)
    pmd_cluster[[sid]] = pmd_cluster[[sid]][ width(pmd_cluster[[sid]]) > 10000 ] 
}

####### basic statistics ########
library(grid)

###### upset plot ##########
dmv_consensus = tapply(seq_along(dmv_list), SAMPLE$subgroup, function(index) {
    gr = common_regions(dmv_list[index], min_recurrency = ifelse(length(index) > 4, 1/3*length(index), 4))
    seqlengths(gr) = NA
    gr
})[c("IDH", "MES", "RTK_I", "RTK_II")]
lmr_consensus = tapply(seq_along(lmr_list), SAMPLE$subgroup, function(index) {
    gr = common_regions(lmr_list[index], min_recurrency = ifelse(length(index) > 4, 1/3*length(index), 4))
    seqlengths(gr) = NA
    gr
})[c("IDH", "MES", "RTK_I", "RTK_II")]
pmd_consensus = tapply(seq_along(pmd_cluster), SAMPLE$subgroup, function(index) {
    gr = common_regions(pmd_cluster[index], min_recurrency = ifelse(length(index) > 4, 1/3*length(index), 4))
    seqlengths(gr) = NA
    gr
})[c("IDH", "MES", "RTK_I", "RTK_II")]


###### enrich to genomic features and chromatin states ####

map = c("E01" = "01_TssActive",
        "E02" = "01_TssActive",
        "E03" = "01_TssActive",
        "E04" = "01_TssActive",
        "E05" = "02_Transcript",
        "E06" = "02_Transcript",
        "E07" = "03_Enhancer",
        "E08" = "03_Enhancer",
        "E09" = "03_Enhancer",
        "E10" = "03_Enhancer",
        "E11" = "03_Enhancer",
        "E12" = "04_Heterochromatin",
        "E13" = "04_Heterochromatin",
        "E14" = "05_TssBiv",
        "E15" = "03_Enhancer",
        "E16" = "06_Repressive",
        "E17" = "06_Repressive",
        "E18" = "07_Quies")
state_col = c("01_TssActive" = "#FF0000", "02_Transcript" = "#008000", "03_Enhancer" = "#C2E105", 
    "04_Heterochromatin" = "#8A91D0", "05_TssBiv" = "#CD5C5C",
    "06_Repressive" = "#808080", "07_Quies" = "#000000")
old_fun = chipseq_hooks$chromHMM
chipseq_hooks$chromHMM = function(sid, chr = CHROMOSOME) {
    gr = old_fun(sid, chr)
    if(length(gr)) {
        gr$states = map[gr$states]
        gr
    } else {
        gr
    }
}

states_normal = get_chromHMM_list(paste0("E0", 67:74), merge = TRUE)

chrom_list = lapply(c("IDH", "MES", "RTK_I", "RTK_II"), function(subgroup) {
    sample_id = rownames(SAMPLE[SAMPLE$subgroup == subgroup, , drop = FALSE])
    get_chromHMM_list(sample_id, merge = TRUE)
})
names(chrom_list) = c("IDH", "MES", "RTK_I", "RTK_II")

pmd_consensus = tapply(seq_along(pmd_cluster), SAMPLE$subgroup, function(index) {
    gr = common_regions(pmd_cluster[index], min_recurrency = ifelse(length(index) > 4, 1/3*length(index), 4))
    seqlengths(gr) = NA
    gr[, NULL]
})[c("IDH", "MES", "RTK_I", "RTK_II")]; pmd_consensus = lapply(pmd_consensus, function(x) x)
lmr_consensus = tapply(seq_along(lmr_list), SAMPLE$subgroup, function(index) {
    gr = common_regions(lmr_list[index], min_recurrency = ifelse(length(index) > 4, 1/3*length(index), 4))
    seqlengths(gr) = NA
    gr[, NULL]
})[c("IDH", "MES", "RTK_I", "RTK_II")]; lmr_consensus = lapply(lmr_consensus, function(x) x)
dmv_consensus = tapply(seq_along(dmv_list), SAMPLE$subgroup, function(index) {
    gr = common_regions(dmv_list[index], min_recurrency = ifelse(length(index) > 4, 1/3*length(index), 4))
    seqlengths(gr) = NA
    gr[, NULL]
})[c("IDH", "MES", "RTK_I", "RTK_II")]; dmv_consensus = lapply(dmv_consensus, function(x) x)

for(i in seq_along(pmd_cluster)) {
    pmd_cluster[[i]] = intersect(pmd_cluster[[i]], gr_fastseg_list[[i]][gr_fastseg_list[[i]]$states == 2])
}

for(i in seq_along(pmd_consensus)) {
    pmd_consensus[[i]] = setdiff(pmd_consensus[[i]], dmv_consensus[[i]])
}

names(pmd_consensus) = paste0("PMD_", names(pmd_consensus))
names(lmr_consensus) = paste0("LMR_", names(lmr_consensus))
names(dmv_consensus) = paste0("DMV_", names(dmv_consensus))

lt = c(pmd_consensus, lmr_consensus, dmv_consensus)

##
for(nm in names(lt)) {
    print(nm)
    subtype = gsub("PMD_|LMR_|DMV_", "", nm)
    if(!grepl("PMD", nm)) next
    mr = lt[[nm]][, NULL]
    mr = get_mean_methylation_in_genomic_features(rownames(SAMPLE[SAMPLE$subgroup %in% subtype, ,drop = FALSE]), mr)
    mr = annotate_to_genomic_features(mr, as.list(split(chrom_list[[subtype]], chrom_list[[subtype]]$states)))
    mr = annotate_to_genomic_features(mr, GENOMIC_FEATURE_LIST[c("cgi", "cgi_shore")])
    mr = annotate_to_gene_models(mr, TXDB)
    lt[[nm]] = mr
}

save(lt, file = qq("@{PROJECT_DIR}/rds/md.RData"))

mean_meth = lapply(lt, function(dmr) {
    mean_meth = mcols(dmr)
    mean_meth = as.matrix(mean_meth[, grepl("mean_meth", colnames(mean_meth))])
    # mean_meth = cbind(rowMeans(mean_meth[, grepl("AK", colnames(mean_meth))]))
    # colMeans(mean_meth)
    mean_meth
})
# colnames(mean_meth) = c("tumor")
n_gr = unlist(lapply(lt, length))
w_gr = sapply(lt, function(x) median(width(x)))
gene_anno = do.call("rbind", lapply(lt, function(dmr) {
        gr = dmr
        c(sum(gr$overlap_to_gene*width(gr))/sum(width(gr)),
          sum(gr$overlap_to_intergenic*width(gr))/sum(width(gr)))
}))
dist_tss = do.call("rbind", lapply(lt, function(dmr) {
        gr = dmr
        x = abs(gr$dist_to_gene_tss)
        c(sum(x < 1000)/length(x),
          sum(x >= 1000 & x < 5000)/length(x),
          sum(x >= 5000 & x < 10000)/length(x),
          sum(x >= 10000)/length(x))
}))
cgi_anno = do.call("rbind", lapply(lt, function(dmr) {
        gr = dmr
        c(sum(gr$overlap_to_cgi*width(gr))/sum(width(gr)),
          sum(gr$overlap_to_cgi_shore*width(gr))/sum(width(gr)))
}))
states_tumor = do.call("rbind", lapply(lt, function(dmr) {
        gr = dmr
        x = c(sum(gr$overlap_to_01_TssActive*width(gr))/sum(width(gr)),
          sum(gr$overlap_to_02_Transcript*width(gr))/sum(width(gr)),
          sum(gr$overlap_to_03_Enhancer*width(gr))/sum(width(gr)),
          sum(gr$overlap_to_04_Heterochromatin*width(gr))/sum(width(gr)),
          sum(gr$overlap_to_05_TssBiv*width(gr))/sum(width(gr)),
          sum(gr$overlap_to_06_Repressive*width(gr))/sum(width(gr)),
          sum(gr$overlap_to_07_Quies*width(gr))/sum(width(gr))
          )
        x/sum(x)
}))

# res1 = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/rds/genomic_regions_correlation_consensus_pmd_vs_gf.rds")
res1 = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/rds/genomic_regions_correlation_consensus_methylseekr_pmd_vs_gf.rds")
res2 = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/rds/genomic_regions_correlation_consensus_lmr_vs_gf.rds")
res3 = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/rds/genomic_regions_correlation_consensus_dmv_vs_gf.rds")
colnames(res1$z) = paste0("PMD_", colnames(res1$z))
colnames(res2$z) = paste0("LMR_", colnames(res2$z))
colnames(res3$z) = paste0("DMV_", colnames(res3$z))
mat_enrich_gf = t(cbind(res1$z, res2$z, res3$z))[names(lt), ]
mat_enrich_gf = mat_enrich_gf[, c("gene", "tss", "exon", "intergenic", "cgi", "cgi_shore", "TFBS", "LINE", "SINE")]
# res1 = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/rds/genomic_regions_correlation_consensus_pmd_vs_simplified_states.rds")
res1 = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/rds/genomic_regions_correlation_consensus_methylseekr_pmd_vs_simplified_states.rds")
res1 = do.call("cbind", lapply(res1, function(x) x$z))
colnames(res1) = paste0("PMD_", colnames(res1))
res2 = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/rds/genomic_regions_correlation_consensus_lmr_vs_simplified_states.rds")
res2 = do.call("cbind", lapply(res2, function(x) x$z))
colnames(res2) = paste0("LMR_", colnames(res2))
res3 = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final_cohort/rds/genomic_regions_correlation_consensus_dmv_vs_simplified_states.rds")
res3 = do.call("cbind", lapply(res3, function(x) x$z))
colnames(res3) = paste0("DMV_", colnames(res3))
mat_enrich_st = t(cbind(res1, res2, res3))[names(lt), ]
colnames(mat_enrich_st) = gsub("^\\d\\d_", "", colnames(mat_enrich_st))

library(ComplexHeatmap)
library(circlize)
###### make a simplified summary plot #########
state_col = c("01_TssActive" = "#FF0000", "02_Transcript" = "#008000", "03_Enhancer" = "#C2E105", 
    "04_Heterochromatin" = "#8A91D0", "05_TssBiv" = "#CD5C5C",
    "06_Repressive" = "#808080", "07_Quies" = "#000000")
split = rep(names(lt))
split_df = data.frame(direction = gsub("_.*$", "", split),
    type = gsub("^(PMD_|LMR_|DMV_)", "", split))
levels(split_df$type) = c("IDH", "MES", "RTK_I", "RTK_II")

pdf(qq("@{PROJECT_DIR}/image/md_annotated_simplified.pdf"), width = 20, height = 6)
# mean methylation
z_score_col_fun = colorRamp2(c(-200, 0, 200), c("green", "white", "red"))
labels = paste(split_df$type, split_df$direction, sep = "_")
ht_list = rowAnnotation(text = row_anno_text(labels, offset = unit(1, "npc"), just = "right"), width = max_text_width(labels))
ht_list = ht_list + 
rowAnnotation("Mean methylation" = row_anno_boxplot(mean_meth, axis = TRUE, ylim = c(0, 1), outline = FALSE), width = unit(4, "cm")) +
rowAnnotation("Number of regions" = row_anno_barplot(n_gr, bar_width = 1, axis = TRUE, baseline = 0), width = unit(4, "cm")) +
rowAnnotation("Median width" = row_anno_barplot(w_gr, bar_width = 1, axis = TRUE, baseline = 0), width = unit(4, "cm")) +
rowAnnotation("Distance to TSS" = row_anno_barplot(dist_tss, bar_width = 1, axis = TRUE, gp = gpar(fill = c("#FF0000FF", "#FF7352FF", "#FFB299FF", "#FFD9CBFF"))), width = unit(4, "cm")) +
rowAnnotation("Gene annotation" = row_anno_barplot(gene_anno, bar_width = 1, axis = TRUE, gp = gpar(fill = c("green", "blue"))), width = unit(4, "cm")) +
rowAnnotation("CGI annotation" = row_anno_barplot(cgi_anno, bar_width = 1, axis = TRUE, gp = gpar(fill = c("#FFA500FF", "#FFD191FF"))), width = unit(4, "cm")) +
Heatmap(mat_enrich_gf, name = "Enrichment to\ngenomic features", col = z_score_col_fun, cluster_columns = FALSE, show_row_names = FALSE,
    width = unit(ncol(mat_enrich_gf)*7, "mm"), column_title = "", combined_name_fun = NULL, show_heatmap_legend = FALSE) +
rowAnnotation("Overlap to\nchromatin states" = row_anno_barplot(states_tumor, bar_width = 1, axis = TRUE, gp = gpar(fill = state_col)), width = unit(4, "cm")) +
Heatmap(mat_enrich_st, name = "Enrichment to\nchromatin states", col = z_score_col_fun, cluster_columns = FALSE, show_row_names = FALSE,
    width = unit(ncol(mat_enrich_st)*7, "mm"), column_title = "", combined_name_fun = NULL,
    column_names_gp = gpar(col = state_col), show_heatmap_legend = FALSE)

lgd_list = list(Legend(at = 2:3, labels = c("gene", "intergenic"), title = "Gene annotation", legend_gp = gpar(fill = c("green", "blue"))),
    Legend(at = 1:4, labels = c("<1kb", "1kb~5kb", "5kb~10kb", ">10kb"), title = "Distance to TSS", legend_gp = gpar(fill = c("#FF0000FF", "#FF7352FF", "#FFB299FF", "#FFD9CBFF"))),
    Legend(at = 1:2, labels = c("CGI", "CGI shore"), title = "CGI annotation", legend_gp = gpar(fill = c("#FFA500FF", "#FFD191FF"))),
    Legend(col_fun = z_score_col_fun, at = c(-200, -100, 0, 100, 200), title = "Z-score"),
    Legend(at = 1:9, labels = gsub("\\d\\d_", "", names(state_col)), title = "Chromatin states", legend_gp = gpar(fill = state_col)))
figureS4_C = grid.grabExpr({
draw(ht_list, padding = unit(c(1, 1, 1, 1), "cm"), split = split_df$direction, cluster_rows = FALSE,
    heatmap_legend_list = lgd_list)
for(an in c("Mean methylation", "Number of regions", "Median width", "Gene annotation", "Distance to TSS", "CGI annotation", "Overlap to\nchromatin states")) {
    decorate_annotation(an, {
        grid.text(an, y = unit(1, "npc") + unit(3, "mm"), just = "bottom")
    })
}
for(an in c("Enrichment to\ngenomic features", "Enrichment to\nchromatin states")) {
    decorate_heatmap_body(an, {
        grid.text(an, y = unit(1, "npc") + unit(3, "mm"), just = "bottom")
    })
}
})
dev.off()

