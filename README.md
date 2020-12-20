# Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype


Yonghe Wu\*, Michael Fletcher\*, Zuguang Gu\*, Qi Wang, Barbara Costa, Anna Bertoni, Ka-Hou Man, Magdalena Schlotter, Jörg Felsberg, Jasmin Mangei, Martje Barbus, Ann-Christin Gaupel, Wei Wang, Tobias Weiss, Roland Eils, Michael Weller, Haikun Liu, Guido Reifenberger, Andrey Korshunov, Peter Angel, Peter Lichter, Carl Herrmann$, Bernhard Radlwimmer$

*\* co-first authors*

*$ co-corresponding authors*

#### Published online (2020-12-18) at Nature Communications: https://doi.org/10.1038/s41467-020-20225-w

## Introduction

Hello! Welcome to the code repository for our paper.

If you've seen the paper, or looked at our GEO/EGA accessions, you'll realise there's a *lot* of data/analysis - please ask us if you have any issues with, want to discuss, or need help troubleshooting the code.

* WGBS, RNAseq data processing; methylation analysis (PMD/LMR/DMV feature calling, basic statistics, etc. as in Figs 1-2) was handled by Zuguang Gu [email: z.gu [AT] dkfz-heidelberg.de]
* ChIPseq data processing; tumour/cell line ChromHMM, tumour CRC analysis (as in Fig 3) was done by Carl Hermann [email: carl.herrmann [AT] bioquant.uni-heidelberg.de ] and Qi Wang
* Mouse RNAseq data processing; RTN Master Regulator (SFig 3), RNAseq subtype genes (SFig 2), tumour methylation microarrays (Fig 7, SFig 1), cell line gene expression microarrays analyses (SFig 5), cell line ATAC/ChIPseq analysis (Fig 4, SFig 6) and the re-analysis of the Darmanis (2017) scRNAseq data (as in Fig 4) were performed by Mike Fletcher [email: m.fletcher [AT] dkfz-heidelberg.de OR sci [AT] dismissed.net.nz]

For general enquiries, please contact the corresponding authors: Carl (email above) or Bernhard Radlwimmer [email: b.radlwimmer [AT] dkfz-heidelberg.de]

## WashU Epigenome Browser Datahub

![WashU Epigenome Browser view of RTK_I tracks at the RTK_I MR MYT1](https://github.com/dkfz-b060/gb_code/blob/master/WashU_MYT1_RTK_I.jpg)
*The complete set of tracks for the RTK_I subtype at the MYT1 gene (as seen in Figure 1); from top to bottom: H3K27ac ChIPseq; subtype superenhancers; H3K4me3, H3K4me1, H3K36me3, H3K27me3, H3K9me3 ChIPseq; RNAseq; methylation (WGBS); MACS2 narrow peak calls for histone mark ChIPseq; subtype methylation features (PMDs, DMVs, LMRs); MACS2 broad peak calls for the histone mark ChIPseq; ChromHMM*

We've created a .json Datahub file for the WashU Epigenome browser: http://epigenomegateway.wustl.edu/browser/

Instructions on how to load the Datahub are available at https://epigenomegateway.readthedocs.io/en/latest/usage.html#adding-a-custom-track-or-data-hub

Information on the format of the Datahub .json (if you'd like to edit it yourself, say) can be found at https://epigenomegateway.readthedocs.io/en/latest/datahub.html

**Get the Datahub**: https://github.com/dkfz-b060/gb_code/blob/master/WashU_Browser_Datahub_all_tracks.json

## Data deposition

**Raw, protected (patient sample)** data has been deposited at the **EGA** at https://www.ebi.ac.uk/ega/studies/EGAS00001003953

**Processed** data is **available at the GEO**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121723 (SuperSeries)

The subseries consist of:

Series Accession|Series Name
----------------|------------
GSE121716|Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype - cell line ChIPseq experiments
GSE121717|Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype - cell line RNAseq experiments
GSE121718|Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype - cell line gene expression microarray data
GSE121719|Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype - tumour ChIPseq data
GSE121720|Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype - tumour rRNA-depleted total ssRNAseq data
GSE121721|Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype - tumour methylation (WGBS) data
GSE121722|Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype - tumour methylation microarray data
GSE133040|Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype - cell line ATACseq
GSE145556|Glioblastoma epigenome profiling identifies SOX10 as a master regulator of molecular tumour subtype - mouse RNAseq experiments

For each sample/data type, you can get various outputs, depending on the data type/analysis: aligned .bams, peak calls, superenhancer analysis results, methylation/gene expression tables, etc. etc. etc. Check it out.
