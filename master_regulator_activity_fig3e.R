## title: "GRN analysis HIPO16"
## author: "Carl Herrmann"
## date: "30/07/2019"

library(pheatmap)
library(RColorBrewer)

## Define the main dir
main.dir='/icgc/dkfzlsdf/analysis/hipo/hipo_016/user_folder/herrmanc/analysis/2019-04-13_VIPER/'
setwd(main.dir)

## load precomputed datasets
crctf = read.table('data/CRC_TF.txt',as.is=TRUE)[,1] # these are the TFs from the CRC analysis
tf.activity.test = readRDS('results/tf.activity.test.rds') # VIPER TF activity based on test (i.e. TCGA) expression datasets
meta = read.table('data/meta_data_samples.csv',header=TRUE,row.names=1,sep=',') # Annotation of the samples
#
subtype = meta[,1,drop=FALSE]
subtype[,1] = factor(as.character(subtype[,1]),levels=c('IDH','MES','RTK_I','RTK_II','normal'))
sub.col = brewer.pal(5,'Set1')
names(sub.col) = levels(subtype$subtype_final)


## Select only the rows corresponding to CRC TFs
tf.act.test.crctf.h16 = tf.activity.test[rownames(tf.activity.test) %in% crctf,]


## Only a subset of the CRC TFs (38) is found in the VIPER activity table:
dim(tf.act.test.crctf.h16)


i = match(rownames(subtype),colnames(tf.act.test.crctf.h16))
i = i[!is.na(i)]
tf.act.test.crctf.h16 = tf.act.test.crctf.h16[,i]

## Computing the per subtype TF activity:
  
tfa.subtype.crctf = sapply(c('IDH','MES','RTK_I','RTK_II'),function(type) {
  s = rownames(subtype)[which(subtype[,1]==type)]
  i = which(colnames(tf.act.test.crctf.h16) %in% s)
  apply(tf.act.test.crctf.h16,1,function(x) {mean(x[i])})
})
tfa.subtype.crctf.norm = t(scale(t(tfa.subtype.crctf)))

### TF activity per subtype

XX = as.data.frame(tfa.subtype.crctf.norm)
col.sub = brewer.pal(4,'Set1')
names(col.sub) = levels(subtype[,1])[1:4]
i.max = apply(XX,1,function(x) {which(x==max(x))})
XX = XX[order(i.max),]
pheatmap(t(XX),
         col=rev(colorRampPalette(brewer.pal(8,'Spectral'))(100)),
         annotation_col = data.frame(subtype=i.max),
         annotation_colors = list(subtype=col.sub),
         cluster_rows = FALSE,
         cluster_cols = FALSE)

