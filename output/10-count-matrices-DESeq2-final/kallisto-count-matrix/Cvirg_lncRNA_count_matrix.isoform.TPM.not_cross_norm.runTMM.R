library(edgeR)

rnaseqMatrix = read.table("/home/shared/8TB_HDD_02/zbengt/github/oyster-lnc/output/10-count-matrices-DESeq2-final/kallisto-count-matrix/Cvirg_lncRNA_count_matrix.isoform.TPM.not_cross_norm", header=T, row.names=1, com='', check.names=F)
rnaseqMatrix = as.matrix(rnaseqMatrix)
rnaseqMatrix = round(rnaseqMatrix)
exp_study = DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))
exp_study = calcNormFactors(exp_study)
exp_study$samples$eff.lib.size = exp_study$samples$lib.size * exp_study$samples$norm.factors
write.table(exp_study$samples, file="/home/shared/8TB_HDD_02/zbengt/github/oyster-lnc/output/10-count-matrices-DESeq2-final/kallisto-count-matrix/Cvirg_lncRNA_count_matrix.isoform.TPM.not_cross_norm.TMM_info.txt", quote=F, sep="\t", row.names=F)
