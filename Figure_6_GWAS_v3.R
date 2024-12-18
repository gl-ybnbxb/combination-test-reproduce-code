library(pheatmap)
library(rmutil)
source('data_scripts/transformation_test_new.R')

load('data_scripts/gwas.rda')

output=matrix(nrow=nrow(G_list), ncol=10)
for(i in 1:nrow(output)){
  if(i%%100==0) cat(i,'\t')
  gene=G_list[i,2]
  gene.id=G_list[i,1]
  # SNP-level p-values
  pvals=gwas.calls[[gene.id]]$gwasP
  pvals.output=c(N.snps=length(pvals),
                 Epic=signif(epic.calls[i],3),  # EPIC p-value
                 Magma=signif(magma.calls[i],3),  # MAGMA p-value
                 Bonferroni=signif(transformation_global_test(pvals, method='Bonferroni'),3),
                 Levy=signif(transformation_global_test(pvals, method='Levy'),3),
                 Cauchy=signif(transformation_global_test(pvals, method='Cauchy'),3),
                 truncated_0.9=signif(transformation_global_test(pvals, method='truncated_0.9'),3),
                 Frechet=signif(transformation_global_test(pvals, method='Frechet'),3),
                 Pareto=signif(transformation_global_test(pvals, method='Pareto'),3),
                 Fisher=signif(transformation_global_test(pvals, method='Fisher'),3))
  output[i,]=pvals.output
  colnames(output)=names(pvals.output)
}
output=as.data.frame(output)
output=cbind(gene.ensembl=G_list[,1],
             gene=G_list[,2],
             output)
summary(output$N.snps)

cbPalette <- c("gray30", "gray80", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Call significant genes after FDR control
# And compare the number of significant genes across methods
for(fdr.thresh in c(0.05, 0.2)){
  overlap.mat=matrix(nrow=ncol(output)-3, ncol=ncol(output)-3)
  rownames(overlap.mat)=colnames(overlap.mat)=colnames(output)[-(1:3)]
  for(i in 1:nrow(overlap.mat)){
    for(j in i:ncol(overlap.mat)){
      i.rej= which(p.adjust(output[,rownames(overlap.mat)[i]], method='fdr') < fdr.thresh)
      j.rej= which(p.adjust(output[,colnames(overlap.mat)[j]], method='fdr') < fdr.thresh)
      overlap.mat[i,j]=length(intersect(i.rej, j.rej))
    }
  }
  overlap.mat.txt=(overlap.mat)
  overlap.mat.txt[is.na(overlap.mat.txt)]=''
  pdf(file=paste0('figure/pval_sig_genes_SCZ_',fdr.thresh,'_v3.pdf'), width=6, height=6)
  pheatmap(log(overlap.mat), cluster_rows=F, cluster_cols=F, na_col="white", 
           display_numbers = overlap.mat.txt, fontsize_number=14,
           legend = FALSE, main=paste('No. of sig. genes after FDR control:',fdr.thresh,'SCZ'))
  dev.off()
}


# For those whose numbers of SNPs are smaller than 50.
# Call significant genes after FDR control
# And compare the number of significant genes across methods
output = output[output$N.snps<=50,]
for(fdr.thresh in c(0.05, 0.2)){
  overlap.mat=matrix(nrow=ncol(output)-3, ncol=ncol(output)-3)
  rownames(overlap.mat)=colnames(overlap.mat)=colnames(output)[-(1:3)]
  for(i in 1:nrow(overlap.mat)){
    for(j in i:ncol(overlap.mat)){
      i.rej= which(p.adjust(output[,rownames(overlap.mat)[i]], method='fdr') < fdr.thresh)
      j.rej= which(p.adjust(output[,colnames(overlap.mat)[j]], method='fdr') < fdr.thresh)
      overlap.mat[i,j]=length(intersect(i.rej, j.rej))
    }
  }
  overlap.mat.txt=(overlap.mat)
  overlap.mat.txt[is.na(overlap.mat.txt)]=''
  pdf(file=paste0('figure/less_snps_pval_sig_genes_SCZ_',fdr.thresh,'_v3.pdf'), width=6, height=6)
  pheatmap(log(overlap.mat), cluster_rows=F, cluster_cols=F, na_col="white", 
           display_numbers = overlap.mat.txt, fontsize_number=14,
           legend = FALSE, main=paste('No. of sig. genes after FDR control:',fdr.thresh,'SCZ'))
  dev.off()
}
