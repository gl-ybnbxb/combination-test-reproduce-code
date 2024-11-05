library(MetaCycle)
library(gap)
library(rmutil)

source('data_scripts/runJTK.R') # Original JTK method
source('data_scripts/runJTK_combination.R') # Modified with combination test
source('data_scripts/transformation_test_new.R')

# Load the circadian gene expression data
hughes=read.csv('data_scripts/hughes.csv')
cycMouseLiverRNA = hughes
rm(hughes)

# Positive and negative controls
load('data_scripts/gene.symbl.p.rda')
load('data_scripts/gene.symbl.n.rda')
length(gene.symbl.p); length(gene.symbl.n)
gene.symbl.n=gene.symbl.n[gene.symbl.n%in%cycMouseLiverRNA$hughes.anno...4.] 
gene.symbl.p=gene.symbl.p[gene.symbl.p%in%cycMouseLiverRNA$hughes.anno...4.]
length(gene.symbl.p); length(gene.symbl.n)

# Focus on the positive and negative control genes (comment out to run on genome wide)
cycMouseLiverRNA=cycMouseLiverRNA[match(c(gene.symbl.p, gene.symbl.n), cycMouseLiverRNA$hughes.anno...4.),]

JTK.results=as.data.frame(runJTK(cycMouseLiverRNA, JTKtime = 18:65))
output=JTK.results[,c(1,3)]
colnames(output)[2]='JTK'
output[,2]=as.numeric(output[,2])

for(test in c("Fisher", "Bonferroni", "Cauchy", "truncated_0.9",
              "Pareto", "Frechet", "Levy")){
  cat(test,'\n')
  JTK.results=as.data.frame(runJTK_combination(cycMouseLiverRNA, JTKtime = 18:65, method=test))
  output=cbind(output, as.numeric(JTK.results$ADJ.P))
  colnames(output)[ncol(output)]=test
  cat('\n')
}
head(output)



output$JTK=as.numeric(output$JTK)
plot(output$JTK, output$Bonferroni) # JTK implements Bonferroni

test.output=output[,c('Bonferroni','Levy','Cauchy','truncated_0.9','Frechet','Pareto','Fisher')]
rownames(test.output)=output$CycID

c.col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
par(mfrow=c(1,2))
par(mar=c(8.1, 4.1, 4.1, 2.1))
boxplot(((test.output[match(gene.symbl.p, rownames(test.output)),])), 
        ylab='p-value', main='Positive control genes',las=2, 
        col=c.col, ylim=c(0,0.005), pch=16, cex=0.5)
boxplot(((test.output[match(gene.symbl.n, rownames(test.output)),])), 
        ylab='p-value', main='Negative control genes', las=2,
        col=c.col, ylim=c(0,1), pch=16, cex=0.5)
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))


