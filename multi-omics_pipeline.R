library(FactoMineR)
library(OmicCircos)
library(mixOmics)

# Load data
load("pollack.RData")
head(pollack.nox)
names(pollack.nox)

## Extract omic layers
cnv <- pollack.nox[, 8:48]
expr <- pollack.nox[, 49:89]

### FactoMinerR analysis----

## Copy Number Variant PCS
head(cnv)
cnv.t<-t(cnv)
# assign names, we include a exp suffix to differentiate genes from expression
colnames(cnv.t)<-paste(pollack.nox$gene,"CN",sep=".")

# Clasify the samples by origin
cond <- c(rep("cell_line", 4), rep("NORWAY", 28), rep("STANFORD", 9))

# Construct data.frame to perform PCA
cnv4pca<-data.frame(cond,cnv.t)
res.pca.cnv<-PCA(cnv4pca,quali.sup=1)
res.pca.cnv

png("results/pca.png", width = 800, height = 600)
plot(res.pca.cnv,habillage=1)
dev.off()

## Expression PCA
colnames(expr) <- colnames(cnv) #To perform later MFA, we need to have the same names
head(expr)

expr.t<-t(expr)
# Assign names and include a exp suffix to differentiate genes from cnv
colnames(expr.t)<-paste(pollack.nox$gene,"exp",sep=".")

# Construct data.frame to perform PCA
expr4pca<-data.frame(cond,expr.t)
res.pca.expr<-PCA(expr4pca,quali.sup=1)
res.pca.expr

png("results/pca_exp.png", width = 800, height = 600)
plot(res.pca.expr,habillage=1)
dev.off()

## Apply MFA data 

expr.l<-nrow(expr)
cnv.l<-nrow(cnv)

dat4Facto<-data.frame(cond=as.factor(cond),expr.t,cnv.t) 
dim(dat4Facto)

png("results/mfa_layer.png", width = 800, height = 600)
es = MFA(dat4Facto, group=c(1,expr.l,cnv.l), type=c("n",rep("c",2)), ncp=5, 
         name.group=c("cond","RNA","CNV"),num.group.sup=c(1)) 
dev.off()

# Obtain the top features that are highly correlated for the first two principal components and identify the origin 
top10.1 <- sort(es$global.pca$var$cor[,"Dim.1"],decreasing=TRUE)[1:10]
top10.2 <- sort(es$global.pca$var$cor[,"Dim.2"],decreasing=TRUE)[1:10]

# Plot CNV and expression data, together with the top correlated genes in a circos plot
options(stringsAsFactors = FALSE);
names(pollack.nox)
pollack.info <- pollack.nox[, c("chr", "start", "gene")]
expr_circos <- cbind(pollack.info, expr)
head(expr_circos)
cnv_circos <- cbind(pollack.info, cnv)
head(cnv_circos)

# pollack.info contains genomic position
head(pollack.info)

# Top correlated genes with first dimension (all of them come from the expression block)
names(top10.1) <- gsub(".exp","",names(top10.1))
top10.1.df <- as.data.frame(top10.1)
top10.14circos <- merge(x=top10.1.df, y=pollack.info,  by.x =0, by.y="gene")
head(top10.14circos)
top10.14circos <- top10.14circos[,c("chr","start","Row.names","top10.1")]
colnames(top10.14circos) <- c("chr","start","gene","value")

# Top correlated genes with second dimension (all of them come from the CN block)
names(top10.2) <- gsub(".CN","",names(top10.2))
top10.2.df <- as.data.frame(top10.2)
top10.24circos <- merge(x=top10.2.df, y=pollack.info,  by.x =0, by.y="gene")
head(top10.24circos)
top10.24circos <- top10.24circos[,c("chr","start","Row.names","top10.2")]
colnames(top10.14circos) <- c("chr","start","gene","value")

# Plot Circos, is recommended to use Human Genome build 18 (hg18), as it is contemporary to the data so the gene coordinates will match
pdf("results/mfa_circos.pdf")
colors <- rainbow(10, alpha=0.5);
par(mar=c(2, 2, 2, 2));

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
circos(R=370, cir="hg18", W=4,   type="chr", print.chr.lab=TRUE, scale=TRUE);
circos(R=280, cir="hg18", W=100, mapping=cnv_circos,  col.v=4,  type="ml3", 
       cluster=FALSE, col.bar=TRUE, lwd=0.1, col="blue");
circos(R=200, cir="hg18", W=100, mapping=expr_circos,  col.v=4,  type="heatmap2", 
       cluster=FALSE, col.bar=TRUE, lwd=0.1, col="blue");
circos(R=100, cir="hg18", W=50, mapping=top10.14circos,  col.v=4,  type="s", 
       cluster=FALSE, col.bar=TRUE, lwd=0.1, col="darkorchid3", cex = 0.3);
circos(R=100, cir="hg18", W=50, mapping=top10.24circos,  col.v=4,  type="s", 
       cluster=FALSE, col.bar=TRUE, lwd=0.1, col="chartreuse3", cex = 0.3);
circos(R=380, cir="hg18", W=20, mapping=top10.14circos, type="label",
       side="out", col= "darkorchid3", cex=0.25);
circos(R=380, cir="hg18", W=20, mapping=top10.24circos, type="label",
       side="out", col= "chartreuse3", cex=0.25);
dev.off()

### mixOmics DIABLO analysis---- 

# Compute Coefficient of Variation (CV) and extract the genes with the highest variation
cv<-apply(expr,1,sd)/abs(apply(expr,1,mean))
top.100<-head(sort(cv,decreasing=TRUE),100)

cv<-apply(cnv,1,sd)/abs(apply(cnv,1,mean))
top.100<-head(sort(cv,decreasing=TRUE),100)

# mixOmics requires that features in different blocks have different names. .m=mRNA & .v=CNV
expr.f <- expr[names(top.100),]
rownames(expr.f) <- paste(pollack.nox[names(top.100),"gene"],"m",sep=".") #names must be different in the blocks

cnv.f <- cnv[names(top.100),]
rownames(cnv.f) <- paste(pollack.nox[names(top.100),"gene"],"v",sep=".") #names must be different in the blocks

# Format data for mixOmics
colnames(expr.f) <- colnames(cnv.f)
str(expr.f)
str(cnv.f)

# Define Y as a the factor the algorithm needs to predict
cond <- c(rep("cell_line", 4), rep("NORWAY", 28), rep("STANFORD", 9))
Y=as.factor(cond)

# Prepare data X
data = list(rna =t(expr.f), dna = t(cnv.f))

# Design the matrix
design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0
design 

## Apply method
# Set initial model
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, design = design)

# Evaluate performance
set.seed(123) 
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 4, nrepeat = 10) # 10 folds gives stability but might be too much for the number of samples/conditions

# Lists the different outputs
png("results/DIABLO_model_perf.png", width = 800, height = 600)
plot(perf.diablo)
dev.off()

perf.diablo$choice.ncomp$WeightedVote

#Although plots on Mahalanobis distance are not behaving as expected, lets take two components
ncomp = 2

## Feature selection tuning
s.keep <- c(1:10) # Number of features that can be used for each omics block
test.keepX = list (rna = s.keep, dna = s.keep)

# Parallelization params to use all but 1 CPU core
BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores()-1)

# Tunning step
t1 = Sys.time() 
# Computes M-fold or Leave-One-Out Cross-Validation scores based on a user-input grid to determine the optimal parsity parameters 
tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 4, nrepeat = 1,
                              BPPARAM = BPPARAM, 
                              dist = "centroids.dist")
t2 = Sys.time()
running_time = t2 - t1
running_time

# Extract optimal number of RNA and CNV features
list.keepX = tune.TCGA$choice.keepX
list.keepX # The parallelization step doesn't set a fix seed per core, so the result can vary between runs


## Apply method with tuned parameters
sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)

# Selected biomarkers (genes the model kept to build the predictive components)
head(selectVar(sgccda.res)$rna$value,10)
head(selectVar(sgccda.res, comp = 2)$rna$value,10) 
head(selectVar(sgccda.res)$dna$value,10)
head(selectVar(sgccda.res, comp = 2)$dna$value,10) 

## Plots obtained from DIABLO model

# 1. Diagnostic and structural plots
pdf("results/DIABLO_plots.pdf")
plotDiablo(sgccda.res, ncomp = 1)
plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
plotVar(sgccda.res, var.names = TRUE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17), cex = c(2,2), col = c('darkorchid', 'lightgreen'))
circosPlot(sgccda.res, cutoff = 0.67, line = TRUE, 
           color.blocks= c('darkorchid', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5, size.variables=0.5)
plotLoadings(sgccda.res, comp = 1, contrib = 'max', method = 'median')
plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')
dev.off()

# 2. Final multi-omics signature heatmap
png("results/DIABLO_heatmap.png", width = 800, height = 800)
cimDiablo(sgccda.res, margins=c(10,20))
dev.off()

## Test model's predictive power by asking it to predict the labels of the data
predicted <- predict(sgccda.res, newdata = data)
confusion.mat = get.confusion_matrix(truth = Y, 
                                     predicted = predicted$WeightedVote$centroids.dist[,2])
confusion.mat

# Error Calculation
sum(c(confusion.mat[upper.tri(confusion.mat)],confusion.mat[lower.tri(confusion.mat)]))/sum(diag(confusion.mat))

