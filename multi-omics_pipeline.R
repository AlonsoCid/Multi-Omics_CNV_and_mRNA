library(FactoMineR)
library(OmicCircos)
library(mixOmics)

### FactoMinerR analysis----

load("pollack.RData")
head(pollack.nox)
names(pollack.nox)

## Copy Number Variant PCS
cnv <- pollack.nox[, 8:48]
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
expr <- pollack.nox[, 49:89]
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
