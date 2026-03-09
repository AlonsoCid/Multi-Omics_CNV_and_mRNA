library(FactoMineR)
library(OmicCircos)

# Load data
load("data/pollack.RData")

# Extract omic layers
cnv <- pollack.nox[, 8:48]
expr <- pollack.nox[, 49:89]

## 1. Data Preparation for MFA

# CNV preparation
cnv.t <- t(cnv)
colnames(cnv.t) <- paste(pollack.nox$gene, "CN", sep=".")

# Expression preparation
colnames(expr) <- colnames(cnv) # Names must match
expr.t <- t(expr)
colnames(expr.t) <- paste(pollack.nox$gene, "exp", sep=".")

# Classify samples by origin
cond <- c(rep("cell_line", 4), rep("NORWAY", 28), rep("STANFORD", 9))

# Construct unified data.frame
expr.l <- nrow(expr)
cnv.l <- nrow(cnv)
dat4Facto <- data.frame(cond=as.factor(cond), expr.t, cnv.t) 

## 2. Apply Multiple Factor Analysis (MFA)

# Calculate MFA (graph=FALSE prevents cluster crashing from X11 popups)
es = MFA(dat4Facto, group=c(1, expr.l, cnv.l), type=c("n", rep("c", 2)), ncp=5, 
         name.group=c("cond", "RNA", "CNV"), num.group.sup=c(1), graph=FALSE) 

# Plot the group representation (omics layer contributions)
png("results/mfa_layer.png", width = 800, height = 600)
plot(es, choix="group", title="MFA: Groups Representation")
dev.off()

## 3. Extract Top Correlated Features

# Obtain top 10 features highly correlated for the first two dimensions
top10.1 <- sort(es$global.pca$var$cor[,"Dim.1"], decreasing=TRUE)[1:10]
top10.2 <- sort(es$global.pca$var$cor[,"Dim.2"], decreasing=TRUE)[1:10]

## 4. Circos Plot Data Formatting

options(stringsAsFactors = FALSE)
pollack.info <- pollack.nox[, c("chr", "start", "gene")]

expr_circos <- cbind(pollack.info, expr)
cnv_circos <- cbind(pollack.info, cnv)

# Format Dimension 1 genes (RNA block)
names(top10.1) <- gsub(".exp", "", names(top10.1))
top10.1.df <- as.data.frame(top10.1)
top10.14circos <- merge(x=top10.1.df, y=pollack.info, by.x=0, by.y="gene")
top10.14circos <- top10.14circos[, c("chr", "start", "Row.names", "top10.1")]
colnames(top10.14circos) <- c("chr", "start", "gene", "value")

# Format Dimension 2 genes (CNV block)
names(top10.2) <- gsub(".CN", "", names(top10.2))
top10.2.df <- as.data.frame(top10.2)
top10.24circos <- merge(x=top10.2.df, y=pollack.info, by.x=0, by.y="gene")
top10.24circos <- top10.24circos[, c("chr", "start", "Row.names", "top10.2")]
colnames(top10.24circos) <- c("chr", "start", "gene", "value")

## 5. Generate Circos Plot

# Using hg18 as it matches the contemporary coordinates of the dataset
pdf("results/mfa_circos.pdf")
colors <- rainbow(10, alpha=0.5)
par(mar=c(2, 2, 2, 2))

plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")
circos(R=370, cir="hg18", W=4, type="chr", print.chr.lab=TRUE, scale=TRUE)
circos(R=280, cir="hg18", W=100, mapping=cnv_circos, col.v=4, type="ml3", 
       cluster=FALSE, col.bar=TRUE, lwd=0.1, col="blue")
circos(R=200, cir="hg18", W=100, mapping=expr_circos, col.v=4, type="heatmap2", 
       cluster=FALSE, col.bar=TRUE, lwd=0.1, col="blue")
circos(R=100, cir="hg18", W=50, mapping=top10.14circos, col.v=4, type="s", 
       cluster=FALSE, col.bar=TRUE, lwd=0.1, col="darkorchid3", cex = 0.3)
circos(R=100, cir="hg18", W=50, mapping=top10.24circos, col.v=4, type="s", 
       cluster=FALSE, col.bar=TRUE, lwd=0.1, col="chartreuse3", cex = 0.3)
circos(R=380, cir="hg18", W=20, mapping=top10.14circos, type="label",
       side="out", col="darkorchid3", cex=0.25)
circos(R=380, cir="hg18", W=20, mapping=top10.24circos, type="label",
       side="out", col="chartreuse3", cex=0.25)
dev.off()