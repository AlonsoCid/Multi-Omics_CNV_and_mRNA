library(FactoMineR)

# Load data
load("data/pollack.RData")

# Extract omic layers
cnv <- pollack.nox[, 8:48]
expr <- pollack.nox[, 49:89]

## 1. Copy Number Variant PCA

cnv.t <- t(cnv)
# Assign names, we include a exp suffix to differentiate genes from expression
colnames(cnv.t) <- paste(pollack.nox$gene, "CN", sep=".")

# Classify the samples by origin
cond <- c(rep("cell_line", 4), rep("NORWAY", 28), rep("STANFORD", 9))

# Construct data.frame to perform PCA
cnv4pca <- data.frame(cond, cnv.t)
res.pca.cnv <- PCA(cnv4pca, quali.sup=1, graph=FALSE) # Suppress auto-plotting

# Save plot
png("results/pca.png", width = 800, height = 600)
plot(res.pca.cnv, habillage=1, title="PCA: Copy Number Variants")
dev.off()

## 2. Expression PCA

colnames(expr) <- colnames(cnv) # To perform later MFA, we need to have the same names
expr.t <- t(expr)

# Assign names and include a exp suffix to differentiate genes from cnv
colnames(expr.t) <- paste(pollack.nox$gene, "exp", sep=".")

# Construct data.frame to perform PCA
expr4pca <- data.frame(cond, expr.t)
res.pca.expr <- PCA(expr4pca, quali.sup=1, graph=FALSE) # Suppress auto-plotting

# Save plot
png("results/pca_exp.png", width = 800, height = 600)
plot(res.pca.expr, habillage=1, title="PCA: Gene Expression")
dev.off()
