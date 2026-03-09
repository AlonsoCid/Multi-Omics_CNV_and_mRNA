library(mixOmics)
library(OmicCircos)
library(BiocParallel)

# Load data
load("data/pollack.RData")

# Extract omic layers
cnv <- pollack.nox[, 8:48]
expr <- pollack.nox[, 49:89]

## 1. Feature Selection (Coefficient of Variation)

# Compute Coefficient of Variation (CV) and extract the genes with the highest variation
cv <- apply(expr,1,sd)/abs(apply(expr,1,mean))
top.100 <- head(sort(cv,decreasing=TRUE),100)
expr.f <- expr[names(top.100),]
rownames(expr.f) <- paste(pollack.nox[names(top.100),"gene"],"m",sep=".") #names must be different in the blocks .m=mRNA

cv <- apply(cnv,1,sd)/abs(apply(cnv,1,mean))
top.100 <- head(sort(cv,decreasing=TRUE),100)
cnv.f <- cnv[names(top.100),]
rownames(cnv.f) <- paste(pollack.nox[names(top.100),"gene"],"v",sep=".") #names must be different in the blocks .v=CNV

# Format data for mixOmics
colnames(expr.f) <- colnames(cnv.f)

# Define Y (Outcome factor)
cond <- c(rep("cell_line", 4), rep("NORWAY", 28), rep("STANFORD", 9))
Y <- as.factor(cond)

# Prepare data X (Blocks transposed so samples are in rows)
data <- list(rna = t(expr.f), dna = t(cnv.f))

# Null design matrix (prioritizes outcome discrimination over block correlation)
design <- matrix(0.1, ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
diag(design) <- 0

## 2. Initial Model & Performance Evaluation

# Set initial model
sgccda.res <- block.splsda(X = data, Y = Y, ncomp = 5, design = design)

# Evaluate performance (10 folds gives stability)
set.seed(123) 
perf.diablo <- perf(sgccda.res, validation = 'Mfold', folds = 4, nrepeat = 10) 

# Save model performance plot
png("results/DIABLO_model_perf.png", width = 800, height = 600)
plot(perf.diablo)
dev.off()

# Take two components based on performance plateau
ncomp <- 2

## 3. Feature Selection Tuning

s.keep <- c(1:10) # Grid of features to test
test.keepX <- list(rna = s.keep, dna = s.keep)

# Parallelization params to use all but 1 CPU core
BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores()-1)

# Tuning step (Computes Cross-Validation scores based on grid)
tune.TCGA <- tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
                               test.keepX = test.keepX, design = design,
                               validation = 'Mfold', folds = 4, nrepeat = 1,
                               BPPARAM = BPPARAM, 
                               dist = "centroids.dist")

# Extract optimal number of RNA and CNV features
list.keepX <- tune.TCGA$choice.keepX

## 4. Final Model Application & Plotting

# Apply method with tuned parameters
sgccda.res <- block.splsda(X = data, Y = Y, ncomp = ncomp, 
                           keepX = list.keepX, design = design)

# Generate Diagnostic and structural plots (Saved as a single multipage PDF)
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

# Generate final multi-omics signature heatmap 
png("results/DIABLO_heatmap.png", width = 800, height = 800)
cimDiablo(sgccda.res, margins=c(10,20))
dev.off()

## 5. Model Predictive Power Test

predicted <- predict(sgccda.res, newdata = data)
confusion.mat <- get.confusion_matrix(truth = Y, 
                                      predicted = predicted$WeightedVote$centroids.dist[,2])
print("Confusion Matrix:")
print(confusion.mat)

# Error Calculation
error_rate <- sum(c(confusion.mat[upper.tri(confusion.mat)], confusion.mat[lower.tri(confusion.mat)])) / sum(diag(confusion.mat))
print(paste("Training Error Rate:", error_rate))
