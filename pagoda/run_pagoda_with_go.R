library("stats")
library("methods")
library("utils")
library("scde") #installed
library("RColorBrewer") #installed
library("gplots") #installed

printf <- function(...) cat(sprintf(...))

min.lib.size <- 0
min.reads <- 10
min.detected <- 0
nAspects <- 7

ncores <- 20

args <- commandArgs(trailingOnly = TRUE)
pid <- args[1]
sample <- args[2]
total_reads <- as.numeric(args[3])
coding_genes <- as.numeric(args[4])
mito_frac <- as.numeric(args[5])
non_genic_frac <- as.numeric(args[6])

printf("PID: %s\n", pid)
printf("Sample: %s\n", sample)

output.dir <- sprintf("pagoda_plots/%s_%s_%.f_with_go", pid, sample, as.numeric(as.POSIXct(Sys.time())))
dir.create(output.dir)
setwd(output.dir)

print("Loading expression matrix...")
load(sprintf("Rdata/expression-%s-%s-%d_%d_%d_%d.Rdata", pid, sample, total_reads, coding_genes, mito_frac, non_genic_frac))
count.matrix <- expression.matrices$count
rownames(count.matrix) <- gene.names
ccounts <- clean.counts(count.matrix, min.lib.size, min.reads, min.detected)

print("Running KNN error modeling...")
# ----KnearestNeighbor Error Models, eval = FALSE-------------------------
ko.ifm <- knn.error.models( counts = ccounts,
                            k = 30,
                            n.cores = ncores,
                            min.count.threshold = 2,
                            min.nonfailed = 5,
                            max.model.plots=20)

print("Normalizing variance...")
# ----Normalize Variance, eval = FALSE------------------------------------
png(file=file.path(output.dir,"Plot_Varinfo.png"),width = 1000, height = 1000)
varinfo <- pagoda.varnorm(  ko.ifm,
                            counts = ccounts,
                            trim = 3/ncol(ccounts),
                            max.adj.var = 5,
                            n.cores = ncores,
                            plot = TRUE,
                            weight.df.power = 2)
dev.off()

print("Selecting gene-sets...")
parse_goterms <- function (f) {
    unlist(lapply(f, function (x) { genes <- unlist(strsplit(x, "\t")); l <- list(genes[-1]); names(l) <- c(genes[1]); l } ), recursive=FALSE)
}

f <- readLines("goterms/hallmark.txt")
go.env <- parse_goterms(f)
f <- readLines("goterms/wholego.txt")
go.env <- c(go.env, parse_goterms(f))

go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env)

print("Controlling sequencing depth...")
# Controlling for sequencing depth
varinfo <- pagoda.subtract.aspect(varinfo,
                                  colSums(ccounts[,rownames(ko.ifm)]>0))

print("Running PwPCA...")
pwpca.nC <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 2, n.cores = ncores)


# Create an Environment with the Gene sets and included genes

#
png(file=file.path(output.dir,"NovelClusters_clpca_Plot.png"),width = 1000, height = 1000)
clpca.nC <- pagoda.gene.clusters( varinfo, trim = 7.1/ncol(varinfo$mat),
                                  n.clusters = 250, n.cores = ncores, plot = TRUE,
                                  n.components = 2, n.samples=100, n.starts=15)
dev.off()


png(file=file.path(output.dir,"NovelClusters_VarianceToSize.png"),width = 1000, height = 1000)
df.nC <- pagoda.top.aspects(pwpca.nC,clpca.nC, return.table = TRUE, plot = TRUE, z.score = 1.96)
dev.off()


png(file=file.path(output.dir,"NovelClusters_Pagoda_TOP_Aspects.png"),width = 1000, height = 1000)
tam.nC <- pagoda.top.aspects(pwpca.nC,clpca.nC, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
dev.off()

print("Clustering cells...")
hc <- pagoda.cluster.cells(tam.nC, varinfo,include.aspects = TRUE)

print("Lodging/reducing redundancy...")
png(file=file.path(output.dir,"NovelClusters_RLoadRedundancy.png"),width = 1000, height = 1000)
tamr.nC <- pagoda.reduce.loading.redundancy(tam.nC, clpca.nC,plot=TRUE,corr.power=7)
dev.off()

colors <- rep("#1B9E77", ncol(ccounts))
col.cols <- rbind(class  = colors)
colnames(col.cols) <- colnames(ccounts)

png(file=file.path(output.dir,"NovelClusters_Pagoda_ReduceRedundancy_Aspects.png"),width = 1000, height = 1000)
tamr2.nC <- pagoda.reduce.redundancy(tamr.nC, top=nAspects, distance.threshold = 0.7, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 3, col.cols=col.cols)
dev.off()

print("Running view aspects...")
png(file=file.path(output.dir,"NovelClusters_Pagoda_View_Aspects.png"),width = 1000, height = 1000)
pagoda.view.aspects(tamr2.nC, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols=col.cols)
dev.off()

library(Rtsne);
# recalculate clustering distance .. we'll need to specify return.details=T
print("Reclustering cells and running t-SNE...")
cell.clustering <- pagoda.cluster.cells(tam.nC,varinfo,include.aspects=TRUE,verbose=TRUE,return.details=T)
#perplexity <- ifelse(nrow(cell.clustering$distance) - 1 < 21, 3, 7)
perplexity <- 3
tSNE.pagoda <- Rtsne(cell.clustering$distance,is_distance=T,initial_dims=100,perplexity=perplexity,max_iter=5000)

png(file=file.path(output.dir,"NovelClusters_Pagoda_tSNE.png"),width = 1000, height = 1000)
par(mfrow=c(1,1), mar = c(2.5,2.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
plot(tSNE.pagoda$Y,col=adjustcolor(col.cols,alpha=0.5),cex=1,pch=19,xlab="",ylab="")
dev.off()

rownames(tSNE.pagoda$Y) <- colnames(tam.nC$xv)

print("Writing to files...")
write.table(unlist(lapply(clpca.nC$cluster, function(x){paste(x, collapse=",")})), "gene_clusters.csv", sep=",", col.names=F)
app <- make.pagoda.app(tamr2.nC, tam.nC, varinfo, go.env, pwpca.nC, clpca.nC, col.cols = col.cols, cell.clustering = hc, title = paste0(pid, sample), embedding = tSNE.pagoda$Y)
save(app, file=file.path(output.dir,paste("NovelClusters_Genesets_app.RData",sep="")))
save.image(file=file.path(output.dir,"NovelClusters.RData"))
# load(file=file.path(output.dir,"NovelClusters.RData"))
