library(scde)
library(ggplot2)
library(biomaRt)
library(GO.db)

GeneName <- read.table("/home/baker/Rna-seq_Data-Analysis/Louisa_Nelson_Single_Cell_Analysis/GeneNames.txt")
head(GeneName)

cdStromal <- read.table("/home/baker/Rna-seq_Data-Analysis/Louisa_Nelson_Single_Cell_Analysis/AllHtseqCountsStromal.tsv", sep=",",header=TRUE, row.names=1)
colnames(cdStromal) <- paste0("Stromal_C",c(1:96))
rownames(cdStromal) <- make.names(GeneName$V1, unique=TRUE)
#colnames(cdStromal) <- sub("_.*", "_Stromal", colnames(cdStromal))
head(cdStromal)


cdTumour <- read.table("/home/baker/Rna-seq_Data-Analysis/Louisa_Nelson_Single_Cell_Analysis/AllHtseqCountsTumour.tsv", sep=",",header=TRUE, row.names=1)
colnames(cdTumour) <- paste0("Tumour_C",c(1:96))
rownames(cdTumour) <- make.names(GeneName$V1, unique=TRUE)
#colnames(cdTumour) <- sub("_.*", "_Tumour", colnames(cdTumour))
head(cdTumour)


cdAll <- cbind(cdStromal, cdTumour)
head(cdAll)


cdCancer <- clean.counts(cdAll, min.lib.size=1000, min.reads = 10, min.detected = 5)
print(dim(cdCancer))

knnCancer <- knn.error.models(cdCancer, k = ncol(cdCancer)/6, n.cores = 2, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)

varinfo <- pagoda.varnorm(knnCancer, counts = cdCancer, trim = 3/ncol(cdCancer), max.adj.var = 5, n.cores = 2, plot = TRUE)


sort(varinfo$arv, decreasing = TRUE)[1:10]

varinfo <- pagoda.subtract.aspect(varinfo, colSums(cdCancer[, rownames(knnCancer)]>0))



# Initialize the connection to the Ensembl BioMart Service
# Available datasets can be listed with 
# listDatasets(useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))
# Use mmusculus_gene_ensembl for mouse
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")

# Constructs a dataframe with two columns: hgnc_symbol and go_id
# If rownames are Ensembl IDs, use ensembl_gene_id as filter value
go <- getBM(attributes = c("hgnc_symbol", "go_id"), filters = "hgnc_symbol", values = rownames(cdCancer), mart = ensembl)

# Use the GO.db library to add a column with the GO-term to the dataframe
go$term <- Term(go$go_id)

# Create a named list of character vectors out of the df
s = split(go$hgnc_symbol, paste(go$go_id,go$term))

# Saves the list as a R environment
go.env <- list2env(s)

# Test
class(go.env)


save.image("Heterogeneity_Analysis.RData")

GOAnnotated <- cbind(go$ensembl_gene_id,go$go_id,go$term)

head(ls(go.env)) # Look at gene set names
head(get(ls(go.env)[1], go.env)) # Look at one gene set
ngenes <- unlist(lapply(as.list(go.env),function(x) sum(x %in% rownames(varinfo$mat))))
head(ngenes,20)

pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 2, verbose=TRUE)
save.image("Heterogeneity_Analysis.RData")


df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)


head(df)

save.image("Heterogeneity_Analysis.RData")



clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = 2, plot = TRUE)

# get full info on the top aspects
tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
# determine overall cell clustering
hc <- pagoda.cluster.cells(tam, varinfo)
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)
tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)

save.image("Heterogeneity_Analysis.RData")


library(Rtsne);
col.cols <- rbind(groups = cutree(hc, 3))
# recalculate clustering distance .. we'll need to specify return.details=T
cell.clustering <- pagoda.cluster.cells(tam,varinfo,include.aspects=TRUE,verbose=TRUE,return.details=T)
# fix the seed to ensure reproducible results
set.seed(0); 
tSNE.pagoda <- Rtsne(cell.clustering$distance,is_distance=T,initial_dims=100,perplexity=10)
par(mfrow=c(1,1), mar = c(2.5,2.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
plot(tSNE.pagoda$Y,col=adjustcolor(col.cols,alpha=0.5),cex=1,pch=19,xlab="",ylab="")
rownames(tSNE.pagoda$Y)=colnames(cdCancer)

save.image("Heterogeneity_Analysis.RData")







