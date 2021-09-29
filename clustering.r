library("rstudioapi")
library(dplyr)
library(edgeR)
library(org.Hs.eg.db)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(factoextra)

wd <- rstudioapi::getSourceEditorContext()$path
wd1 <- strsplit(wd, "/")
wd1 <- paste0(wd1[[1]][1:lengths(wd1)-1], collapse = "/")
knitr::opts_knit$set(root.dir = setwd(wd1))


filename1 <- paste0(wd1, '/' , 'pcoa.png')
filename2 <- paste0(wd1, '/' , 'pca.png')
filename3 <- paste0(wd1, '/' , 'pca.scree.png')

# BEGIN Block of code is from the edgeR user guide case study 4.1
rawdata <- read.delim("supporting.files/TableS1.txt", check.names=FALSE, stringsAsFactors=FALSE)

y <- DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])

idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
y <- y[idfound,]
dim(y)

egREFSEQ <- toTable(org.Hs.egREFSEQ)

head(egREFSEQ)

m <- match(y$genes$RefSeqID, egREFSEQ$accession)
y$genes$EntrezGene <- egREFSEQ$gene_id[m]

egSYMBOL <- toTable(org.Hs.egSYMBOL)
head(egSYMBOL)

m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol <- egSYMBOL$symbol[m]
head(y$genes)

o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Symbol)
y <- y[!d,]
nrow(y)

y$samples$lib.size <- colSums(y$counts)

rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL

y <- calcNormFactors(y)
y$samples

plotMDS(y)
# END Block of code is from the edgeR user guide case study 4.1

# Get data for the PCoA and the PCA
pcoa.data <- plotMDS(y)
pcoa.data <- pcoa.data[[3]]

pca.counts <- cpm(y, log = TRUE)
pca.counts <- pca.counts %>%
  t()
pca <- prcomp(pca.counts)
pca.data <- pca$x

# formatting the PCoA data
colnames(pcoa.data) <- c("x", "y")
pcoa.data <- as.data.frame(pcoa.data)
# convert rownames to a column
pcoa.data <- tibble::rownames_to_column(pcoa.data, "patient.cell")

pcoa.data <- pcoa.data %>%
  separate(col = "patient.cell", into = c('patient', 'cell'), sep = -1) %>%
  mutate(cell.de = case_when(cell =="N" ~ "Normal", 
                               cell == "T" ~ "Cancer"))

#Formatting the PCA data
pca.data <- as.data.frame(pca.data)
pca.data <- tibble::rownames_to_column(pca.data, "patient.cell")
pca.data <- pca.data %>%
  separate(col = "patient.cell", into = c('patient', 'cell'), sep = -1) %>%
  mutate(cell.de = case_when(cell =="N" ~ "Normal", 
                             cell == "T" ~ "Cancer"))

# We need to obtain the percentages for each of the principal components
percentage <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
#give column names to the percentage
percentage <- paste(colnames(pca$x), "(", 
                    paste(as.character(percentage), "%", ")", sep="") )


png(filename1, height = 1200, width = 1200)
p <- ggplot(pcoa.data1, aes(x = x, y = y, colour = cell.de)) +
  geom_point(aes(shape=cell.de, color=cell.de), size = 12) + 
  geom_text_repel(aes(label = patient), size=10, 
                  max.overlaps = Inf, show.legend = FALSE, 
                  colour = "black", point.padding = 20) + 
  xlab("Leading logFC dim1") + 
  ylab("Leading logFC dim2") + 
  theme(plot.title = element_text(size = rel(6), hjust = 0.5), 
        axis.title = element_text(size = rel(5)), 
        axis.text = element_text(size = rel(3)), 
        legend.key.size = unit(1, 'cm'), legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_blank(),
        legend.text = element_text(size=25))
print(p)
dev.off()


png(filename2, height = 1200, width = 1200)
p <- ggplot(pca.data1, aes(x = PC1, y = PC2, colour = cell.de)) +
  geom_point(aes(shape=cell.de, color=cell.de), size = 12) + 
  geom_text_repel(aes(label = patient), size=10, max.overlaps = Inf, 
                  show.legend = FALSE, colour = "black", 
                  point.padding = 20) + 
  xlab(percentage[1]) + 
  ylab(percentage[2]) + 
  theme(plot.title = element_text(size = rel(6), hjust = 0.5), 
        axis.title = element_text(size = rel(5)), 
        axis.text = element_text(size = rel(3)), 
        legend.key.size = unit(1, 'cm'), legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_blank(), 
        legend.text = element_text(size=25)) 
print(p)
dev.off()


png(filename3, height = 1200, width = 1200)
p <- fviz_eig(pca) + 
  theme(text=element_text(size=35), axis.title = element_text(size=42), 
        axis.text = element_text(size=35))
print(p)
dev.off()
