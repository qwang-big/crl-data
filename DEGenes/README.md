# Append DEGs to the ECDF results
```r
library(DESeq2)
countdata <- read.table("CLLcounts.txt", header=TRUE, row.names=1)
countname <- rownames(countdata)
countdata["CD74",]=countdata["CD74",]/100
countdata <- as.matrix(apply(countdata[,c(grep("cancer",colnames(countdata)),grep("normal",colnames(countdata)))],2,function(d) as.integer(round(d))))
rownames(countdata) <- countname
condition <- factor(c(rep("exp", 10), rep("ctl", 5)))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, fitType='local')
res <- results(dds)
res <- res[order(res$padj), ]
res <- res[rownames(res) %in% cll$pg$prom, ]
df  <- rbind(df,data.frame(value=sort(match(markers$CLL, rownames(res))),type="DEG",case="CLL"))

countdata <- read.table("CRCcounts.txt", header=TRUE, row.names=1)
countname <- rownames(countdata)
countdata["B2M",]=countdata["B2M",]/100
countdata <- as.matrix(apply(countdata[,c(grep("cancer",colnames(countdata)),grep("normal",colnames(countdata)))],2,function(d) as.integer(round(d))))
rownames(countdata) <- countname
condition <- factor(c(rep("exp", 6), rep("ctl", 15)))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, fitType='local')
res <- results(dds)
res <- res[order(res$padj), ]
res <- res[rownames(res) %in% crc$pg$prom, ]
df  <- rbind(df,data.frame(value=sort(match(markers$CRC, rownames(res))),type="DEG",case="CRC"))

countdata <- read.table("PTCcounts.txt", header=TRUE, row.names=1)
countname <- rownames(countdata)
countdata["FN1",]=countdata["FN1",]/100
countdata["CLU",]=countdata["CLU",]/100
countdata["TG",]=countdata["TG",]/100
countdata <- as.matrix(apply(countdata[,c(grep("cancer",colnames(countdata)),grep("normal",colnames(countdata)))],2,function(d) as.integer(round(d))))
rownames(countdata) <- countname
condition <- factor(c(rep("exp", 3), rep("ctl", 7)))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, fitType='local')
res <- results(dds)
res <- res[order(res$padj), ]
res <- res[rownames(res) %in% ptc$pg$prom, ]
df  <- rbind(df,data.frame(value=sort(match(markers$PTC, rownames(res))),type="DEG",case="PTC"))

countdata <- read.table("mCLLcounts.txt", header=TRUE, row.names=1)
countname <- rownames(countdata)
countdata <- as.matrix(apply(countdata[,c(grep("CLL",colnames(countdata)),grep("NC",colnames(countdata)))],2,function(d) as.integer(round(d))))
rownames(countdata) <- countname
condition <- factor(c(rep("exp", 2), rep("ctl", 3)))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, fitType='local')
res <- results(dds)
res <- res[order(res$padj), ]
res <- res[rownames(res) %in% mcll$pg$prom, ]
df  <- rbind(df,data.frame(value=sort(match(markers$mCLL, rownames(res))),type="DEG",case="mCLL"))

countdata <- read.table("LGGcounts.txt", header=TRUE, row.names=1)
countname <- rownames(countdata)
countdata <- as.matrix(apply(countdata[,c(grep("cancer",colnames(countdata)),grep("normal",colnames(countdata)))],2,function(d) as.integer(round(d))))
rownames(countdata) <- countname
condition <- factor(c(rep("exp", 3), rep("ctl", 2)))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, fitType='local')
res <- results(dds)
res <- res[order(res$padj), ]
res <- res[rownames(res) %in% glm$pg$prom, ]
df  <- rbind(df,data.frame(value=sort(match(markers$LGG, rownames(res))),type="DEG",case="LGG"))

countdata <- read.table("EScounts.txt", header=TRUE, row.names=1)
countname <- rownames(countdata)
countdata <- as.matrix(apply(countdata[,c(1,2,7,8)],2,function(d) as.integer(round(d))))
rownames(countdata) <- countname
condition <- factor(c(rep("exp", 2), rep("ctl", 2)))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, fitType='local')
res <- results(dds)
res <- res[order(res$padj), ]
res <- res[rownames(res) %in% tsc$pg$prom, ]
df  <- rbind(df,data.frame(value=sort(match(markers$TSC, rownames(res))),type="DEG",case="TSC"))

countdata <- read.table("EScounts.txt", header=TRUE, row.names=1)
countname <- rownames(countdata)
countdata <- as.matrix(apply(countdata[,c(3,4,7,8)],2,function(d) as.integer(round(d))))
rownames(countdata) <- countname
condition <- factor(c(rep("exp", 2), rep("ctl", 2)))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, fitType='local')
res <- results(dds)
res <- res[order(res$padj), ]
res <- res[rownames(res) %in% msc$pg$prom, ]
df  <- rbind(df,data.frame(value=sort(match(markers$MSC, rownames(res))),type="DEG",case="MSC"))

countdata <- read.table("EScounts.txt", header=TRUE, row.names=1)
countname <- rownames(countdata)
countdata <- as.matrix(apply(countdata[,c(5,6,7,8)],2,function(d) as.integer(round(d))))
rownames(countdata) <- countname
condition <- factor(c(rep("exp", 2), rep("ctl", 2)))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, fitType='local')
res <- results(dds)
res <- res[order(res$padj), ]
res <- res[rownames(res) %in% npc$pg$prom, ]
df  <- rbind(df,data.frame(value=sort(match(markers$NPC, rownames(res))),type="DEG",case="NPC"))
```