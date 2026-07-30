exprese <- read.table(
    "featurecounts-ALL.txt",
    sep="\t",
    header=TRUE,
    check.names=FALSE
)

colnames(exprese)[1] <- "Geneid"

ensg <- read.table(
    "Homo_sapiens.GRCh38.95-gene-engs.txt",
    sep="\t",
    header=TRUE,
    check.names=FALSE
)

exprese$gene <- ensg$gene[
    match(exprese$Geneid, ensg$Geneid)
]

write.table(
    exprese,
    "featurecounts-ALL-genes.txt",
    row.names=FALSE,
    sep="\t",
    quote=FALSE
)
