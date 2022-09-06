## Code for converting read counts into expression values using TMM normalization
## Files required are read count results and annotation file of the sample annotation
# file must have a column Lab_Code same as colnames of read count table to make groups
# It uses biomart for annotation of the ensembl_gene_id
# Ensembl gene ids supported only, but can be updated in code manually

library(edgeR)
normalized_using_TMM <- function(readCount,annotation_file,atria=T,Lab_Code="Lab_Code",
                                 output_file="Expression data normalized annotated.txt"){
## Reading and correcting count table and colnames ##
x = read.csv(file = readCount,sep="\t",header=T,dec=".")

## If atria is used for trimming
if(atria == TRUE)
  colnames(x) = gsub(".atria_nat","",colnames(x))

## Reading Annotation File ##
anno_x = read.csv(file = annotation_file,sep="\t",header=T)
## Lab code required and should be same as colnames of read counts data
anno_x = anno_x[order(anno_x[,Lab_Code]),]

## Creating group based on sampleName ##
group <- anno_x[,Lab_Code]   ## Changed to individual sample

## Create DGEList with groups ##
y <- DGEList(counts=x, group=group)

## Filteration ##
keep <- filterByExpr(y, group=group,min.count=5,min.total.count=10)
y2 = y[which(keep),]

## TMM Normalization ##
y2 = calcNormFactors(y2,method="TMM")
exprs_data = cpm(y2,log=T)

## Annotation of variables ##
library(biomaRt)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

transcript_ids <- rownames(exprs_data)

res <- getBM(attributes = c('ensembl_gene_id','chromosome_name',"start_position",
                            'external_gene_name'),
             filters = 'ensembl_gene_id', 
             values = transcript_ids,
             mart = mart)
exprs_res = data.frame(exprs_data)
exprs_res = cbind(ensembl_gene_id=rownames(exprs_res),exprs_res)
exprs_res <- merge(res,exprs_res,by="ensembl_gene_id",all.y=T)
write.table(exprs_res,file = output_file,sep="\t",row.names=F,quote=F,dec=",",na = "")
}