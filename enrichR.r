##### Perform enrichment analysis using enrichR package in R #####



#install.packages("enrichR")
library(enrichR)
setEnrichrSite("Enrichr")

dbs <- listEnrichrDbs()
dbs <- c("KEGG_2019_Human","GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021","KEGG_2021_Human")

### Knockout database ###
x = read.csv(filename,sep="\t",header=T,dec=",")   #### Change filename to significant gene file name###


## Note: Change your gene column name to external_gene_name and p-value column to p.value
enriched <- enrichr(subset(x,p.value <= 0.01)$external_gene_name, dbs) ## Subset a p < 0.01


tiff(filename,width=8,height=5,units="in",res=300)
plotEnrich(enriched$KEGG_2019_Human, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",title = "Top 20 KEGG Pathways (p < 0.01)")
dev.off()

### You can write other dbs also
write.table(enriched$KEGG_2021_Human,"KeggTable.txt",sep=",",row.names=F,quote=F,dec=",")

