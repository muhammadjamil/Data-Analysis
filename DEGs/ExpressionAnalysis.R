rm(list=ls())
asCN <- function(x){as.numeric(as.character(x))}
## reading tables ##
exprsData = read.csv("TMM Normalized Expression.txt",sep="\t",header=T,dec=",")

### Creating expression table and variable annotation table ####
exprs <- exprsData[,c(4:13)]  ## Check for the column
anno <- exprsData[,c(1:3)] ## check for the column
design <- c(rep(1,5),rep(2,5)) ## Recheck the design based on column sorting

    dat2 = cbind(V1=design,t(exprs))
    model = sapply(2:ncol(dat2),function(i) {
      tmp <- data.frame(V1=dat2[,1],V2=dat2[,i])
      summary(lm(V2~V1,tmp))$coefficients[2,4]
    })
    
    ### Creating dataframe Res as final output
    res = data.frame(anno,
                     logFC = colMeans(dat2[which(dat2[,1] == 1),2:ncol(dat2)])-colMeans(dat2[which(dat2[,1] == 2),2:ncol(dat2)]),
                     P.Value = model,
                     q.value = p.adjust(model,"BH"))
   