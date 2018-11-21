library(ProCoNA)


subsetModCors <- function (modCors, modules) {
### subsets the module-phenotype correlation matrix which has funny rownames
##<< The matrix of module-phenotype correlations
##<< Which modules are desired.
  rows <- rownames(modCors)  # like ME3
  ms <- sapply(modules, function(x) paste("ME", x, sep=""))
  modCors[which(rows %in% ms),]
}

# definde a function to output particular modules info
moduleData <- function(pepnet, pepcors, module, pepinfo, fileprefix) {
  moduleX <- peptides(pepnet)[which(mergedColors(pepnet)==module)]
  moduleInfo <- pepinfo[which(pepinfo$Mass_Tag_ID %in% moduleX),]
  moduleCors <- pepcors[which(pepcors$Module==module),]
  corname <- paste(fileprefix, "_correlations.csv", sep="")
  write.table(moduleCors, file=corname, sep=",", row.names=F)
  infoname <- paste(fileprefix, "_peptide_info.csv", sep="")
  write.table(moduleInfo, file=infoname, sep=",", row.names=F)
}

data(ProCoNA_Data)  # import the example data (rows = samples; column = peptides(proteins))
peptideData <- subsetPeptideData(peptideData, percentageNAsAllowed=0.2)

# setting parameters used in this script #
ALLOW_WGCNA_THREADS <- 1 # set the multi-threads of WGCNA
textSize <- 0.7        # set text size in heat map
title <- "Heatmap of correlation between modules & phenotypes" # heatmap title
plotName <- "heatmap.pdf"

#dim(peptideData)   # the dimension peptideData
#----- choose the softthreshold beta ----------------------#
# set the candidate soft threshold
powers <- c(c(1:10), seq(from=12, to=20, by=2))
# use WGCNA packages to pick soft threshold beta
sft <- pickSoftThreshold(peptideData, networkType="signed", RsquaredCut=0.8, powerVector=powers)

# set the graph
sizeGrWindow(9, 5)
par(mfrow=c(1,2))
cex1 <- 0.9

# plotting R^2 graph and mean connectivity to find soft threshold
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",
	   ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")
    )
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
	   main=paste("Mean connectivity")
	 )
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# beta$powerEstimate # this recommands the proper soft threshold beta

#--- Build the network ------- #
peptideNetwork <- buildProconaNetwork(networkName="my network",
                                      pepdat=peptideData,
                                      networkType="signed",
                                      pow=sft$powerEstmate,
                                      pearson=FALSE,
                                      toPermTestPermutes=1000)
net <- peptideNetwork
# plot the hclust graph
pdf("treegraph.pdf")
plotNet(peptideNetwork) 
dev.off()

if (nrow(phenotypes) != nrow(mergedMEs(net))) {
  stop("Net must have merged eigenvectors with dimensions matching the phenotype data! See mergedMEs(net).")
}

modules <- unique(mergedColors(net))
modCors <- modulePhenotypeCorrelations(net, phenotypes)
modCors <- subsetModCors(modCors, modules)
mms <- 1:(ncol(modCors)/2)
mmps <- ((ncol(modCors)/2)+1):ncol(modCors)
modCorMM <- as.matrix(modCors[,mms])
modCorPS <- as.matrix(modCors[,mmps])
# pastes into a vector down columns
textMatrix <- paste(signif(modCorMM, 2), "\n(", sep="")
textMatrix <- paste(textMatrix, paste(signif(modCorPS, 1), ")", sep = ""), sep="")
textMatrix <- matrix(textMatrix, byrow=F, nrow=nrow(modCorMM))  
# Display the correlation values within a heatmap plot
par(mar = c(6, 8.5, 3, 3))
if(!is.null(plotName)) {
  pdf(plotName)
}
colnames(phenotypes) <- c("1","2","3","4","5","6","7","8","9","10")
labeledHeatmap(Matrix = modCors[,mms], 
               xLabels = colnames(phenotypes), 
               yLabels = rownames(modCors),
               ySymbols = rownames(modCors),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = textSize,
               zlim = c(-1,1),
               main = title)
if(!is.null(plotName)) {
  dev.off()
}

pepcor <- moduleMemberCorrelations(pnet=peptideNetwork,
                                   pepdat=peptideData,
                                   phenotypes=phenotypes)

for (i in 0:max(pepcor$Module)) {
  moduleData(peptideNetwork, pepcor, i, masstagdb, paste("Module_",as.numeric(i)))
}

#Prepare csv file for gephi visualization 
TOMmatrix <- net@TOM
pairs <- (dim(TOMmatrix)[1]*(dim(TOMmatrix)[1]-1)/2)
nodename1 <- vector(mode = "character",length = pairs)
nodename2 <- vector(mode = "character",length = pairs)
TOMdis <- vector(mode = "numeric",length = pairs)
count = 1
for (j in 1:dim(TOMmatrix)[1]) {
  if (j < dim(TOMmatrix)[1]) {
    for (eachnode in (j+1):dim(TOMmatrix)[1]){
      nodename1[count] <- rownames(TOMmatrix)[j]
      nodename2[count] <- rownames(TOMmatrix)[eachnode]
      TOMdis[count] <- TOMmatrix[j,eachnode]
      count=count + 1
    }
  }
}

network <- data.frame(Source=nodename1, Target=nodename2, Weight=TOMdis, 
                      Type=rep(c("undirected"), times=pairs))
write.table(network, file="Network.csv", sep=",", row.names=F, quote=F)
nodes <- net@peptides
table <- data.frame(ID = nodes, 
                    Type = rep(c("proteins"),
                    times = length(nodes)), 
                    Description = " ",  # Description can be added as protein name in faa file etc. (To be added afterwards)
                    Module = net@mergedColors 
                    )                    
write.table(table, file="Table.csv", sep=",", row.names=F, quote=F)
