library(motmot)
library(BSDA)

pTreeText <- "(((oryLat04:0.51012,(((parOli02:0.160689,cynSem:0.433672)0:0.0437825,
              serDum01:0.107517)0:0.0352359,monAlb01:0.188107)0:0.0202366)0:0.0117874,(gasAcu14:
              0.324389,(((tetNig2:0.460829,molMol01:0.180985)1:0.0734705,larCro01:0.15854)0:
              0.00820025,dicLab01:0.119023)0:0.018837)0:0.0188741)0:0.116492,(synSco01:0.170947,
              (hipEre01:0.0331536,hipCom01:0.0298633)1:0.12246)1:0.392849)0;"

pTree <- read.tree(text = pTreeText)
finData <- read.table("finData.tab", sep="\t", header=TRUE)
row.names(finData) <- finData$species

pRateMatrix <- as.rateMatrix(phy=pTree, x="split", data=finData)
pCaudalRateData <- as.rateData(y="num_caudal_rays", x="split", 
                               rateMatrix=pRateMatrix, phy=NULL, data=finData, log.y=FALSE)
groupMeans <- phyloMean(rateData=pCaudalRateData, common.mean = FALSE)

outgroupMean <- groupMeans[1,1]
targetMean <- sum(groupMeans)

targetData <- finData[finData$phenotypic_group == "target",]
targetData$squares <- (targetData$num_caudal_rays - targetMean)^2
targetNminus1 <- length(targetData$species)-1
targetSD <- sqrt(sum(targetData$squares)/targetNminus1)

outgroupData <- finData[finData$phenotypic_group == "outgroup",]
outgroupData$squares <- (outgroupData$num_caudal_rays - outgroupMean)^2
outgroupNminus1 <- length(outgroupData$species)-1
outgroupSD <- sqrt(sum(outgroupData$squares)/outgroupNminus1)

# run Welch-modified T Test
tsum.test(mean.x = outgroupMean, s.x = outgroupSD, n.x=length(outgroupData$species), 
          mean.y = targetMean, s.y = targetSD, n.y = length(targetData$species))
