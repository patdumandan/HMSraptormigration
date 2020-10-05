#phylogenetic signal test

require(ape)
require(phytools)
phylotree=ape::read.nexus("output.nex")
cons=averageTree(phylotree)
print(spp)
sppint=as.vector(c(3.1, -0.38, -0.81, -0.85, -0.56, 2.78, 1.71, 1.79, 2.3, 5.34, -0.4, 5.95, -0.45, 4.96, 2.15, 0.14))
sppsl1=as.vector(c(2.04, -0.72, -1.52, -2.62,-4.3, 0.75, 3.85, -0.74, 0.42, 1.29, -2.53,0.32,1.72, 0.78, -0.45,-0.98 ))
sppsl2=as.vector(c(-1.02, 1.78, 1.39, 2.01, 3.47, -0.52, -2.28, 1.04, -1.34, -1.07, 2.82, -0.99, -3.14, -1.24,-0.84, 1.05))
sppbkpt=as.vector(c(-0.23,-0.36, -0.16, -0.2, -0.23, -0.1, -0.24, -0.38, -0.09, -0.1, -0.17, -0.26, -0.1, -0.15, -0.14, -0.24))

phylosig(cons, sppbkpt, method="lambda", test=TRUE)
phylosig(cons, sppint, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list())
lambda
> # get likelihood for bm model
cons$map<-matrix(cons$edge.length, nrow(cons$edge),1,dimnames=list(NULL,"1"))
bm.logL<-brownie.lite(cons,sppint)$logL1
> # conduct hypothesis test using chi-square
LR<--2*(bm.logL-lambda$logL)
P<-pchisq(LR,df=1,lower.tail=F)
P
