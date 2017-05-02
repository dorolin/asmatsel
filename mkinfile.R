## create file that gives marker positions and positions of loci under
##   divergent selection and involved in assortative mating for 'asmatsel'

## locus controlling assortative mating must be tagged with '1'
## loci under divergent selection must be tagged with
##   2:(1 + number of loci under selection)

## ---------------------------------------------------------------------

## map for Zg and using bins (too large to run off cluster)
map2<-seq(0.006,1.505,0.001)
## add unlinked loci
for(i in 1:10){
    map2<-c(map2,seq(max(map2)+1.001,max(map2)+1.010,0.001))
}

map2<-cbind(map2,rep(0,length(map2)))

map2[map2[,1]==0.65,2]<-1 ## assortative mating locus
map2[map2[,1]==0.75,2]<-2 ## selected locus 1
map2[map2[,1]==1.25,2]<-3 ## selected locus 2

## plot(map2[,1],type="l")
## abline(v=which(map2[,2]==1),lty=3)
## abline(v=which(map2[,2]>=2),lty=2)

write.table(map2,file="map2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)


## -------------------------------

## shorter test map for Zg and using bins
map2c<-seq(0.606,1.305,0.001)

map2c<-cbind(map2c,rep(0,length(map2c)))

map2c[map2c[,1]==0.65,2]<-1 ## assortative mating locus
map2c[map2c[,1]==0.75,2]<-2 ## selected locus 1
map2c[map2c[,1]==1.25,2]<-3 ## selected locus 2

## plot(map2c[,1],type="l")
## abline(v=which(map2c[,2]==1),lty=3)
## abline(v=which(map2c[,2]>=2),lty=2)

write.table(map2c,file="map2c.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)


## -------------------------------

## map for Zg and using bins (too large to run off cluster)
## no assortative mating locus but instead third locus under selection
map2S<-seq(0.006,1.505,0.001)
## add unlinked loci
for(i in 1:10){
    map2S<-c(map2S,seq(max(map2S)+1.001,max(map2S)+1.010,0.001))
}

map2S<-cbind(map2S,rep(0,length(map2S)))

map2S[map2S[,1]==0.65,2]<-4 ## selected locus 3
map2S[map2S[,1]==0.75,2]<-2 ## selected locus 1
map2S[map2S[,1]==1.25,2]<-3 ## selected locus 2

## plot(map2S[,1],type="l")
## abline(v=which(map2S[,2]==1),lty=3)
## abline(v=which(map2S[,2]>=2),lty=2)

write.table(map2S,file="map2S.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)


## -------------------------------

## shorter test map for Zg and using bins
## no assortative mating locus but instead third locus under selection
map2cS<-seq(0.606,1.305,0.001)

map2cS<-cbind(map2cS,rep(0,length(map2cS)))

map2cS[map2cS[,1]==0.65,2]<-4 ## selected locus 3
map2cS[map2cS[,1]==0.75,2]<-2 ## selected locus 1
map2cS[map2cS[,1]==1.25,2]<-3 ## selected locus 2

## plot(map2cS[,1],type="l")
## abline(v=which(map2cS[,2]==1),lty=3)
## abline(v=which(map2cS[,2]>=2),lty=2)

write.table(map2cS,file="map2cS.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

