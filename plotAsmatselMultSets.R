## plot results from asmatsel simulations for multiple settings
## ------------------------------------------------------------------

## ## ## concatenate output from batch runs (single simulation per run) in bash
## runs=run2S_s05
## > $runs"_all.out"
## for i in {10..59}
## do
## j=$(( $i - 9 ))
## paste -d ' ' <(cut -d ' ' -f1 $runs"-"$i".out" | sed "s/1/$j/g") <(cut -d ' ' --complement -f1 $runs"-"$i".out") >> $runs"_all.out"
## done


## ------

load("run2_run2S.RData")

## >>> skip this when using load()
## -----------------------------------------------------------------------
runnames<-c("run2_a09_all","run2_a099_all","run2S_s05_all","run2S_s08_all")
mapnames<-c("map2.txt","map2.txt","map2S.txt","map2S.txt")

strengthA<-c("0.9","0.99","0.0","0.0")
strengthS<-c("0.5","0.5","0.5","0.8")
xaxlablim<-c(0.0,1.5) ## specify labels to appear on x axis in plot

## color gradient for four simulations (from, to, background)
coltab<-rbind(c("gray90","#084081","white"), 
              c("gray90","#084081","white"), 
              c("gray90","#00441B","white"), 
              c("gray90","#00441B","white")) 

nsets<-length(runnames)

## assuming these basic things are identical among runs
Rtag<-1 ## set to 1 to summarize replicates
Ztag<-1 ## set to 0 if it wasn't calculated
bs<-10 ## binsize
nULM<-100 ## number of unlinked loci in map (nULM/bs unlinked bins)
dULM<-1.000 ## distance among bins of unlinked loci on map

Dat<-list("vector",nsets)
Map<-list("vector",nsets)

for (i in 1:nsets){
    Dat[[i]]<-read.table(paste(runnames[i],".out",sep=""))
    Map[[i]]<-read.table(mapnames[i])
}

## assuming the following things are identical among runs
nloc<-dim(Dat[[1]])[2]-3
nlocM<-nrow(Map[[1]])

nrep<-length(unique(Dat[[1]][,1]))
noutgen<-length(unique(Dat[[1]][,2]))
nphase<-length(unique(Dat[[1]][,3]))

pos<-Map[[1]][,1]
if(nloc-nlocM > 0){ ## unlinked loci not in map
    pos<-c(Map[[1]][,1],seq(Map[[1]][nlocM,1]+0.1,
                       Map[[1]][nlocM,1]+(nloc-nlocM)*0.1,0.1))
}
if(nULM > 0){ ## unlinked loci in map
    sp<-0.02 ## set spacer for unlinked loci in map
    ## adjust positions to have map not stretched in plot
    for(i in 1:(nULM/bs)){
        idx<-nlocM-((nULM-((i-1)*bs+1)):(nULM-(i*bs)))
        pos[idx]<-(pos[idx] - i*dULM) + sp
    }
}



Selpos<-list("vector",nsets) ## get selected positions
allNNpos<-numeric() ## get all non-neutral positions
for (i in 1:nsets){
    tmp<-which(Map[[i]][,2]>=2)
    tmp2<-which(Map[[i]][,2]>=1)
    Selpos[[i]]<-Map[[i]][tmp[order(Map[[i]][tmp,2])],1]
    allNNpos<-c(allNNpos,Map[[i]][tmp2[order(Map[[i]][tmp2,2])],1])
}
allNNpos<-sort(unique(allNNpos))
allNNnames<-LETTERS[1:length(allNNpos)]


## summarize replicates
if(Rtag==1){
    for (i in 1:nsets){
        all<-array(NA,dim=c(noutgen*(3+Ztag),nlocM+3,nrep))
        for (j in 1:nrep){
            all[,,j]<-as.matrix(Dat[[i]][(noutgen*(3+Ztag)*(j-1)+1):(noutgen*(3+Ztag)*j),])
        }
        tmp<-apply(all[,-c(1:3),], c(1,2), function(x) mean(x, na.rm=T))
        Dat[[i]]<-cbind(Dat[[i]][,1:3],tmp)
    }
}


## colors for color gradient
color.gradient<-function(colors=c("aquamarine1","aquamarine4"),
                         colsteps=50) {
    return(colorRampPalette(colors) (colsteps))
}




myxlim<-range(pos)
myylim<-c(0,0)
for (i in 1:nsets){
    myylim[2]<-max(c(myylim[2],max(Dat[[i]][seq(3,nrow(Dat[[i]]),3+Ztag),-c(1:3)])))
}

if(Ztag == 1){
    myylimZ<-c(0,0)
    for (i in 1:nsets){
        myylimZ[2]<-max(c(myylimZ[2],max(Dat[[i]][seq(4,nrow(Dat[[i]]),3+Ztag),-c(1:3)])))
    }
}




save.image("run2_run2S.RData")
## ----------------------------------------------------------------
## <<< end of skip this




## plot one phase only
p<-3

pdf(paste("Phase",p,".pdf",sep=""),width=4.5,height=3.5,pointsize=9)

par(mfrow=c(nsets,1+Ztag),mar=c(3.0,4.0,1.2,0.5))

for(i in 1:nsets){
    tmp<-Dat[[i]][which(Dat[[i]][,1]==1 & Dat[[i]][,3]==p),-c(1:3)]
    fst<-tmp[seq(3,nrow(tmp),3+Ztag),]
    subsM<-seq(ceiling(bs/2),length(pos)-nULM,bs)
    if(nULM>0){
        subsULM<-seq(nlocM-nULM+ceiling(bs/2),length(pos),bs)
    }

    mycol<-color.gradient(colors=coltab[i,1:2],colsteps=nrow(fst))
    mylwd<-rev(log(seq(2,10,length.out=nrow(fst))))
    myframecol<-coltab[i,3]
    setting<-paste(c("d = ","s = "),c(strengthA[i],strengthS[i]),sep="")

    ## plot Zg
    if(Ztag == 1){
        ld<-tmp[seq(4,nrow(tmp),3+Ztag),]
        plot(pos[subsM],ld[1,subsM],type="n",xlab="",ylab="",
             ylim=myylimZ,xlim=myxlim,axes=F)
        rect(xleft=myxlim[1]-1,xright=myxlim[2]+1,ybottom=myylimZ[1]-1,
             ytop=myylimZ[2]+1,col=myframecol,border=NA)
        box(col="black")
        axis(1,padj=-0.6,hadj=0.5,cex.axis=1.0,lwd=0,lwd.ticks=0.9,
             at=xaxlablim)
        axis(1,padj=-0.6,hadj=0.5,cex.axis=1.0,lwd=0,lwd.ticks=0.9,
             at=allNNpos[seq(1,length(allNNpos),2)],
             labels=allNNnames[seq(1,length(allNNpos),2)])
        axis(1,padj=-0.6,hadj=0.5,cex.axis=1.0,lwd=0,lwd.ticks=0.9,
             at=allNNpos[seq(2,length(allNNpos),2)],
             labels=allNNnames[seq(2,length(allNNpos),2)])
        axis(2,hadj=0.8,las=1,cex.axis=1.0,lwd=0,lwd.ticks=0.9)
        for(j in seq(1,nrow(ld),1)){
            valM<-which(!is.na(ld[j,subsM]))
            lines(pos[subsM[valM]],ld[j,subsM[valM]],lwd=mylwd[j],col=mycol[j])
            if(nULM>0){
                rect(xleft=pos[subsM[length(subsM)]]+0.006,
                     ybottom=-1, xright=pos[subsULM[1]]-0.008,
                     ytop=2,col="gray60",border=NA)
                points(pos[subsULM],ld[j,subsULM],
                       pch=16,cex=mylwd[j]/3,col=mycol[j])
            }
        }
        ## indicate barrier loci
        tmpypos<-diff(myylimZ)/20+myylimZ[1]
        if(p>=nphase-1 & sum(Map[[i]][,2]==1)) {
            ## abline(v=pos[Map[[i]][,2]==1],col="blue",lty=2,lwd=1.5)
            points(pos[Map[[i]][,2]==1],tmpypos,pch=23,lwd=1.2,col="black",bg="gray30",cex=1.0)
        }
        if(p<length(Selpos[[i]])){
            for(j in 1:p){
                ## abline(v=Selpos[[i]][j],col="red",lty=3,lwd=1.5)
                points(Selpos[[i]][j],tmpypos,pch=24,lwd=1.2,col="black",bg="white")
            }
        }
        else{
            ## abline(v=Selpos[[i]],col="red",lty=3,lwd=1.5)
            points(Selpos[[i]],rep(tmpypos,length(Selpos[[i]])),pch=24,lwd=1.2,col="black",bg="white")
        }

        ## legend
        if(i == 1){
            ## legend("topleft",legend=c("Assortative mating locus","Locus under divergent selection"),col=c("blue","red"),lty=c(2,3),y.intersp=0.9,bty="n")
            legend("topleft",legend=c("Assortative mating locus","Locus under divergent selection"),col=c("black","black"),pch=c(23,24),pt.bg=c("gray30","white"),pt.lwd=c(1.2,1.2),pt.cex=c(1.0,1.0),y.intersp=0.9,bty="n")
        }

        ## strengths of d and s
        mtext(paste(setting[1],setting[2],sep=", "),cex=0.8)

        if(i == nsets){
            mtext("Position on map [M]",side=1,line=1.8,cex=0.8)
        }
        mtext(expression('Z'[g]),side=2,line=2.4,cex=0.8)
    }

    ## plot Fst
    plot(pos[subsM],fst[1,subsM],type="n",xlab="",ylab="",
         ylim=myylim,xlim=myxlim,axes=F)
    rect(xleft=myxlim[1]-1,xright=myxlim[2]+1,ybottom=myylimZ[1]-1,
         ytop=myylimZ[2]+1,col=myframecol,border=NA)
    box(col="black")
    axis(1,padj=-0.6,hadj=0.5,cex.axis=1.0,lwd=0,lwd.ticks=0.9,
         at=xaxlablim)
    axis(1,padj=-0.6,hadj=0.5,cex.axis=1.0,lwd=0,lwd.ticks=0.9,
         at=allNNpos[seq(1,length(allNNpos),2)],
         labels=allNNnames[seq(1,length(allNNpos),2)])
    axis(1,padj=-0.6,hadj=0.5,cex.axis=1.0,lwd=0,lwd.ticks=0.9,
         at=allNNpos[seq(2,length(allNNpos),2)],
         labels=allNNnames[seq(2,length(allNNpos),2)])
    axis(2,hadj=0.7,las=1,cex.axis=1.0,lwd=0,lwd.ticks=0.9,
         at=seq(0.0,0.6,0.2))
    for(j in seq(1,nrow(fst),1)){
        valM<-which(!is.na(fst[j,subsM]))
        lines(pos[subsM[valM]],fst[j,subsM[valM]],lwd=mylwd[j],col=mycol[j])
        if(nULM>0){
            rect(xleft=pos[subsM[length(subsM)]]+0.006,
                 ybottom=-1, xright=pos[subsULM[1]]-0.008,
                 ytop=2,col="gray60",border=NA)
            points(pos[subsULM],fst[j,subsULM],
                   pch=16,cex=mylwd[j]/3,col=mycol[j])
        }
    }
    ## indicate barrier loci
    tmpypos<-diff(myylim)/20+myylim[1]
    if(p>=nphase-1 & sum(Map[[i]][,2]==1)) {
        ## abline(v=pos[Map[[i]][,2]==1],col="blue",lty=2,lwd=1.5)
        points(pos[Map[[i]][,2]==1],tmpypos,pch=23,lwd=1.2,bg="black",col="gray30",cex=1.0)
    }
    if(p<length(Selpos[[i]])){
        for(j in 1:p){
            ## abline(v=Selpos[[i]][j],col="red",lty=3,lwd=1.5)
            points(Selpos[[i]][j],tmpypos,pch=24,lwd=1.2,col="black",bg="white")
        }
    }
    else{
        ## abline(v=Selpos[[i]],col="red",lty=3,lwd=1.5)
        points(Selpos[[i]],rep(tmpypos,length(Selpos[[i]])),pch=24,lwd=1.2,col="black",bg="white")
    }

    ## strengths of d and s
    mtext(paste(setting[1],setting[2],sep=", "),cex=0.8)

    if(i == nsets){
        mtext("Position on map [M]",side=1,line=1.8,cex=0.8)
    }
    mtext(expression(paste(italic('F')[ST])),side=2,line=2.2,cex=0.8)
}

dev.off()
