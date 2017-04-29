
## ## ## concatenate output from batch runs (single simulation per run) in bash
## runs=run2S_s08
## > $runs"_all.out"
## for i in {10..59}
## do
## j=$(( $i - 9 ))
## paste -d ' ' <(cut -d ' ' -f1 $runs"-"$i".out" | sed "s/1/$j/g") <(cut -d ' ' --complement -f1 $runs"-"$i".out") >> $runs"_all.out"
## done


## ------
setwd("~/Documents/Rosalia")

runname<-"run2_a099_all"
map<-read.table("map2.txt")
Rtag<-0 ## set to 1 to summarize replicates
Ztag<-1 ## set to 0 if it wasn't calculated
bs<-10 ## binsize
nULM<-100 ## number of unlinked loci in map (nULM/bs unlinked bins)
dULM<-1.000 ## distance among bins of unlinked loci on map
dat<-read.table(paste(runname,".out",sep=""))

nloc<-dim(dat)[2]-3
nlocM<-nrow(map)

nrep<-length(unique(dat[,1]))
noutgen<-length(unique(dat[,2]))
nphase<-length(unique(dat[,3]))

pos<-map[,1]
if(nloc-nlocM > 0){ ## unlinked loci not in map
    pos<-c(map[,1],seq(map[nlocM,1]+0.1,
                       map[nlocM,1]+(nloc-nlocM)*0.1,0.1))
}
if(nULM > 0){ ## unlinked loci in map
    sp<-0.02 ## set spacer for unlinked loci in map
    ## adjust positions to have map not stretched in plot
    for(i in 1:(nULM/bs)){
        idx<-nlocM-((nULM-((i-1)*bs+1)):(nULM-(i*bs)))
        pos[idx]<-(pos[idx] - i*dULM) + sp
    }
}

selpos<-which(map[,2]>=2)
selpos<-map[selpos[order(map[selpos,2])],1]

## summarize replicates
if(Rtag==1){
    all<-array(NA,dim=c(noutgen*(3+Ztag),nlocM+3,nrep))
    for (i in 1:nrep){
        all[,,i]<-as.matrix(dat[(noutgen*(3+Ztag)*(i-1)+1):(noutgen*(3+Ztag)*i),])
    }
    tmp<-apply(all[,-c(1:3),], c(1,2), function(x) mean(x, na.rm=T))
    dat<-cbind(dat[,1:3],tmp)
}




## colors for color gradient
color.gradient<-function(colors=c("aquamarine1","aquamarine4"),
                         colsteps=50) {
    return(colorRampPalette(colors) (colsteps))
}

coltab<-rbind(c("aquamarine1","aquamarine4"),
              c("cyan1","cyan4"),
              c("coral1","coral4"),
              c("gray90","gray10"))





pdf(paste(runname,"_Fst.pdf",sep=""))
## plot Fst
## myxlim<-c(0.5,1.0)
myxlim<-range(pos)
myylim<-range(dat[seq(3,nrow(dat),3+Ztag),-c(1:3)])
par(mfrow=c(nphase,1),mar=c(3.5,3.5,0.5,0.5))
for(p in 1:nphase){
    tmp<-dat[which(dat[,1]==1 & dat[,3]==p),-c(1:3)]
    fst<-tmp[seq(3,nrow(tmp),3+Ztag),]
    subsM<-seq(ceiling(bs/2),length(pos)-nULM,bs)
    if(nULM>0){
        subsULM<-seq(nlocM-nULM+ceiling(bs/2),length(pos),bs)
    }

    mycol<-color.gradient(colors=coltab[p,],colsteps=nrow(fst))
    mylwd<-rev(log(seq(2,10,length.out=nrow(fst))))

    plot(pos[subsM],fst[1,subsM],type="n",xlab="",ylab="",
         ylim=myylim,xlim=myxlim)
    for(i in seq(1,nrow(fst),1)){
        valM<-which(!is.na(fst[i,subsM]))
        lines(pos[subsM[valM]],fst[i,subsM[valM]],lwd=mylwd[i],col=mycol[i])
        if(nULM>0){
            rect(xleft=pos[subsM[length(subsM)]]+0.006,
                 ybottom=-1, xright=pos[subsULM[1]]-0.008,
                 ytop=2,col="gray80",border=NA)
            points(pos[subsULM],fst[i,subsULM],
                   pch=16,cex=mylwd[i]/3,col=mycol[i])
        }
    }
    if(p>=nphase-1) abline(v=pos[map[,2]==1],col="red",lty=3,lwd=1.5)
    if(p<length(selpos)){
        for(i in 1:p){
            abline(v=selpos[i],col="red",lty=2,lwd=1.5)
        }
    }
    else abline(v=selpos,col="red",lty=2,lwd=1.5)

    mtext("Position on map",side=1,line=2.2,cex=0.9)
    mtext("Fst",side=2,line=2.2,cex=0.9)
}

dev.off()



if(Ztag == 1){
    pdf(paste(runname,"_Zg.pdf",sep=""))
    ## plot Zg
    ## myxlim<-c(0.5,1.0)
    myxlim<-range(pos)
    myylim<-range(dat[seq(4,nrow(dat),3+Ztag),-c(1:3)])
    par(mfrow=c(nphase,1),mar=c(3.5,3.5,0.5,0.5))
    for(p in 1:nphase){
        tmp<-dat[which(dat[,1]==1 & dat[,3]==p),-c(1:3)]
        ld<-tmp[seq(4,nrow(tmp),3+Ztag),]
        subsM<-seq(ceiling(bs/2),length(pos)-nULM,bs)
        if(nULM>0){
            subsULM<-seq(nlocM-nULM+ceiling(bs/2),length(pos),bs)
        }

        mycol<-color.gradient(colors=coltab[p,],colsteps=nrow(ld))
        mylwd<-rev(log(seq(2,10,length.out=nrow(ld))))

        plot(pos[subsM],ld[1,subsM],type="n",xlab="",ylab="",
             ylim=myylim,xlim=myxlim)
        for(i in seq(1,nrow(ld),1)){
            valM<-which(!is.na(ld[i,subsM]))
            lines(pos[subsM[valM]],ld[i,subsM[valM]],lwd=mylwd[i],col=mycol[i])
            if(nULM>0){
                rect(xleft=pos[subsM[length(subsM)]]+0.006,
                     ybottom=-1, xright=pos[subsULM[1]]-0.008,
                     ytop=2,col="gray80",border=NA)
                points(pos[subsULM],ld[i,subsULM],
                       pch=16,cex=mylwd[i]/3,col=mycol[i])
            }
        }
        if(p>=nphase-1) abline(v=pos[map[,2]==1],col="red",lty=3,lwd=1.5)
        if(p<length(selpos)){
            for(i in 1:p){
                abline(v=selpos[i],col="red",lty=2,lwd=1.5)
            }
        }
        else abline(v=selpos,col="red",lty=2,lwd=1.5)

        mtext("Position on map",side=1,line=2.2,cex=0.9)
        mtext("Zg",side=2,line=2.2,cex=0.9)
    }

    dev.off()
}




## --------------------------------------------------------------
## plot replicates for one phase, Fst only
if(Rtag==0){

    ## plot one phase only
    p<-3

    par(mfrow=c(5,2),mar=c(2.0,3.5,0.1,0.1))

    myxlim<-range(pos)
    myylim<-range(dat[seq(3,nrow(dat),3+Ztag),-c(1:3)],na.rm=T)

    for(r in 1:nrep){
        tmp<-dat[which(dat[,1]==r & dat[,3]==p),-c(1:3)]
        fst<-tmp[seq(3,nrow(tmp),3+Ztag),]
        subsM<-seq(ceiling(bs/2),length(pos)-nULM,bs)
        if(nULM>0){
            subsULM<-seq(nlocM-nULM+ceiling(bs/2),length(pos),bs)
        }

        mycol<-color.gradient(colors=coltab[p,],colsteps=nrow(fst))
        mylwd<-rev(log(seq(2,10,length.out=nrow(fst))))

        plot(pos[subsM],fst[1,subsM],type="n",xlab="",ylab="",
             ylim=myylim,xlim=myxlim)
        for(i in seq(nrow(fst),nrow(fst),1)){ ## only last
            valM<-which(!is.na(fst[i,subsM]))
            lines(pos[subsM[valM]],fst[i,subsM[valM]],lwd=mylwd[i],col=mycol[i])
            if(nULM>0){
                rect(xleft=pos[subsM[length(subsM)]]+0.006,
                     ybottom=-1, xright=pos[subsULM[1]]-0.008,
                     ytop=2,col="gray80",border=NA)
                points(pos[subsULM],fst[i,subsULM],
                       pch=16,cex=mylwd[i]/3,col=mycol[i])
            }
        }
        if(p>=nphase-1) abline(v=pos[map[,2]==1],col="red",lty=3,lwd=1.5)
        if(p<length(selpos)){
            for(i in 1:p){
                abline(v=selpos[i],col="red",lty=2,lwd=1.5)
            }
        }
        else abline(v=selpos,col="red",lty=2,lwd=1.5)

        if(r %% 5 == 0){ mtext("Position on map",side=1,line=2.2,cex=0.9)}
        mtext("Fst",side=2,line=2.2,cex=0.9)
    }



}
## --------------------------------------------------------------
## --------------------------------------------------------------
## --------------------------------------------------------------
## plot frequency of assortative mating locus
plot(dat[seq(1,nrow(dat),3+Ztag),which(map[,2]==1)+3],ylim=c(0,1))
points(dat[seq(2,nrow(dat),3+Ztag),which(map[,2]==1)+3],col="seagreen")

## plot frequency of selected locus #1
plot(dat[seq(1,nrow(dat),3+Ztag),which(map[,2]==2)[1]+3],ylim=c(0,1))
points(dat[seq(2,nrow(dat),3+Ztag),which(map[,2]==2)[1]+3],col="seagreen")

## plot frequency of selected locus #2
plot(dat[seq(1,nrow(dat),3+Ztag),which(map[,2]==2)[2]+3],ylim=c(0,1))
points(dat[seq(2,nrow(dat),3+Ztag),which(map[,2]==2)[2]+3],col="seagreen")


## --------------
## test if Zg calculations make sense

## make bins
bins<-rep(1:(length(pos)/10),each=10)

tmp<-dat[which(dat[,1]==1 & dat[,3]==1),-c(1:3)]
## last gen in phase
tmp<-tmp[(nrow(tmp)-3):nrow(tmp),]

afs1<-as.vector(t(tmp[1,]))
afs2<-as.vector(t(tmp[2,]))


calcZg<-function(afs1,afs2,nhap1=1,nhap2=1,bins){
    ## nhap1<-1 and nhap2<-1 assuming equal pop sizes

    ## set invariant loci to NA
    excl<-which(afs1+afs2 <= 0 | afs1+afs2 >= 2)
    if(length(excl)>0){
        afs1[excl]<-NA
        afs2[excl]<-NA
    }

    nhapt<-nhap1+nhap2 ## total # haplotypes
    afst<-(afs1+afs2)/2

    allbins<-unique(bins)
    nallsnps<-length(bins)

    Zg<-numeric(length(allbins))
    k<-1 ## position of snp
    n<-k

    for(i in 1:length(allbins)){ ## for all unique bins
        tmpbin<-allbins[i]
        while((k<nallsnps) & (bins[k]==bins[k+1])){
            k<-k+1
        }
        nx<-length(n:k) ## snps subset for bin
        zloc<-numeric()
        if(nx<1){
            print(paste("Couldn't find a SNP in bin ",tmpbin,sep=""))
            Zg[i]<-NA
        }
        if(nx>1){ ## don't calculate stats for bins with one variant only
            ## all pairwise LDs
            for (l in n:(k-1)){
                for (m in (l+1):k){
                    num1<-(nhap1/nhapt)*(afs1[l]*afs1[m] - afst[l]*afst[m])
                    num2<-(nhap2/nhapt)*(afs2[l]*afs2[m] - afst[l]*afst[m])
                    denom<-afst[l]*(1-afst[l])*afst[m]*(1-afst[m])
                    zloc<-c(zloc,((num1+num2)^2)/denom)
                    ## can be zero for fixed alleles or
                    ## > 1 or Inf with some rounding issues
                }
            }
            keep<-zloc[!is.na(zloc) & zloc<=1]
            nkeep<-length(keep) ## not >1 or NA
            if (nkeep>2){ ## calculate bin summary stats
                Zg[i]<-sum(keep)/nkeep ## mean
            }
        }
        k<-k+1
        n<-k
    }
    return(Zg)
}


Zg<-calcZg(afs1=as.vector(t(tmp[1,])),
           afs2=as.vector(t(tmp[2,])),
           bins=bins)

plot(as.vector(t(tmp[4,seq(1,length(pos),10)])),Zg)

plot(pos[seq(1,length(pos),10)],Zg,type="l")
lines(pos[seq(1,length(pos),10)],as.vector(t(tmp[4,seq(1,length(pos),10)])),
      col="blue")



