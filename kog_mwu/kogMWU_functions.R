
mwut=function(gos,Alternative) {
	terms=unique(gos$term)
	gos$seq=as.character(gos$seq)
	nrg=gos[!duplicated(gos$seq),2]
	names(nrg)=gos[!duplicated(gos$seq),1]
#	nrg=nrg+rnorm(nrg,sd=0.01) # to be able to do exact wilcox test
	rnk=rank(nrg)
	names(rnk)=names(nrg)
	pvals=c();drs=c();nams=c();levs=c();nseqs=c()
	for (t in terms){
		got=gos[gos$term==t,]
		got=got[!duplicated(got$seq),]
		ngot=gos[gos$term!=t,]
		ngot=ngot[!duplicated(ngot$seq),]
		ngot=ngot[!(ngot$seq %in% got$seq),]
		sgo.yes=got$seq
		n1=length(sgo.yes)
		sgo.no=ngot$seq
		n2=length(sgo.no)
		if (n2 < n1) {
			print(paste("skipping",t,"nseqs =",n1))
			next
		}
		wi=wilcox.test(nrg[sgo.yes],nrg[sgo.no],alternative=Alternative)	# removed correct=FALSE 
		r1=sum(rnk[sgo.yes])/n1
		r0=sum(rnk[sgo.no])/n2
		dr=r1-r0
		drs=append(drs,round(dr,0))
		levs=append(levs,got$lev[1])
#		nams=append(nams,as.character(got$term[1]))
		pvals=append(pvals,wi$p.value)
		nseqs=append(nseqs,n1)	
	}
	res=data.frame(cbind(nseqs,"delta.rank"=drs,"pval"=pvals))
	res=cbind("term"=as.character(terms),res)
	res$pval=as.numeric(as.character(res$pval))
	res$delta.rank=as.numeric(as.character(res$delta.rank))
#	res$level=as.numeric(as.character(res$level))
	res$nseqs=as.numeric(as.character(res$nseqs))
	res=res[order(res$pval),]
	res$padj=p.adjust(res$pval,method="BH")
	return(res)
}

ft=function(gos) {
gos=annotated
	terms=unique(gos$term)
	gos$seq=as.character(gos$seq)
	pft=c();nam=c();lev=c();nseqs=c()
	for (t in terms) {
		got=gos[gos$term==t,]
		got=got[!duplicated(got$seq),]
		ngot=gos[gos$term!=t,]
		ngot=ngot[!duplicated(ngot$seq),]
		ngot=ngot[!(ngot$seq %in% got$seq),]
		go.sig=sum(got$value)
		go.ns=length(got[,1])-go.sig
		ngo.sig=sum(ngot$value)
		ngo.ns=length(ngot[,1])-ngo.sig
		sig=c(go.sig,ngo.sig) # number of significant genes belonging and not belonging to the tested GO category
		ns=c(go.ns,ngo.ns) # number of not-significant genes belonging and not belonging to the tested GO category
		mm=matrix(c(sig,ns),nrow=2,dimnames=list(ns=c("go","notgo"),sig=c("go","notgo")))
		ff=fisher.test(mm,alternative="greater")
		pft=append(pft,ff$p.value)
		nam=append(nam,as.character(got$name[1]))
		lev=append(lev,got$lev[1])
		nseqs=append(nseqs,length(got[,1]))
	}
	res=data.frame(cbind("term"=as.character(terms),nseqs,"pval"=pft))
	res$pval=as.numeric(as.character(res$pval))
	res$nseqs=as.numeric(as.character(res$nseqs))
	res=res[order(res$pval),]
	res$padj=p.adjust(res$pval,method="BH")
	return(res)
}

#---------------------

# gene2kog="~/Documents/annotations_april2014/amillepora_2014/amil_apr2014_iso2kogClass.tab"
# inname="~/Documents/bay_davies_matz_mortality2012/larvae_logPvalSigned_apr2014.csv"
# Alternative="t"

kog.mwu=function(inname,gene2kog,Alternative) {
	rsq=read.csv(inname)
	names(rsq)=c("seq","value")
	bads=which(rsq$value==Inf | rsq$value==(-Inf) | is.na(rsq$value))
	if (length(bads)>0) { rsq=rsq[-bads,]}
	kogs=read.table(gene2kog,sep="\t")
	annotated=rsq[rsq$seq %in% kogs$V1,]
	kogs=kogs[kogs$V1 %in% rsq$seq,]
	
	kogrows=match(kogs$V1,annotated$seq)
	annotated=annotated[kogrows,]
	annotated$term=as.character(kogs$V2)
	
	mwut.t=TRUE
	if (length(levels(as.factor(annotated$value)))==2) {
		print("Binary classification detected; will perform Fisher's test");
		mwut.t=F
		rr=ft(annotated)
	} else {
		print("Continuous measure of interest: will perform MWU test");		
		rr=mwut(annotated,Alternative)
	}
	return(rr)
}
