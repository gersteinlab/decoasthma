##--------------------------------------------------
##
## LDA-link analysis of asthma patient
##
##--------------------------------------------------




##################################################
##
## data -1.simple_cor.r  data processing
##################################################


###the last one 
wd = "" ###change here to data folder

setwd(wd)
####setwd("/ysm-gpfs/pi/gerstein/from_louise/sl2373/asthma/p152ds")


###exo.genus, genus level abundance of microbes
load("../exogenous/exo.genus.wRowName.rdata")

###clinical expression
load("../exogenous/ExoAsthma_humanGenes_clinical.RData")

###update Jan 31

load("counts.rpm.protein.rpkm.clinical.Rdata")

bulk=all.mats.protein$rpm
colnames(bulk)=gsub("\\.fq","", colnames(bulk))

exo<-exo.genus

colnames(exo)=gsub("\\.fq","", colnames(exo))
exo=exo[,colnames(bulk)]


bulk=apply(bulk, 2, function(x){ x[is.na(x)]=0.0; x;})

exo=apply(exo, 2, function(x){ x[is.na(x)]=0.0; x;})

exo=exo[-which(rowSums(exo)==0),]


##################################################
## get correlation
bce.cor=cor(t(bulk), t(exo))
bce.cor.vec=as.vector(bce.cor)

bce.cor.log2=cor(log2(t(bulk)+1), log2(t(exo)+1))
bce.corlog2.vec=as.vector(bce.cor.log2)
##17201 x  498

if (FALSE){
###Fig 4A
pdf("bulk_exo_cor_distr.pdf")
plot(density(as.vector(bce.cor)), xlab="correlation", ylab="density")
dev.off()
}


bulk.quant=apply(bulk,2, function(x){ y=rank(x); })
quartnum= floor(nrow(bulk)*0.5)

bulk.exprcnt= apply(bulk.quant, 1, function(x){ sum(x>quartnum);})
bulk.filt=bulk[which(bulk.exprcnt > 70), ]


bulk.quant2=apply(bulk,2, function(x){ y=rank(x); })
quartnum= floor(nrow(bulk)*0.5)
bulk.exprcnt2= apply(bulk.quant2, 1, function(x){ sum(x>quartnum);})

###filter exo data
exo.gt0 = apply(exo, 1, function(x){ sum(x>0);})
exo.filt=exo[which(exo.gt0>10),]   ##130x172
exo.filt2 = exo.filt[,which(colSums(exo.filt)>0)]  ##130x170


cutoff=30 ##for p152ds
bulk.filt2=bulk[which(bulk.exprcnt2 > cutoff), ]
###correlation 
bce.corfilt2=cor(t(bulk.filt2), t(exo.filt))


b.cols=split(t(bulk.filt2), col(t(bulk.filt2)))
e.cols= split(t(exo.filt), col(t(exo.filt)))



be.cortest=outer(b.cols, e.cols,Vectorize(function(x,y){cor.test(x,y)$p.value;}))


cor.vec=as.vector(bce.corfilt2)
mx = "fdr"  ###can also try BH, hochberg, hommel etc
be.cor.padj_fdr = p.adjust(as.vector(be.cortest), method = mx)
#summary(be.cor.padj_fdr[which( abs(cor.vec)> 0.4)])
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 1.369e-07 3.203e-05 2.608e-04 3.756e-04 1.599e-03



div=10
data=ceiling(bulk.filt2/div)
data=as.matrix(t(data))
data[which(data>1000)]=1000

ntopic=10
library("topicmodels")
lda.bulk10= LDA(data, k = ntopic, method = "Gibbs")
#lda.d2t.bulkfilt2= lda.noCtrl.model@gamma

save(lda.bulk10, file="bulk_ldaout10.rdata")


fin="exo_signal2_ldad10.txt"
SEED=123
data = read.table(fin, sep="\t", header=T, row.names=1, stringsAsFactors=F)

data.matrix=as.matrix(data)


lda.exo10 = LDA(data, k = ntopic, method = "Gibbs", control = list(seed = SEED, burnin = 1000,
  thin = 100, iter = 1000))


save(lda.exo10, file="exo_ldaout10.rdata")




###this is the final version

bulk.ldaout="0rpmlogF_noCtrl/bulk.filt2.lda10.txt"
exo.ldaout="0rpmlogF_noCtrl/exo.filt2.lda10.txt"

load("bulk_ldaout10.rdata")
load("exo_ldaout10.rdata")

####################################################################################################
## exogerous 
## lda.bulk10=lda.noCtrl.model from rpm1logF.rpm.10.expc10.log_FALSE1000.lda.model.rdata
##link from exogenous
load("ExoAsthma_humanGenes_clinical.RData")
sev = df[,c(1,46)]
lda.d2t.bulkfilt2= lda.bulk10@gamma
rownames(lda.d2t.bulkfilt2) = gsub("\\.", "-", lda.bulk10@documents)
sample.comm=intersect(rownames(lda.d2t.bulkfilt2), sev[,1])

sev=sev[which(sev[,1] %in% sample.comm),]
lda.d2t.bf2sev=lda.d2t.bulkfilt2[sev[,1],]

#sev = df[,c(1,46)]
lda.d2t.exofilt2= lda.exo10@gamma
rownames(lda.d2t.exofilt2) = gsub("\\.", "-", lda.exo10@documents)

lda.d2t.ef2sev=lda.d2t.exofilt2[sev[,1],]



###distribution of cor-  beta=log(prob) or ln(prob)
bulk.filt2.lda10 = lda.bulk10@beta
colnames(bulk.filt2.lda10) = lda.bulk10@terms
##write.table(bulk.filt2.lda10, file=bulk.ladout,sep="\t", quote=F, row.names=F, col.names=T)

##10x149
exo.filt2.lda10=lda.exo10@beta
colnames(exo.filt2.lda10)=lda.exo10@terms
##write.table(exo.filt2.lda10, file=exo.ldaout,sep="\t", quote=F, row.names=F, col.names=T)


###Fig 4D-F and Fig S6, S7

pdf("top20_bulk_gene2topic_dist.pdf", width=12,height=14)
par(mfrow=c(4,3))

par(cex=0.7)
par(lend=2)
par(tcl= -0.15)   
par(mar=c(3,3,3,1)+0.1)
par(mgp=c(1.1, 0.15, 0))

for ( i in 1:10){

  i.order=order(bulk.filt2.lda10[i,], decreasing=T)
  yy=exp(bulk.filt2.lda10[i,i.order[1:20]])
  xx=gsub("\\.protein_coding","",colnames(bulk.filt2.lda10)[i.order[1:20]])
  
  xbar=barplot(yy, xaxt="n", xlab="", ylab="G2T Prob",main=paste("topic",i),col="#C0C0C0")


  axis(1, labels=xx, at=xbar,cex=0.5,las=2)

  
}
dev.off()




pdf("top10_exo_microbe2topic_dist.pdf", width=12,height=14)
par(mfrow=c(4,3))

par(cex=0.7)
par(lend=2)
#par(tcl= -0.15)   
par(mar=c(7,3,3,1)+0.1)
#par(mgp=c(1.1, 0.15, 0))

for ( i in 1:10){

  i.order=order(exo.filt2.lda10[i,], decreasing=T)
  yy=exp(exo.filt2.lda10[i,i.order[1:10]])
  xx=gsub("X(\\d+)","",colnames(exo.filt2.lda10)[i.order[1:10]])
  
  xbar=barplot(yy, xaxt="n", xlab="", ylab="M2T Prob",main=paste("topic",i),col="#C0C0C0")


  axis(1, labels=xx, at=xbar,cex=0.5,las=2)

  
}
dev.off()


##########################################################################################
##                                      #
## randomForest
##########################################################################################
##define positive and negative based on correlation


##bulk.filt2.lda10=read.table(bulk.ldaout, sep="\t", header=T)
##exo.filt2.lda10=read.table(exo.ldaout, sep="\t", header=T)

cor.cut=0.4 #0.4, 0.5, 0.6, 0.7
###log2cor test #cor.cut=0.35



pos.idx =  which(abs(bce.corfilt2) > cor.cut)
pos.b= pos.idx %% nrow(bulk.filt2)
pos.b[which(pos.b == 0)] = nrow(bulk.filt2)
pos.e = ceiling(pos.idx /nrow(bulk.filt2))

neg.tmp.idx = which(abs(bce.corfilt2) < 0.05)
neg.tmp.b = neg.tmp.idx %% nrow(bulk.filt2)
neg.tmp.b[which(neg.tmp.b == 0 )] = nrow(bulk.filt2)

neg.tmp.e  = ceiling(neg.tmp.idx/ nrow(bulk.filt2))

neg.b = neg.tmp.b[which( !neg.tmp.b %in% pos.b & !neg.tmp.e %in% pos.e)]
neg.e = neg.tmp.e[which( !neg.tmp.b %in% pos.b & !neg.tmp.e %in% pos.e)]

#colnames(bce.corfilt2)=paste("X", colnames(bce.corfilt2),sep="")
neg.rand.idx = sample(length(neg.b), length(pos.b), rep=F)
pos.ds = cbind(   t(bulk.filt2.lda10[ , pos.b]) , t(exo.filt2.lda10[, pos.e]))
neg.ds = cbind( t(bulk.filt2.lda10 [, neg.b[neg.rand.idx]]) , t(exo.filt2.lda10[, neg.e[neg.rand.idx]]))

ds.xx = rbind(pos.ds, neg.ds)
ds.xx =cbind(data.frame(y=rep(c("pos","neg"), c(nrow(pos.ds), nrow(neg.ds)))), ds.xx)
ds.xx$y= as.factor(ds.xx$y)
colnames(ds.xx)[-1] = paste("X", colnames(ds.xx)[-1], sep="")

##################################################
####todo: here for cross-validation
##################################################


###where is the model for randomForest
###need to remove the overlapping idx
## pos.b unique, then remove pos.e dup, similar to negative set

pos.b.grp0=split(pos.b, as.factor(pos.e))
pos.b.grp = pos.b.grp0[sample(length(pos.b.grp0), length(pos.b.grp0), rep=F)]
pos.b.grplen= unlist(lapply(pos.b.grp, length))
pos.b.grpcumsum = cumsum(pos.b.grplen)

neg.b.grp0 = split(neg.b, as.factor(neg.e))
neg.b.grp = neg.b.grp0[sample(length(neg.b.grp0), length(neg.b.grp0),rep=F)]
neg.b.grplen =unlist(lapply(neg.b.grp, length))
neg.b.grpcumsum = cumsum(neg.b.grplen)


k0=10
psize = ceiling(sum(pos.b.grplen)/k0)
psize_neg = ceiling(sum(neg.b.grplen)/k0)

method="rf"
res = NULL
ntree=500
uprate=1
for (k in 1:k0){
  cat(k,'\r')
  ###positive
  cid = ((k-1)*psize+1):min(k*psize, sum(pos.b.grplen))
  start=which(pos.b.grpcumsum >= cid[1])[1]
  end = which(pos.b.grpcumsum <=cid[length(cid)])
  end = end[length(end)]
  
  pos.tmp = pos.b.grp[start:end]
  pos.b0id = unique(unlist(pos.tmp))
  pos.e0id = as.integer(names(pos.tmp))

  te.idx = which(pos.b %in% pos.b0id & pos.e %in% pos.e0id)
  tr.idx = which(!pos.b %in% pos.b0id & !pos.e %in% pos.e0id)
  
  pte = cbind ( t(bulk.filt2.lda10[ , pos.b[te.idx]]), t(exo.filt2.lda10[, pos.e[te.idx]]))
  ptr = cbind ( t(bulk.filt2.lda10[ , pos.b[tr.idx]]), t(exo.filt2.lda10[, pos.e[tr.idx]]))

  if(FALSE){
  te.idx = sample(length(neg.b), nrow(pte),rep=F)
  neg.b.trtmp = neg.b[-te.idx]
  neg.e.trtmp = neg.e[-te.idx]

  
  neg.b.trtmp0 = neg.b.trtmp[which(!neg.b.trtmp %in% neg.b[te.idx] & !neg.e.trtmp %in% neg.e[te.idx])]
  neg.e.trtmp0 = neg.e.trtmp[which(!neg.b.trtmp %in% neg.b[te.idx] & !neg.e.trtmp %in% neg.e[te.idx])]

  if (length(neg.b.trtmp0) < 1/2 * nrow(ptr)){
    next
  }
  tr.idx = sample(length(neg.b.trtmp0), nrow(ptr), rep=F)
  

  nte = cbind ( t(bulk.filt2.lda10[ , neg.b[te.idx]]), t(exo.filt2.lda10[, neg.e[te.idx]]))
  ntr = cbind ( t(bulk.filt2.lda10[ , neg.b[tr.idx]]), t(exo.filt2.lda10[, neg.e[tr.idx]]))
}

  ###neg
  cid = ((k-1)*psize_neg+1):min(k*psize_neg, sum(neg.b.grplen))
  start=which(neg.b.grpcumsum >= cid[1])[1]
  end = which(neg.b.grpcumsum <=cid[length(cid)])
  end = end[length(end)]
  
  neg.tmp = neg.b.grp[start:end]
  neg.b0id = unique(unlist(neg.tmp))
  neg.e0id = as.integer(names(neg.tmp))
  te.idx = which(neg.b %in% neg.b0id & neg.e %in% neg.e0id)
  tr.idx = which(!neg.b %in% neg.b0id & !neg.e %in% neg.e0id)
  
  nte = cbind ( t(bulk.filt2.lda10[ , neg.b[te.idx]]), t(exo.filt2.lda10[, neg.e[te.idx]]))
  ntr = cbind ( t(bulk.filt2.lda10[ , neg.b[tr.idx]]), t(exo.filt2.lda10[, neg.e[tr.idx]]))

  ptr = ptr[sample(nrow(ptr), ceiling(uprate*nrow(ntr)),rep=T),]
  pte = pte[sample(nrow(pte), ceiling(uprate*nrow(nte)),rep=T),]

  print(paste('nrow pte, nte', nrow(pte),nrow(nte), 'nrow ptr,ntr', nrow(ptr), nrow(ntr),'\r'))
  te = rbind(pte,nte)
  tr = rbind(ptr, ntr)
  colnames(te)=paste("X", 1:ncol(te),sep="")
  colnames(tr)=paste("X", 1:ncol(tr), sep="")


  te<-cbind(data.frame(y=rep(c(1,0), c(nrow(pte), nrow(nte)))), te)
  tr<-cbind(data.frame(y=rep(c(1,0), c(nrow(ptr), nrow(ntr)))), tr)

  te$y=as.factor(te$y)
  tr$y= as.factor(tr$y)

  tmp<-NULL
  if ( method=="rf"){
    ##stratify
    rfmodel<- randomForest(y~., data=tr, ntree=ntree) #,

  tmp<- predict(rfmodel, te[,-1],type="prob")
  tmp = data.frame(te[,1], tmp[,2])
  ##print(head(tmp))
}else if (method=="logit"){
  logist = glm(y~., data=tr, family=binomial(link='logit'))
  tmp= predict(logist, te[,-1], type="response")
  print(sum(tmp[1:nrow(pte)]>0.5)/nrow(pte))
  tmp = data.frame(te[,1], tmp) 
}else if (method=="svm"){
  svm.model=svm(y~., data=tr,probability=TRUE, scale=F)
  tmp=predict(svm.model, newdata=te[,-1], probability = TRUE)
  tmp=data.frame(te[,1], attr(tmp,"probabilities")[,2])
}else if (method=="lasso"){

  las.model=glmnet(as.matrix(tr[,-1]), tr[,1], alpha=1, family="binomial", intercept=TRUE, lambda=lbd0)
  tmp=predict(las.model, newx=as.matrix(te[,-1]), type="response")
  tmp=data.frame(te[,1], tmp)

}else if (method=="svr"){
  svr.model=svm(y~., data=tr, scale=T)
  tmp=predict(svr.model, newdata=te[,-1])
  tmp=data.frame(te[,1], tmp)
}


res=rbind(res,tmp)




}
library("PRROC")
roc1=roc.curve(scores.class0=res[which(res[,1]==1),2], scores.class1=res[which(res[,1]==0),2])$auc
prc1=pr.curve(scores.class0=res[which(res[,1]==1),2], scores.class1=res[which(res[,1]==0),2])


roc1
prc1



library("randomForest")
library("glmnet")
library("e1071")


##################################################
##cross validation

rf.b2e = randomForest(y~., data=ds.xx)
rf.b2e.impt=importance(rf.b2e)

pdf("0rpmlogF_noCtrl/rf_varimp.pdf")
varImpPlot(rf.b2e)
dev.off()



summary(as.vector(b2e.prob))

test.p = c(1:length(bce.corfilt2)) %% nrow(bulk.filt2)
test.p[which(test.p==0)] = nrow(bulk.filt2)
test.e = ceiling(c(1:length(bce.corfilt2)) / nrow(bulk.filt2))

test.ds = cbind( t(bulk.filt2.lda10[, test.p]), t(exo.filt2.lda10[, test.e]))
test.ds = as.data.frame(test.ds)

colnames(test.ds)=colnames(ds.xx)[-1]
#test.pred = predict(rf.b2e, newdata=test.ds, type="prob")
##Error in z[keep, ] <- out.class.votes :
##  NAs are not allowed in subscripted assignments
test.pred=rbind(predict(rf.b2e, newdata=test.ds[c(1:1200000),], type="prob"), predict(rf.b2e, newdata=test.ds[-c(1:1200000),], type="prob"))
b2e.prob= matrix(test.pred[,2], nrow = nrow(bce.corfilt2), ncol=ncol(bce.corfilt2),byrow=F)
#b2e.prob=as.matrix(b2e.prob)



###todo: double check 
prob.cut=0.9
#suf2="rpmlas1se"


suf2 = "corcut0.8.rpm.las1se"
b2e.links.idx=which(b2e.prob>prob.cut)  #0.95
###0.7 for lasso

b2e.links.p = b2e.links.idx %% nrow(bulk.filt2)
b2e.links.p[which(b2e.links.p ==0)]=nrow(bulk.filt2)

b2e.links.e = ceiling(b2e.links.idx / nrow(bulk.filt2))

length(unique(b2e.links.p))
length(unique(b2e.links.e))


b2e.links.c0.9=data.frame(gene=gsub("([^:]+):(.+)","\\1", rownames(bulk.filt2)[b2e.links.p]), exo=gsub("(\\d+)(.+)", "\\1:\\2", rownames(exo.filt2)[b2e.links.e]))

write.table(unique(gsub("([^:]+):(.+)","\\1", rownames(bulk.filt2)[b2e.links.p])), file=paste("0rpmlogF_noCtrl/b2e.links.c",prob.cut,"_",suf2, ".gene.txt",sep=""),sep="\t", quote=F, row.names=F, col.names=F)

write.table(b2e.links.c0.9, file=paste("0rpmlogF_noCtrl/b2e.links.probcut",prob.cut,"_",suf2,".txt",sep=""),sep="\t", quote=F, row.names=F, col.names=F)

save(list=ls(),file="rf.link.2ndRun.rdata")




###fig 4B

GAD.raw = read.table("fig4a_corlinks.gene.david.GAD.txt", 
                     header=T,stringsAsFactors=F, sep="\t")

terms = c("Asthma", "Arrhythmias, Cardiac|Long QT Syndrome", "Celiac Disease|", "Atrial Fibrillation|",
          "respiratory syncytial virus bronchiolitis", "gastrointestinal symptoms", "Bone Density",
          "Asthma|Bronchial Hyperreactivity|Hypersensitivity, Immediate", "depression | long QT syndrome",
          "Long QT Syndrome|Sinus Tachycardia|Tachycardia, Sinus" , "Sudden Infant Death")
terms.id = rev(which(GAD.raw$Term %in% terms))

pdf("fig4B.pdf", height = 4, width = 6)
par(mar=c(3,22,1,1))
xbar=barplot(-log(GAD.raw[terms.id,"PValue"]),horiz=T,col="orange", xlim=c(2,6), xpd=FALSE)
axis(2, at=xbar, labels=GAD.raw[terms.id, "Term"], las=2)
dev.off()

pdf("Fig4c.pdf")

##x1-x10: bulk 10topic; x11-x20: exo 10topics
varImpPlot(rfmodel,type=2)
dev.off()


####tool for cross-validation
ctool<-function(pdata0, ndata0, k0,lbd0, kk=1, method=c("rf","logit","svm","lasso","softmax"), ntree=500){


require("randomForest")
require("e1071")
require("PRROC")
require("ROCR")
require("glmnet")
require("AUC")


na.p= apply(pdata0[,-1],1,function(x){ y=sum(is.na(x)); y;})
na.n= apply(ndata0[,-1],1,function(x){y=sum(is.na(x)); y;})

nap.idx= which(na.p==0)
nan.idx=which(na.n==0)
aucs=matrix(rep(0,kk*3),ncol=3);
prcs=matrix(rep(0,kk*3),ncol=3);


if (colnames(pdata0)[1] != "y" || colnames(ndata0)[1] != "y"){
print("The first colnames should be the target variable and named as y")
return;
}

if(sum(colnames(pdata0) == colnames(ndata0)) != ncol(pdata0) || ncol(pdata0) != ncol(ndata0)){
print("the pdata and ndata has different number of column")
return
}


xx.col=grep("^[1-9]", colnames(pdata0[,-1]))
if(length(xx.col)>0){
print("The colnames of input file should not start with a number, change it and try again");
return;
}
rocs=NULL

if(length(nap.idx)==0 | length(nan.idx)==0){

print("No data left after removing missing data")
return
}


pdata0=pdata0[nap.idx,]
ndata0=ndata0[nan.idx,]

psize=ceiling(nrow(pdata0)/k0)
nsize=ceiling(nrow(ndata0)/k0)


res.out=NULL
for(k1 in 1:kk){  ##n time

#pdata[,1]=as.factor(pdata[,1])
#ndata[,1]=factor(ndata[,1])
pdata=pdata0[sample(nrow(pdata0),rep=F),]
ndata=ndata0[sample(nrow(ndata0),rep=F),]


res=NULL

for (k in 1:k0){  ## n fold
cat(k,'\r')

cid =((k-1)*psize+1):min(k*psize, nrow(pdata))
ptr=pdata[-cid,]
pte=pdata[cid,]

cid=((k-1)*nsize+1):min(k*nsize, nrow(ndata))
ntr=ndata[-cid,]
nte=ndata[cid,]

tr=rbind(ptr,ntr)
te=rbind(pte,nte)
####
#tr$y=as.factor(as.character(tr$y))
#te$y=as.factor(as.character(te$y))
tmp<-NULL
if ( method=="rf"){
##stratify
rfmodel<- randomForest(y~., data=tr, ntree=ntree) #, sampsize=c(nrow(ptr), nrow(ptr)), strata=tr$y)

tmp<- predict(rfmodel, te[,-1],type="prob")
tmp = data.frame(te[,1], tmp[,2])
#print(head(tmp))
}else if (method=="logit"){
logist = glm(y~., data=tr, family=binomial(link='logit'))
tmp= predict(logist, te[,-1], type="response")
print(sum(tmp[1:nrow(pte)]>0.5)/nrow(pte))
tmp = data.frame(te[,1], tmp) 
}else if (method=="svm"){
svm.model=svm(y~., data=tr,probability=TRUE, scale=F)
tmp=predict(svm.model, newdata=te[,-1], probability = TRUE)
tmp=data.frame(te[,1], attr(tmp,"probabilities")[,2])
}else if (method=="lasso"){

las.model=glmnet(as.matrix(tr[,-1]), tr[,1], alpha=1, family="binomial", intercept=TRUE, lambda=lbd0)
tmp=predict(las.model, newx=as.matrix(te[,-1]), type="response")
tmp=data.frame(te[,1], tmp)

}else if (method=="svr"){
svr.model=svm(y~., data=tr, scale=T)
tmp=predict(svr.model, newdata=te[,-1])
tmp=data.frame(te[,1], tmp)
}


res=rbind(res,tmp)


}

print(dim(res))


#print(AUC::auc(roc(res[,2], res[,1])))
###convert class as 1, 0 integer


res.out=rbind(res.out, cbind(res, data.frame(run=rep(k1,nrow(res)))))


}


out=list()
#out[["auc"]]=aucs
#out[["roc"]]=rocs

#out[["res"]]=res.out
#out
res.out
}




####end of to del
