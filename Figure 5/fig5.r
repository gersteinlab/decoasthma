
##fig5A

all.links=read.table("b2e.links_rf.txt",sep="\t", header=F, stringsAsFactors=F)
load("fig5_pin2gene.rdata")
all.links=all.links[order(all.links[,2]),]
#rownames(all.links)=gsub("::","", all.links[,2])
all.links=cbind(all.links, gsub("::","", all.links[,2]))
colnames(all.links)[4]="V4"
all.links.rle=rle(as.character(all.links[,4]))
all.links.rle=data.frame(val=all.links.rle$values,len= all.links.rle$lengths)



g.ds=all.links[grep("Haemophilus|Bifidobacterium|Candida|Shigella|Saccharomyces|Lactobacillus|Aeromonas|Pasteurella|Proteus|Burkholderia",all.links$V4),]
g.ds$V4=as.character(g.ds$V4)
g.ds$V4=gsub("^\\d+","", g.ds$V4)

g= graph.data.frame(g.ds[,c(1,4)],directed=F)
V(g)$type=1
V(g)[name %in% g.ds[,1]]$type=2
V(g)$label=V(g)$name


l <- layout_with_sugiyama(g, V(g)$type)
#l$layout[,1] = (l$layout[,1] - min(l$layout[,1]))

pdf("fig5_kegg_bipitite.pdf",height=4,width=8)
par(mar=c(1,0,1,0))
plot.igraph(g, layout=-l$layout[,1:2], asp=0, margin=0,

     vertex.color = xx.rcol[V(g)$type],
     vertex.shape = "none",
      vertex.label.color=xx.rcol[V(g)$type],
	vertex.label.dist=0,	
	vertex.label.cex=0.6,srt=315,
           # rescale=FALSE,xlim=range(l$layout[,1]), ylim=range(l$layout[,2])
           # rescale=TRUE
)
dev.off()



library("igraph")
library("rTRM")

g1 = graph.data.frame(bot.b2e.linkds, directed=F)
g2= graph.data.frame(msign.path,directed=F)
g = g1+g2
    xx.rcol=c("gray","orange","steelblue")
    xx.shape=c("circle","circle", "circle")

##pathway
V(g)$type=1
##gene
V(g)[name %in% c(
"TAGAP","IL1B","NFKBIZ","DAW1","LCE3D","SBSN","SPRR2F","OSBPL6","DCLK2","GCSAML","STK32A","KIAA1024L","NPY1R","SEC14L6","SPATA31D1","GFAP","RNASE9","HES7","PROK1","SKOR2","HORMAD1","SCN5A","SPATA31D1","RNASE9","GFI1","HES7","IRX3","PROK1","RAB40A","GNB5","GSTM2","AC006026.2","CDC34","NAALADL2","LTF","CHIT1","SARDH","STK32A","VASH2","SAMD14","TAGAP","IL1B","MUC6","ALKBH6","DAAM2","DOK6","BCL2L14","HESX1","NPY1R","SEC14L6","TCN1","KRT24","NPY1R","SEC14L6","GNLY","EGFR","GPR15","BCL2L14","HESX1","TAGAP","IL1B","KCNJ2","CACNA1E","SCN1B","HPN","CSF3","FAM13C","B4GALNT3","GABRG3","FAM174B","NPY1R","SEC14L6","KIAA1024L","GABRG3","STK32A","MKRN3","ABCA4","IMMP2L","GABRG3","KLHL13","RAB40A"

)]$type=2

V(g)[name %in% c(
"Haemophilus","Haemophilus","Haemophilus","Janthinobacterium","Lactobacillus","Lactobacillus","Lactobacillus","Propionibacterium","Propionibacterium","Candida","Gemella","Rothia","Megasphaera","Megasphaera","Megasphaera","Megasphaera","Megasphaera","Megasphaera","Megasphaera","Megasphaera","Bifidobacterium","Bifidobacterium","Bifidobacterium","Bifidobacterium","Bifidobacterium","Bifidobacterium","Bifidobacterium","Bifidobacterium","Enterobacter","Cronobacter","Cronobacter","Cronobacter","Cronobacter","Cronobacter","Aeromonas","Dietzia","Dietzia","Streptomyces","Streptomyces","Shigella","Brenneria","Brenneria","Burkholderia","Azotobacter","Xanthomonas","Neorhizobium","Ensifer","Ensifer","Cellulosimicrobium","Cellulosimicrobium","Mannheimia","Mannheimia","Exiguobacterium","Exiguobacterium","Catonella","Proteus","Proteus","Sulfurospirillum","Sulfurospirillum","Pasteurella","Pasteurella","Pasteurella","Pasteurella","Pasteurella","Pasteurella","Pasteurella","Dickeya","Cupriavidus","Cupriavidus","Geobacter","Geobacter","Geobacter","Geobacter","Geobacter","Geobacter","Geobacter","Trabulsiella","Methylibium","Methylibium","Anaeroglobus","Saccharomyces"

)]$type=3
V(g)$label=V(g)$name




l <- layout_with_sugiyama(g, V(g)$type)
#l$layout[,1] = (l$layout[,1] - min(l$layout[,1]))

pdf("fig5_kegg_tripitite.pdf",height=8,width=4)
par(mar=c(1,0,1,0))
plot.igraph(g, layout=-l$layout[,1:2], asp=0, margin=0,

     vertex.color = xx.rcol[V(g)$type],
     vertex.shape = "none",
      vertex.label.color=xx.rcol[V(g)$type],
	vertex.label.dist=0,	
	vertex.label.cex=0.6,srt=315,
           # rescale=FALSE,xlim=range(l$layout[,1]), ylim=range(l$layout[,2])
           # rescale=TRUE
)
dev.off()






      if(xg %in% rownames(lm22)){
###draw lm figure
        #xx.lm22.expr=lm2[xg,]
        #xx.lm22.freq=
      }





#####figure 5b

##genepool=unique(all.links[,1])
genepool=c("IL1B","CACNA1E","GCSAML")
microbepool=c("Haemophilus","Pasteurella","Candida")

ncut=2000
date0="180718."
for (i in 1:nrow(all.links.rle)){
  if (all.links.rle[i, 2] < ncut){  ###only keep the 
    ##3mm = as.character(all.links.rle[i, 1])  ##microbe
    ##
    mm=
    cat(mm,'\r')
    mm.expr=exo[mm,cc]
#### mm >0 ; mm
    mm.up=names(mm.expr[which(mm.expr>0)])
    mm.dw=names(mm.expr[which(mm.expr==0)])


    xx.genes=all.links[which(all.links$V4 == mm ),1]
    xx.genes=xx.genes[which(xx.genes %in% genepool)]
### gene list
    if(length(mm.up)>0 & length(mm.dw)>0){
      for(xg in xx.genes ){
        xx.bulk=bulk[xg,cc] ###xg
        ##nm=paste(mm,".", xg,sep="")

 
        xx.glm22=NULL
        xx.lm22.up=NULL
        xx.lm22.dw=NULL
        pv.list=NULL
        xx.order=NULL
        xx.lm22.list=NULL
      if(xg %in% rownames(lm22)){
###draw lm figure
        #xx.lm22.expr=lm2[xg,]
        #xx.lm22.freq=
        ###sum( colnames(lm22)==rownames(lm22freq)) ==22
###freq >0 and expr > 0, top7 significance
        xx.glm22=lm22[xg,]

        xx.lm22.up = split(lm22freq[, mm.up[which(mm.up %in% colnames(lm22freq))] ], seq(nrow(lm22freq)))
        xx.lm22.dw = split(lm22freq[, mm.dw[which(mm.dw %in% colnames(lm22freq))] ], seq(nrow(lm22freq)))

        pv.list=as.numeric(mapply(function(x,y){ sprintf("%.3f",wilcox.test(x,y)$p.value);},xx.lm22.up, xx.lm22.dw))

        xx.order=order(pv.list)[1:7]

        xx.lm22.list=c(xx.lm22.up[xx.order], xx.lm22.dw[xx.order])[c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]

        
      }else if(xg %in% rownames(lm22.deconv)){

######need positive value
        ##sum(colnames(lm22.deconv) == rownames(lm22freq))==22
        xx.glm22=lm22.deconv[xg,]

        xx.lm22.up = split(lm22freq[, mm.up[which(mm.up %in% colnames(lm22freq))] ], seq(nrow(lm22freq)))
        xx.lm22.dw = split(lm22freq[, mm.dw[which(mm.dw %in% colnames(lm22freq))] ], seq(nrow(lm22freq)))

        pv.list=as.numeric(mapply(function(x,y){ sprintf("%.3f",wilcox.test(x,y)$p.value);},xx.lm22.up, xx.lm22.dw))

        xx.order=order(pv.list)[1:7]

        xx.lm22.list=c(xx.lm22.up[xx.order], xx.lm22.dw[xx.order])[c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]




      }

                ####modifiy here 
        pdf(paste(date0, ".lm22.", mm,".", xg, ".pdf",sep=""), width=6, height=6)

        par(cex=1)
        par(lend=2)
        par(tcl= -0.15)   
        par(mar=c(8.5,3,3,2)+0.1)
        par(mgp=c(1.1, 0.15, 0))
        boxplot(xx.lm22.list, xaxt="n", outline=F, at=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12, 13.5, 14.5, 16,17),col=rep(c("gray","cyan"),7),xlim=c(0,18),ylab="Freq")

        par(new = T)
        pv.list2=pv.list[xx.order]
        xx.glm22=xx.glm22[xx.order]
      
        xx.lm22.expr=log2(xx.glm22+2)
        plot(x=c(1.5, 4, 6.5, 9, 11.5, 14, 16.5)+1.1,y=xx.lm22.expr,xlim=c(0,18), xlab="", xaxt='n', yaxt='n',ylab="",ylim=range(xx.lm22.expr)*c(0.8,1.2), pch=22, bg="red",main=paste("lm22", gsub("^\\d","",mm), xg,sep=" "))
        axis(4, at=seq(as.integer(min(xx.lm22.expr)), as.integer(max(xx.lm22.expr)),length.out=5),labels=seq(as.integer(min(xx.lm22.expr)), as.integer(max(xx.lm22.expr)),length.out=5))
        mtext("Log(Expr+2)", side=4, line=1, cex.lab=1,las=0, col="black")

        axis(1, at=c(1.5, 4, 6.5, 9, 11.5, 14, 16.5), labels= FALSE) # rownames(lm22freq)[xx.order],las=2)
        text(x=c(0.5, 1.5, 4, 6.5, 9, 11.5, 14, 16.5), y=max(xx.lm22.expr),labels=c("P:",pv.list2),col="orange",cex=1)
        legend("topright", c("with microbe", "No microbe", "Expr"), fill=c("gray", "cyan", "red"))

        text(x=c(1.5, 4, 6.5, 9, 11.5, 14, 16.5), y=par()$usr[3]-0.02*(par()$usr[4]-par()$usr[3]),labels=rownames(lm22freq)[xx.order], srt=45, adj=1, xpd=TRUE)
  
        dev.off()
        ### enod of modify
      


      }
    


    }
  }

}



