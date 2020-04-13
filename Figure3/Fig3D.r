



##Fig3D  Fig S5
exo.corfilt2=cor(t(exo.filt2))
exo.corfilt2.lda=cor(exo.filt2.lda10)

# Make an Igraph object from this matrix:

plot(network)

library(RColorBrewer)
coul = brewer.pal(nlevels(as.factor(mtcars$cyl)), "Set2")
# Map the color to cylinders
my_color=coul[as.numeric(as.factor(mtcars$cyl))]
 
# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)
plot(network,
    vertex.size=12,
    vertex.color=my_color,
    vertex.label.cex=0.7,
    vertex.label.color="white",
    vertex.frame.color="transparent"
    )
# title a



rownames(exo.corfilt2)=gsub("\\d+:","", rownames(exo.corfilt2))
colnames(exo.corfilt2)=gsub("\\d+:","", colnames(exo.corfilt2))

rownames(exo.corfilt2.lda)=gsub("X\\d:","", rownames(exo.corfilt2.lda))
colnames(exo.corfilt2.lda)=gsub("X\\d+","", colnames(exo.corfilt2.lda))



#topn = 10

toprank.lda = t(apply(exo.filt2.lda10, 1 , function(x,y){x=as.numeric(x); yy=y[order(x,decreasing=T)[1:10]]; yy;}, y=gsub("X\\d+", "", colnames(exo.filt2.lda10))))

toprank.cut = t(apply(exo.filt2.lda10, 1 , function(x,y){x=as.numeric(x); xx=x[order(x,decreasing=T)[1:10]]; xx;}, y=gsub("X\\d+", "", colnames(exo.filt2.lda10))))

top.dup=as.vector(toprank.lda)[duplicated(as.vector(toprank.lda))]
top.uniq=as.vector(toprank.lda)[-which(duplicated(as.vector(toprank.lda)))]


member=rep("",length(unique(as.vector(toprank.lda))))
names(member)=unique(as.vector(toprank.lda))
for(i in 1:10){
  member[as.character(toprank.lda[i,])] = paste(member[as.character(toprank.lda[i,])], i,sep=",")

}

member=gsub("^,","",member)

for (x in names(member)){
  asmember[x]

}

toprank.lda
##cluster edge betweenness
cutoff = 0.3

exo.cf2.nw =abs( exo.corfilt2)
exo.cf2.nw[exo.cf2.nw < cutoff]=0
cnet=graph_from_adjacency_matrix( exo.cf2.nw, weighted=T, mode="undirected", diag=F)
V(cnet)$label=V(cnet)$name
                                        #l =
cnet=delete.vertices(cnet, degree(cnet)==0)
V(cnet)[name %in% names(member)]$label= paste(V(cnet)[name %in% names(member)]$label,"(", member[V(cnet)[name %in% names(member)]$name],")",sep="")


#ceb = cluster_edge_betweenness(cnet)
ceb = cluster_label_prop(cnet)
#ceb=cluster_fast_greedy(cnet)

plot(ceb, cnet, edge.arrow.mode=0,  vertex.label=V(cnet)$label, vertex.label.cex=0.5, vertex.size=1) 

dev.copy2pdf(file="microbe.corRAW.network.c0.2lp.pdf")




cutoff=0.8

 exo.cf2.lda = abs(exo.corfilt2.lda)
exo.cf2.lda[exo.cf2.lda < cutoff]=0
cnet2=graph_from_adjacency_matrix( exo.cf2.lda, weighted=T, mode="undirected", diag=F)
 V(cnet2)$label=V(cnet2)$name

#g=cnet                                      #l =
cnet2=delete.vertices(cnet2, degree(cnet2)==0)
V(cnet2)[name %in% names(member)]$label= paste(V(cnet2)[name %in% names(member)]$label,"(", member[V(cnet2)[name %in% names(member)]$name],")",sep="")

ceb = cluster_edge_betweenness(cnet2)
#ceb = cluster_label_prop(cnet2)
#ceb=cluster_fast_greedy(cnet)

plot(ceb, cnet2, edge.arrow.mode=0,  vertex.label=V(cnet2)$label, vertex.label.cex=0.5, vertex.size=1) 
dev.copy2pdf(file="microbe.corLDA.network.c0.8.pdf")

