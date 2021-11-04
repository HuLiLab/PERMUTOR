#REQUIRED PACKAGES 
library(plyr)
library(igraph)
library(rlist)
library(rapport)
library(stringi)
library(stringr)

'%ni%' <- Negate('%in%')

# LOAD DATA
#absolute fold change per patient |log2(diseased/normal)|
fcrnaseq=read.csv("FCRNAseq.csv",row.names = 1)
# matrix indicating which genes are mutated in which patients 0= non-mutated 1=mutated
patmut=read.csv("patmutations.csv",row.names = 1)
# disease context information- higher numbers indicate stronger disease drivers 
context=read.csv("diseasecontext.csv",row.names = 1)
# starting protien-protein interaction network in edgelist format
genedgelist=read.csv("generalppi.csv",row.names = 1)   
# target therapy gene targets 
ttherapy=as.vector(unlist(read.csv("target_therapy_genes.csv",row.names = 1)))

#SETTING WEIGHTS
# path score
const1 <- 1 #RNAseq FC weight
const2<- 10 #tier 1 disease gene weight
const3 <- 2 # tier 2 disease gene weight
const4 <- 0 # not a diseased gene
const5 <- 5 # mutation weight

# P-VALUE CUT OFF FOR SIGNIFICANT PATHS
pvalcut=0.01

# IDGI
const6<-1 #edge weight
const7<-10 #subnetwork weight 
const8<-5 #hub weight
const9<-6 #highly trafficked hub weight
const10<-3 #bottleneck weight

#PREPROCESSING
ppigenes <- unique(c(genedgelist[,1],genedgelist[,2]))
ppi <- graph_from_edgelist(as.matrix(genedgelist), directed= FALSE)

#PATIENT SELECTION-include the column number of the patients you want to run below
selected<-c(2)

#PER PATIENT LOOP
for (p in selected) {
  # initialize storage variables 
  allrandscores<-vector()
  allrandpaths<-list()
  patient<-list()
  
  patientcharacteristics<-data.frame(matrix(nrow=1,ncol=9))
  colnames(patientcharacteristics)<-c("Mutated Genes","Ind Network Edges","Ind Network Nodes","Ind Network Mutated Genes","Ind Network Components","Hubs","Highly Trafficked Hubs", "Bottlenecks", "Highly Trafficked Bottlenecks")
  patientcharacteristics[1,1]<-sum(patmut[,p])
  mutlist<-rownames(patmut)[which(patmut[,p]==1)]
  mut<-combn(mutlist,2)
  

  #FINDING ALL SHORTEST PATHS
  for (a in 1:length(mut[1,])) {
    pathlist<-all_shortest_paths(ppi,mut[1,a],mut[2,a])
    #make sure all of these paths are valid paths 
    if(length(pathlist[[1]])!= 0){
      #loop through all shortest paths between each gene pair and add them to patient list 
      for (b in 1:length(pathlist[[1]])) {
        patient<- list.append(patient,names(pathlist[[1]][[b]]))
        names(patient)[[length(patient)]]<-paste(mut[1,a],mut[2,a],sep=" ")
      }
    }
  }
  save(patient,file=paste(p,"allpatientpaths.Rdata"))
  
  #GENERATE RANDOM PATHS 
  #find max and min path length 
  maxl<-max(unlist(lapply(patient,length)))
  minl<-min(unlist(lapply(patient,length)))
  # number of random paths to be generated per patient per length
  randnum<-1000
  for (c in minl:maxl) {
    set.seed(as.double(paste0(p,c)))
    tempscore<-vector(mode="integer",length = randnum)
    templist<-vector("list",length=randnum)
    #rand path creation and scoring loop
    for (d in 1:randnum) {
      newpath<-vector(mode="character",length = c)
      newpath[1]<-sample(mutlist,1)
      newpath[c]<-sample(mutlist,1)
      if("" %in% newpath) {for(e in 2:(c-1)){newpath[e]<-sample(ppigenes,1)}}
      templist[[d]]<-unlist(newpath)
      #initialize scoring variables 
      rtotal<-0
      for (f in 1:length(newpath)) {
        gene<-newpath[[f]]
        if( gene %in% rownames(fcrnaseq)){rtotal<-rtotal+abs(const1*fcrnaseq[unlist(gene),p])}
        if( gene %in% rownames(context)){
          if(context[unlist(gene),1]==1){rtotal<-rtotal+const2}
          else if(context[unlist(gene),1]==2){rtotal<-rtotal+const3}
        }
        else {rtotal<-rtotal+const4}
        if(gene %in% mutlist & f != 1 & f != length(newpath)) {rtotal<-rtotal+const5}
      }
      rtotal<-rtotal/length(newpath)
      tempscore[d]<-rtotal
    }
    assign(paste("randscore",as.character(c)),tempscore)
    assign(paste("randlist",as.character(c)),templist)
    save(tempscore,file = paste(as.character(p),"randscore",as.character(c),".Rdata"))
    allrandscores<-c(allrandscores,tempscore)
    allrandpaths<-list.append(allrandpaths,templist)
  }
  save(allrandpaths,file = paste(as.character(p),"allrandpaths",".Rdata"))
  save(allrandscores,file = paste(as.character(p),"allrandscores",".Rdata"))
  
  #SCORING REAL PATHS
  totalscorelist<-vector()
  pvaluelist<-vector()
  scorematrix<-data.frame(matrix(ncol = 3, nrow = length(patient)))
  colnames(scorematrix)<-c("score","length","p-value")
  for (g in 1:length(patient)) {
    total<-0
    rownames(scorematrix)[g]<-paste(lapply(patient[[g]], function(x) paste(x," ")), collapse = '', sep='')
    scorematrix[g,2]<-length(patient[[g]])
    for (h in 1:length(patient[[g]])) {
      gene<-patient[[g]][h]
      if (gene %in% rownames(fcrnaseq)){ total<-total+abs(const1*fcrnaseq[gene,p])}
      if (gene %in% rownames(context)){
        if(context[unlist(gene),1]==1){total<-total+const2}
        else if(context[unlist(gene),1]==2) {total+const3}
      }
      else {total<-total+const4}
      if(unlist(gene) %in% mutlist & h != 1 & h != length(patient[[g]])) {total<-total+const5}
    }
    total<-total/length(patient[[g]])
    scorematrix[g,1]<-total
    totalscorelist<-c(totalscorelist,total)
    randtotalscorelist<-get(paste("randscore",as.character(length(patient[[g]]))))
    
    # CALCULATING EMPIRICAL P-VALUES
    pvalue<-(length(subset(randtotalscorelist, total <=randtotalscorelist))/randnum)
    pvaluelist<-c(pvaluelist,pvalue)
    scorematrix[g,3]<-pvalue
    patient[[g]]<-c(patient[[g]],pvalue)
  }
  save(pvaluelist,file = paste(as.character(p),"pvalues",".Rdata"))

  # FIND SIGNIFICANT PATHS
  indpath<-list()
  for (i in 1:(length(patient))) {
    if (patient[[i]][[length(patient[[i]])]]<=pvalcut){indpath<-list.append(indpath,patient[[i]][1:(length(patient[[i]])-1)])}
  }
  save(indpath,file=paste(as.character(p),"individualnetworkpaths",".Rdata"))
  # CREATE EDGELIST FROM SIGNIFICANT PATHS 
  indedges<-data.frame()
  for (j in 1:length(indpath)) {
    for(k in 1:(length(indpath[[j]])-1)) {
      indedges<-rbind(indedges,data.frame(indpath[[j]][k],indpath[[j]][k+1]))}
  }
  
  # REMOVING DUPLICATE EDGES 
  indedges<-unique(t(apply(indedges, MARGIN = 1, sort)))
  save(indedges,file=paste(as.character(p),"indedges",".Rdata"))

  indgenes<-unique(c(indedges[,1],indedges[,2]))
  save(indgenes,file=paste(as.character(p),"indnodes",".Rdata"))
  # ADD NETWORK INFORMATION TO PAT. CHARACTERISTICS 
  patientcharacteristics[1,2]<-length(indedges[,1])
  patientcharacteristics[1,3]<-length(indgenes)
  patientcharacteristics[1,4]<-length(indgenes[indgenes %in% mutlist])
  
  # GENERATE INDIVIDUAL DISEASE MODULE 
  indppi<-graph_from_edgelist(indedges,directed = FALSE)
  # BASELINE SCORING OF INDIVIDUALIZED DISEASE NETWORK 
  basecomp<-count_components(indppi)
  patientcharacteristics[1,5]<-basecomp
  
  hubbase<-0
  bottlebase<-0
  thubbase<-0
  tbottlebase<-0
  
  edgelist<-sapply(indgenes, function(x) sum(str_count(indedges,x)))
  betweenbase<-betweenness(indppi,indgenes)
  hubcut<-sort(edgelist)[length(edgelist)-round(length(edgelist)*.05)]
  if(sort(edgelist)[round(length(edgelist)*.05)] > 1){bottlecut<-sort(edgelist)[round(length(edgelist)*.05)]}
  else{bottlecut<-2}
  betweencut<-order(betweenbase)[length(betweenbase)-round(length(betweenbase)*.05)]
  # calculating actual hub/bottleneck numbers 
  for(i in 1:length(indgenes)){
    if(edgelist[i] >= hubcut & betweenbase[i] >= betweencut){thubbase<-thubbase+1}
    else if (edgelist[i] >= hubcut){hubbase<-hubbase+1}
    else if (edgelist[i] <= bottlecut & betweenbase[i] >= betweencut) {tbottlebase<-tbottlebase+1}
    else if (edgelist[i] <=bottlecut){bottlebase<-bottlebase+1}
  }
  patientcharacteristics[1,6]<-hubbase
  patientcharacteristics[1,7]<-thubbase
  patientcharacteristics[1,8]<-bottlebase
  patientcharacteristics[1,9]<-tbottlebase
  
  # INITIALIZE IDGI SCORING STRUCTURE
  indmut<-indgenes[indgenes %in% mutlist]
  inddata <- data.frame(matrix(ncol = 6, nrow = length(indmut)))
  colnames(inddata)<- c( "Mut Network Edges","Mut Network Components","Hubs","Highly Trafficked Hubs", "Bottlenecks",  "Total Score")
  rownames(inddata)<- indmut
  # MUTATION SCORING LIST
  for (j in 1:length(indmut)) {
    gene<- indmut[j]
    mutpath<-list()
    for(x in 1:length(indpath)){
      if( gene != indpath[[x]][1] & gene != indpath[[x]][length(indpath[[x]])]) {
        mutpath<-list.append(mutpath,indpath[[x]])
      }
    }
    # rescore all paths to see if considering gene unmutated makes the path nonsignificant
    remlist<-vector()
    for (z in 1:length(mutpath)) {
      if(gene %in% unlist(mutpath[[z]])) {
        ctotal<-0
        for (w in 1:(length(mutpath[[z]])-1)) {
        cgene<-mutpath[[z]][w]
        if(cgene %in% rownames(fcrnaseq)){ctotal<-ctotal+abs(const1*fcrnaseq[gene,p])}
        if(cgene %in% rownames(context)){
          if(context[unlist(cgene),1]==1){ctotal<-ctotal+const2}
          if(context[unlist(cgene),1]==2){ctotal<-ctotal+const3}
        }
        else {ctotal<-ctotal+const4}
        if(unlist(cgene) %in% mutlist & w!=1 & w!= length(mutpath[[z]]) & cgene != gene){ctotal<-ctotal+const5}
        }
        ctotal<-ctotal/length(patient[[z]])
        randtotalscorelist<-get(paste("randscore",as.character(length(mutpath[[z]]))))
        pvalue<-(length(subset(randtotalscorelist, ctotal <=randtotalscorelist))/randnum)
        if(pvalue >= 0){remlist<-c(remlist,z)}  
      }
    }
    if(length(remlist) != 0){mutpath<-list.remove(mutpath,remlist)}
    
    mutedges<-data.frame()
    for (v in 1:(length(mutpath))) {
      for(t in 1:(length(mutpath[[v]])-1)){mutedges<-rbind(mutedges,data.frame(mutpath[[v]][t],mutpath[[v]][t+1]))}
    }
    mutedges<-as.data.frame(mutedges)
    mutedges<-unique(t(apply(mutedges, MARGIN = 1, sort)))
    mutppi<-graph_from_edgelist(mutedges)
    
    
    
    
    mutgenes<-unlist(get.vertex.attribute(mutppi))
    mutcomp<-count_components(mutppi)
    hub<-0
    thub<-0
    bottle<-0
    tbottle<-0
    
    mutedgelist<-sapply(mutgenes, function(x) sum(str_count(mutedges,x)))
    between<-betweenness(mutppi,mutgenes)
    for (k in 1:length(mutgenes)) {
      if(mutedgelist[k] >= hubcut & between[k] >= betweencut){thub<-thub+1}
      else if (mutedgelist[k] >= hubcut){hub<-hub+1}
      else if (mutedgelist[k] <= bottlecut & between[k] >= betweencut) {tbottle<-tbottle+1}
      else if (mutedgelist[k] <=bottlecut){bottle<-bottle+1}
    }
    inddata[j,1]<-length(mutedges[,1])-length(indedges[,1])
    inddata[j,2]<-mutcomp-basecomp
    inddata[j,3]<-hub-hubbase
    inddata[j,4]<-thub-thubbase
    inddata[j,5]<-bottle-bottlebase
    inddata[j,6]<-(abs(const6*inddata[j,1])+abs(const7*inddata[j,2])+abs(const8*inddata[j,3])+abs(const9*inddata[j,4])+abs(const10*inddata[j,5]))/length(indmut)
  }
  # Saving IDGI scoring
  save(inddata,file=paste(as.character(p),"IDGI",".Rdata"))
  
  #ITGI scoring 
  drugtarg<-ttherapy[ttherapy %in% indgenes]
  drugtarg<-unique(drugtarg)
  
  targscore<-data.frame(matrix(ncol=6,nrow=length(drugtarg)))
  colnames(targscore)<-c( "Targ Network Edges","Targ Network Components","Hubs","Highly Trafficked Hubs", "Bottlenecks",  "Total Score")
  rownames(targscore)<-drugtarg
  
  for (l in 1:length(drugtarg)) {
    targ<-drugtarg[l]
    targppi<-delete.vertices(indppi,targ)
    targedges<-as_edgelist(targppi)
    targgenes<-unlist(get.vertex.attribute(targppi))
    targcomp<-count_components(targppi)
    hubtarg<-0
    thubtarg<-0
    bottletarg<-0
    tbottletarg<-0
    targedgelist<-sapply(targgenes, function(x) sum(str_count(targedges,x)))
    targbetween<-betweenness(targppi,targgenes)
    for (m in 1:length(targgenes)) {
      if(targedgelist[m]>= hubcut & targbetween[m]>= betweencut){thubtarg<-thubtarg+1}
      else if (targedgelist[m]>= hubcut) {hubtarg<-hubtarg+1}
      else if (targedgelist[m]<= bottlecut & targbetween[m]>= betweencut){tbottletarg<-tbottletarg+1}
      else if(targedgelist[m]<= bottlecut){bottletarg<-bottletarg+1}
    }
    targscore[l,1]<-length(targedges[,1])-length(indedges[,1])
    targscore[l,2]<-targcomp-basecomp
    targscore[l,3]<-hubtarg-hubbase
    targscore[l,4]<-thubtarg-thubbase
    targscore[l,5]<-bottletarg-bottlebase
    targscore[l,6]<-(abs(const6*targscore[l,1])+abs(const7*targscore[l,2])+abs(const8*targscore[l,3])+abs(const9*targscore[l,4])+abs(const10*targscore[l,5]))/length(targgenes)
  }
  
  save(targscore,file=paste(as.character(p),"ITGI",".Rdata"))
  
  #Combinatorial ITGI scoring 
  
  # if the patient has more than 20 target scoring genes- check the top 20 ITGI scoring genes
  if (length(targscore[,1]) > 20) {
    combgenes<-rownames(targscore[order(targscore$`Total Score`, decreasing = TRUE),])[1:20]
  }
  else {combgenes<-rownames(targscore)}
  
  drugcomb<-combn(combgenes,2)
  
  combscore<-data.frame(matrix(ncol=6,nrow=length(drugcomb[1,])))
  colnames(combscore)<-c( "Targ Network Edges","Targ Network Components","Hubs","Highly Trafficked Hubs", "Bottlenecks", "Total Score")
  rownames(combscore)<-paste(drugcomb[1,],drugcomb[2,])
  
  for (l in 1:length(drugcomb[1,])) {
    targ<-c(drugcomb[1,l],drugcomb[2,l])
    targppi<-delete.vertices(indppi,targ)
    targedges<-as_edgelist(targppi)
    targgenes<-unlist(get.vertex.attribute(targppi))
    targcomp<-count_components(targppi)
    hubtarg<-0
    thubtarg<-0
    bottletarg<-0
    tbottletarg<-0
    targedgelist<-sapply(targgenes, function(x) sum(str_count(targedges,x)))
    targbetween<-betweenness(targppi,targgenes)
    for (m in 1:length(targgenes)) {
      if(targedgelist[m]>= hubcut & targbetween[m]>= betweencut){thubtarg<-thubtarg+1}
      else if (targedgelist[m]>= hubcut) {hubtarg<-hubtarg+1}
      else if (targedgelist[m]<= bottlecut & targbetween[m]>= betweencut){tbottletarg<-tbottletarg+1}
      else if(targedgelist[m]<= bottlecut){bottletarg<-bottletarg+1}
    }
    combscore[l,1]<-length(targedges[,1])-length(indedges[,1])
    combscore[l,2]<-targcomp-basecomp
    combscore[l,3]<-hubtarg-hubbase
    combscore[l,4]<-thubtarg-thubbase
    combscore[l,5]<-bottletarg-bottlebase
    combscore[l,6]<-(abs(const6*combscore[l,1])+abs(const7*combscore[l,2])+abs(const8*combscore[l,3])+abs(const9*combscore[l,4])+abs(const10*combscore[l,5]))/length(targgenes)
  }
  
  save(combscore,file=paste(as.character(p),"ITGIcombinations",".Rdata"))
}
save(patientcharacteristics,file=paste(as.character(p),"patientcharacteristics",".Rdata"))
  
