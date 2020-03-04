EDec_stage2_aroundGOI <- function(goi,Tum,propTum,k=2){


require(EDec)

k = 2


groups = colnames(Tum)
groups = as.data.frame(groups)
colnames(groups) = "Sample"
groups$group = sample(1:k,nrow(groups),replace = T)


G = list()

for(i in 1:k){ 
  G[[i]] = Tum[,as.vector(groups[groups$group == i,1])] 
}
rm(i)

D = vector("list",k)



for(i in 1:k){
  P = propTum[as.vector(groups[groups$group == i,1]),]
  D[[i]] = run_edec_stage_2(gene_exp_bulk_samples = G[[i]],cell_type_props = as.matrix(P))
}

rm(i)

groups$newgroup = 0

for(grp in 1:length(G)){
  cur = G[[grp]]
  sam = colnames(cur)
  for(s in sam){
    ms = list()
    mixturesOverCTSLoci = as.matrix(cur[goi,s])
    colnames(mixturesOverCTSLoci) = s
    prop = as.matrix(propTum[s,])
    me = D[[grp]][[1]]
    me = as.matrix(me[goi,])
    colnames(me) = goi
    prop = t(prop)
    rss = norm(mixturesOverCTSLoci - prop%*%me,"F")^2
    ms[[grp]] = rss
    
    oth = seq(1:length(G))
    oth = oth[!oth == grp]
    for(o in oth){
      me = D[[o]][[1]]
      me = as.matrix(me[goi,])
      colnames(me) = goi
      rss = norm(mixturesOverCTSLoci - prop%*%me,"F")^2
      ms[[o]] = rss
      rw = which(grepl(s,groups$Sample))
      groups[rw,"newgroup"] = which.min(ms)
      
    }
  }
}
rm(cur)
rm(me)
rm(mixturesOverCTSLoci)
rm(prop)
rm(grp)
rm(ms)
rm(o)
rm(oth)
rm(rss)
rm(rw)
rm(s)
rm(sam)


num_reas = list()
s_reas= nrow(groups[groups$group != groups$newgroup,])
num_reas[length(num_reas)+1] = s_reas
groups$group = groups$newgroup

it = 0

while(it < 20 & s_reas > 0){

  G = list()
  
  for(i in 1:k){ 
    G[[i]] = Tum[,as.vector(groups[groups$group == i,1])] 
  }
  rm(i)
  
  D = vector("list",k)
  
  
  
  for(i in 1:k){
    P = propTum[as.vector(groups[groups$group == i,1]),]
    D[[i]] = run_edec_stage_2(gene_exp_bulk_samples = G[[i]],cell_type_props = as.matrix(P))
  }
  
  rm(i)
  rm(P)
  
  groups$newgroup = 0
  
  for(grp in 1:length(G)){
    cur = G[[grp]]
    sam = colnames(cur)
    for(s in sam){
      ms = list()
      mixturesOverCTSLoci = as.matrix(cur[goi,s])
      colnames(mixturesOverCTSLoci) = s
      prop = as.matrix(propTum[s,])
      me = D[[grp]][[1]]
      me = as.matrix(me[goi,])
      colnames(me) = goi
      prop = t(prop)
      rss = norm(mixturesOverCTSLoci - prop%*%me,"F")^2
      ms[[grp]] = rss
      
      oth = seq(1:length(G))
      oth = oth[!oth == grp]
      for(o in oth){
        me = D[[o]][[1]]
        me = as.matrix(me[goi,])
        colnames(me) = goi
        rss = norm(mixturesOverCTSLoci - prop%*%me,"F")^2
        ms[[o]] = rss
        rw = which(grepl(s,groups$Sample))
        groups[rw,"newgroup"] = which.min(ms)
        
      }
    }
  }
  rm(cur)
  rm(me)
  rm(mixturesOverCTSLoci)
  rm(prop)
  rm(grp)
  rm(ms)
  rm(o)
  rm(oth)
  rm(rss)
  rm(rw)
  rm(s)
  rm(sam)
  
  s_reas= nrow(groups[groups$group != groups$newgroup,])
  num_reas[length(num_reas)+1] = s_reas
  groups$group = groups$newgroup
  
  it = it+1
}

rm(s_reas)
rm(k)

return(D)

}

