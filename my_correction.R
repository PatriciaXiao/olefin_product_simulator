#R version 3.3.2 
  
rxnmix<- function(molecules, nrrxn){
    dimstart<-length(molecules)
    moleculestart<-sample(molecules,1)
    ndoublebond<-moleculestart-1
    doublebond<-sample(c(1:ndoublebond),1)
    moleculecleaved<-c(doublebond,moleculestart-doublebond)
    rest<-sample(moleculecleaved,1)
    ###################### Start of the reaction
    indexstartmolecule<-which(molecules==moleculestart)[1]
    moleculesreaction<-molecules[-indexstartmolecule]
    nmolecules<-dimstart-1
    indexmoleculeselected<-sample(c(1:nmolecules),nrrxn, replace=T)
    i=0
    for(i in 1:nrrxn){
        moleculeselected<-moleculesreaction[indexmoleculeselected[i]]
        ndoublebond<-moleculeselected-1
        doublebond<-sample(c(1:ndoublebond),1)
        resttemp<-sample(c(doublebond,moleculeselected-doublebond),1)
        moleculesreaction[indexmoleculeselected[i]]<-(rest+moleculeselected-resttemp)
        rest<-resttemp
        i+1
    }
    return(moleculesreaction)
}

repeats=10 # 0
x=c(rep(0,110)) # 110 times 0 array


for(j in 1:repeats) {
    molecules12<-c(rep(12,200))
    # print(c(rep(6,100),rep(18,100))) # 100 times 6 and 100 times 18
    mix12<-rxnmix(molecules12, 2000) # 200 or 1000 or 2000 or 20000
    densdat12<-hist(mix12,breaks=c(0:110),plot=FALSE)$density
    # densdat12<-hist(mix12,breaks=c(0:110),plot=TRUE)$density
    x=cbind(x,densdat12)
}

y=rowSums(x)/repeats

# print(y)
