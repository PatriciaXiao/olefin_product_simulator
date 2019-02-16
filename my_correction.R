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

for(j in 1:repeats) {
    #
}