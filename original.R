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

##############
# Simulation with 200 dodecene molecules (C12), 200/1,000/2,000/20,000 steps, average of 100 runs (Figure 11 in the main text)


repeats=100
x=c(rep(0,110))
for(j in 1:repeats) {
    source("D:/simulation_rxnmix.R")
    molecules12<-c(rep(12,200))
    mix12<-rxnmix(molecules12,200 or 1000 or 2000 or 20000)
    densdat12<-hist(mix12,breaks=c(0:110),plot=FALSE)$density
    x=cbond(x,densdat12)
    j+1
}

y=rowSums(x)/repeats
write.table(y,file="D:/C12_200_1000_100.txt",sep="\t",dec =",")

# Simulation with 200 dodecene molecules (C12), 2,000 steps, average of 100 runs (Figures 12 and 13 in the main text)

repeats=100
x=c(rep(0,110))
for(j in 1:repeats) {
    source("D:/simulation_rxnmix.R")
    molecules12<-c(rep(12,200))
    mix12<-rxnmix(molecules12,2000)
    densdat12<-hist(mix12,breaks=c(0:110),plot=FALSE)$density
    x=cbond(x,densdat12)
    j+1
}
y=rowSums(x)/repeats
write.table(y,file="D:/C12_200_1000_100.txt",sep="\t",dec =",")

# Simulation with 100 hexene and 100 octadecene molecules (C6/C18 1:1), 2,000 steps, average of 100 runs (Figure 12 in the main text)

repeats=100
x=c(rep(0,110))
for(j in 1:repeats) {
    source("D:/simulation_rxnmix.R")
    molecules618<-c(rep(6,100),rep(18,100))
    mix618<-rxnmix(molecules618,2000)
    densdat618<-hist(mix618,breaks=c(0:110),plot=FALSE)$density
    x=cbond(x,densdat618)
    j+1
}
y=rowSums(x)/repeats
write.table(y,file="D:/C12_200_1000_100.txt",sep="\t",dec =",")


# Simulation with 200 propene molecules (C3), 2,000 steps, average of 100 runs (Figure 13 in the main text)
repeats=100
x=c(rep(0,110))
for(j in 1:repeats) {
    source("D:/simulation_rxnmix.R")
    molecules3<-c(rep(3,200))
    mix3<-rxnmix(molecules3,2000)
    densdat3<-hist(mix3,breaks=c(0:110),plot=FALSE)$density
    x=cbond(x,densdat3)
    j+1
}
y=rowSums(x)/repeats
write.table(y,file="D:/C12_200_1000_100.txt",sep="\t",dec =",")

# Simulation with 200 hexene molecules (C6), 2,000 steps, average of 100 runs (Figure 13 in the main text)
repeats=100
x=c(rep(0,110))
for(j in 1:repeats) {
    source("D:/simulation_rxnmix.R")
    molecules6<-c(rep(6,200))
    mix6<-rxnmix(molecules6,2000)
    densdat6<-hist(mix6,breaks=c(0:110),plot=FALSE)$density
    x=cbond(x,densdat6)
    j+1
}
y=rowSums(x)/repeats
write.table(y,file="D:/C12_200_1000_100.txt",sep="\t",dec =",")
