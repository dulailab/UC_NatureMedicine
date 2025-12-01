# Objective: quick example using gbmt to cluster patient trajectories (Adalimumab/Vedolizumab/GEMINI)
# across weeks and compare information criteria for different cluster counts.
# Package "gbmt" is used for trajectory analysis

library(gbmt)
#loading VARSITY data for ada and vedo
pro2_wide <- read.csv("C:\\Users\\ibi9245\\OneDrive - Takeda\\Projects\\VARSITY\\pro2_wide.csv")
ada <- pro2_wide[pro2_wide$Treatment=="Adalimumab",]
vedo <-pro2_wide[pro2_wide$Treatment=="Vedolizumab",]
#define feature set for clustering
varNames <- c("RECBL", "STFRQ")
#generate models of with 3-6 clusters
ada6 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=6, data=ada, scaling=2)
ada5 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=5, data=ada, scaling=2)
ada4 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=4, data=ada, scaling=2)
ada3 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=3, data=ada, scaling=2)
vedo6 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=6, data=vedo, scaling=2)
vedo5 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=5, data=vedo, scaling=2)
vedo4 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=4, data=vedo, scaling=2)
vedo3 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=3, data=vedo, scaling=2)

rbind(ada3$ic,ada4$ic,ada5$ic,ada6$ic )
rbind(vedo3$ic,vedo4$ic,vedo5$ic,vedo6$ic )
#output cluster assignment for a particular model
varsity_vedoRst <- varsity_vedo6$assign.list

#loading gemini dataset
gemini <- read.csv("C:\\Users\\ibi9245\\OneDrive - Takeda\\Projects\\Entyvio\\C13006\\longitudinal PO2 VDZ.csv")
#generate models with 3-6 clusters
gemini_vedo3 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=3, data=gemini, scaling=2)
gemini_vedo4 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=4, data=gemini, scaling=2)
gemini_vedo5 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=5, data=gemini, scaling=2)
gemini_vedo6 <- gbmt(x.names=varNames, unit="SUBJID", time="Week", d=2, ng=6, data=gemini, scaling=2)

rbind(gemini_vedo3$ic,gemini_vedo4$ic,gemini_vedo5$ic,gemini_vedo6$ic )

#output cluster assignment
gemini_vedoRst <- gemini_vedo4$assign.list


