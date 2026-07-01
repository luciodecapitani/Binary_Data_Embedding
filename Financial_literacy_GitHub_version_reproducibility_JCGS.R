#code for reproducing the results on Financial Literacy in paper
#" Bidimensional Embedding and Visualization of Multidimensional Ordinal Binary Data"


#install.packages("poseticDataAnalysis")
library(poseticDataAnalysis)
source(file="preliminary_function_JCGS.R")

#loading dataset
DT<-read.csv("iacofi.csv")

#selecting the variables of interest
SUBDT<-DT[,c("qk4", "qk5", "qk6", "qk7_1","qk7_1alt", "qk7_2", "qk7_3", "qk7_3alt","qd1","qd7","area5","reddito","edu","capofam")]

#Union of the answers to questions qk7_1 e qk7_3: they are substantially the same question
SUBDT[which(is.na(SUBDT[,7])),7]<-SUBDT[which(is.na(SUBDT[,7])),8]
SUBDT[which(is.na(SUBDT[,4])),4]<-SUBDT[which(is.na(SUBDT[,4])),5]
#removing qk7_1alt and qk7_3alt from the dataset after merging the answers
SUBDT$qk7_1alt<-SUBDT$qk7_3alt<-NULL
#specifying the coding of correct answers
#qk4:0; qk5: 102; qk6: 1; qk71(merged): 1;  qk72: 1; qk73 (merged): 1
corrette<-c(0,102,1,1,1,1)

#vector of sample weights
weights<-DT$wght
weights<-weights/sum(weights)
#sample size
N<-nrow(SUBDT)


MIS<-matrix(-2, nrow=N, ncol=6) #creating the dataframe with binary responses (-2 is a fictitious value to initialize the data matrix)
MIS<-as.data.frame(MIS)
colnames(MIS)<-c("qk4", "qk5", "qk6", "qk7_1", "qk7_2", "qk7_3")
for(h in 1:nrow(SUBDT)){
  MIS[h,]<-(SUBDT[h,1:6]==corrette)
}

############################################################################
############################################################################
#building the frequency distribution with and without classification for
#educational degree
############################################################################
############################################################################

#without classification
B<-profiles_freq_dist_weighted(MIS,weights)

#with classification: educational degree
edu_labels<-c("Bachelor/master or higher ","High school", "Middle school", "primary school", "No degree")
edu<-vector(length=N)
edu[which((SUBDT$edu==1))]<-edu_labels[1]
edu[which((SUBDT$edu>=2) & (SUBDT$edu<=3))]<-edu_labels[2]
edu[which((SUBDT$edu>=4) & (SUBDT$edu<=5))]<-edu_labels[3]
edu[which((SUBDT$eduZ>=6) & (SUBDT$edu<=7))]<-edu_labels[4]
edu[which((SUBDT$edu==8))]<-edu_labels[5]
B_edu<-profiles_freq_dist_weighted_classified(MIS,weights,edu,edu_labels)

############################################################################
############################################################################
#embedding analisys with no classification
############################################################################
############################################################################

obs_Profiles<-B[[1]]
freq_W<-B[[3]]
freq<-B[[2]]/sum(B[[2]])

#computing optimal embedding without classificaiton
result <- OptimalBidimensionalEmbedding(obs_Profiles, freq_W, output_every_sec=1)

############################################################################
############################################################################
#Description of embedding results
############################################################################
############################################################################

#Figure 2 top left panel
x11()
plot(sort(result$allLoss), main ="Global loss", xlab="2-D map index", ylab="Loss", ylim=c(0,max(result$allLoss)), pch=20)

#Figure 2 top right panel
x11()
plot(sort(result$bestRepresentation[,5]), main ="Profile's approximation error", xlab="Unit index", ylab="Approximation error", ylim=c(0,max(result$bestRepresentation[,5])), pch=20)

#Figure 2 bottom left panel
npoints<-ceiling(N*result$bestRepresentation[,4])
errorep<-rep(result$bestRepresentation[,5],npoints)
x11()
boxplot(errorep, boxwex=.2, range=0, cex.axis=1,  cex.lab=1, col="white", ylab="Approximation error", main="Units' approximation error distribution")#cex.main=1.8,

#Figure 2 bottom right panel
npoints<-ceiling(N*result$bestRepresentation[,4])
freq<-npoints/sum(npoints)
x11()
plot(freq, result$bestRepresentation[,5], main ="Profile's pproximation error vs Population share", xlab="Population share", ylab="Approximation error", ylim=c(0,max(result$bestRepresentation[,5])), pch=20)


############################################################################
############################################################################
#Embedding graphical representation
############################################################################
############################################################################

#Figure 3 right panel
x11()
POsetDR.plot(result, N=2304,k=6, gV_x=3, gV_y=3,Grid = T, Grid_lty = c(1,5,2,4,3)[1:3],Grid_lwd = c(2,2,2,1,1)[1:3],Grid_col = "darkgray", 
Profile = T, Profiles_col = "black", varNames = c(expression(K[1]),expression(K[2]),expression(K[3]),expression(K[4]),expression(K[5]),expression(K[6])), 
Profiles_expand=c(0.5,0.5),title = " ", type= "bubbles_and_approximation",expandX = c(0,0), expandY = c(0.1,0.1), points_col="black", 
center_col= "black",bubbles_col = "gray92",center_pch = 16, points_cex = 0.3,center_cex = 0.4,jitter_amount = 3,expand_radius = 30,cex.profile=1)

#Figure 3 left panel
X11()
POsetDR.plot(result, N=3000,k=6,gV_x=3,gV_y=3,Grid = TRUE,Grid_lty = c(1,5,2,4,3)[1:3], Grid_lwd = c(3,2.5,2,1.5,1)[1:3], Grid_col = "darkgray", 
Profile = TRUE, Profiles_col = "black", varNames = c(expression(K[1]),expression(K[2]),expression(K[3]),expression(K[4]),expression(K[5]),expression(K[6])),
Profiles_expand=c(0.5,0.5),title = " ", type="jitter_circle",expandX = c(0,0), expandY = c(0,0), points_col="black",center_col= "red", 
points_pch = 16,center_pch = 16, points_cex = 0.4,center_cex = 0.7, jitter_amount = 2,rampPalette = "viridis",bubble_sizes = NULL,cex.profile=1)

############################################################################
############################################################################
#embedding analisys with classification based on educational level
############################################################################
############################################################################
result_edu1<-
  result_edu2<-
  result_edu3<-
  result_edu4<-
  result_edu5<-result

result_edu1$bestRepresentation$weights<-B_edu[[5]][,1]
result_edu2$bestRepresentation$weights<-B_edu[[5]][,2]
result_edu3$bestRepresentation$weights<-B_edu[[5]][,3]
result_edu4$bestRepresentation$weights<-B_edu[[5]][,4]
result_edu5$bestRepresentation$weights<-B_edu[[5]][,5]

#edu 1: "Bachelor/master or higher "
#Figure 5 top right panel
x11()
POsetDR.plot(result_edu1, N=2304, k=6, gV_x=3,gV_y=3, Grid = T, Grid_col = "darkgray",Grid_lty = c(1,5,2,4,3)[1:3], Grid_lwd = c(2,2,2,1,1)[1:3], 
             Profile = F, Profiles_col = "black",
             varNames = c(expression(K[1]),expression(K[2]),expression(K[3]),expression(K[4]),expression(K[5]),expression(K[6])),
             Profiles_expand=c(0.5,0.5), title = "Bachelor/master or higher ", type= "bubbles", expandY = c(3,3), expandX = c(3,3),
             points_col="black", center_col= "black",bubbles_col = "gray92", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 48,cex.profile=1)
#edu 2:High school"
#Figure 5 top left panel
x11()
POsetDR.plot(result_edu2, N=2304, k=6,  gV_x=3,gV_y=3,  Grid = T, Grid_col = "darkgray", Grid_lty = c(1,5,2,4,3)[1:3], Grid_lwd = c(2,2,2,1,1)[1:3], 
             Profile = F, Profiles_col = "black",
             varNames = c(expression(K[1]),expression(K[2]),expression(K[3]),expression(K[4]),expression(K[5]),expression(K[6])),
             Profiles_expand=c(0.5,0.5), title = "High school", type= "bubbles",  expandY = c(3,3), expandX = c(3,3),
             points_col="black", center_col= "black",bubbles_col = "gray92", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 34,cex.profile=1)
#edu 3: Middle school
#Figure 5 bottom right panel
x11()
POsetDR.plot(result_edu3, N=2304, k=6,  gV_x=3,gV_y=3,  Grid = T, Grid_col = "darkgray", Grid_lty = c(1,5,2,4,3)[1:3], Grid_lwd = c(2,2,2,1,1)[1:3], 
             Profile = F, Profiles_col = "black",
             varNames = c(expression(K[1]),expression(K[2]),expression(K[3]),expression(K[4]),expression(K[5]),expression(K[6])),
             Profiles_expand=c(0.5,0.5), title = "Middle school", type= "bubbles",  expandY = c(3,3), expandX = c(3,3),
             points_col="black", center_col= "black",bubbles_col = "gray92", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 55,cex.profile=1)
#edu 4: no statistical unit with this education level
#x11()
#POsetDR.plot(result_edu4, N=2304, k=6, gV=3, Grid = T, Grid_col = "darkgray", Profile = T, Profiles_col = "black",
#             varNames = c(expression(K[1]),expression(K[2]),expression(K[3]),expression(K[4]),expression(K[5]),expression(K[6])),
#             Profiles_expand=0.5, title = "elementari", type= "bubbles", expandY = c(0,0), expandX = c(0,1.1),
#             points_col="black", center_col= "black",bubbles_col = "gray92", center_pch = 16, points_cex = 0.3,
#             center_cex = 0.4, jitter_amount = 3, expand_radius = 1.3)
#edu 5: no degree
#Figure 5 bottom left panel
x11()
POsetDR.plot(result_edu5, N=2304, k=6, gV_x=3,gV_y=3,  Grid = T, Grid_col = "darkgray", Grid_lty = c(1,5,2,4,3)[1:3], Grid_lwd = c(2,2,2,1,1)[1:3], 
             Profile = F, Profiles_col = "black",
             varNames = c(expression(K[1]),expression(K[2]),expression(K[3]),expression(K[4]),expression(K[5]),expression(K[6])),
             Profiles_expand=c(0.5,0.5), title = "No degree", type= "bubbles", expandY = c(3,3), expandX = c(3,3),
             points_col="black", center_col= "black",bubbles_col = "gray92", center_pch = 16, points_cex = 0.3,
             center_cex = 0.4, jitter_amount = 3, expand_radius = 135,cex.profile=1)

#########################################
##Contingency analisys
#########################################
#1) marginal frequencies definition
freqM<-matrix(0,6,4)
colnames(freqM)<-c("0_count","1_count","0_weigth","1_weigth")
rownames(freqM)<-c("K_1","K_2","K_3","K_4","K_5","K_6")
for(i in 1:6){
  freqM[i,1]<-sum((B[[1]][,i]==0)*B[[2]])
  freqM[i,2]<-sum((B[[1]][,i]==1)*B[[2]])
  freqM[i,3]<-sum((B[[1]][,i]==0)*B[[3]]) #pesi dello 0
  freqM[i,4]<-sum((B[[1]][,i]==1)*B[[3]]) #pesi dell'1
}
#2) complete profiles list
n_var<-6
complete_matrix_profile<-matrix(0,2^n_var,n_var+5)
for(i in 1:n_var){
  complete_matrix_profile[,i]<-rep(c(0,1), each=2^(n_var-i), times=2^(i-1))
}
colnames(complete_matrix_profile)<-c("K_1","K_2","K_3","K_4","K_5","K_6","indep","observed","relative_contingencies","x","y")
#3)theoretical and observed frequency computation 
pow<-c(2^3,2^5,2^4,2^2,1,2)
pow_r<-c(2^2,1,2,2^3,2^5,2^4)
for(i in 1:(2^n_var)){
  complete_matrix_profile[i,7]<-prod((complete_matrix_profile[i,1:6]==0)*freqM[,3]+(complete_matrix_profile[i,1:6]==1)*freqM[,4])
  complete_matrix_profile[i,10]<-sum(complete_matrix_profile[i,1:6]*pow)+1
  complete_matrix_profile[i,11]<-sum(complete_matrix_profile[i,1:6]*pow_r)+1
}
complete_matrix_profile[result$bestRepresentation$profiles+1,8]<-result$bestRepresentation$weights
#4) relative contingencies computation
complete_matrix_profile[,9]<-(complete_matrix_profile[,8]-complete_matrix_profile[,7])/complete_matrix_profile[,7]



variable_priority<-c(2,3,1,4,6,5)
k<-6

#Figure 4
Contingency.plot(complete_matrix_profile, variable_priority, k=6, gV=3,Grid = TRUE,
                       varNames=c("K1","K2","K3","K4","K5","K6"),
                       Grid_col = "black", Grid_lty = c(1,5,2,4,3)[1:3], Grid_lwd = c(3,2.5,2,1.5,1)[1:3],
                       Profile = TRUE, Profiles_col = "black",
                       Profiles_expand=1, title = NULL, expandX = c(0,0),
                       expandY = c(0,0), rampPalette1 = c("white", "darkred"),
                       rampPalette2 = c("white", "darkgreen"))


