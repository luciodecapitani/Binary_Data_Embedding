############################################################################
############################################################################
# Preliminary_Functions_Financial_literacy.r
#
# Packages and functions used by Financial_literacy_JCGS.r, which reproduces the results of the
# Financial knowledge example of the paper
# "Bidimensional Embedding and Visualization of Multidimensional Ordinal Binary Data".
# This file is sourced at the beginning of Financial_literacy_JCGS.r.
#
# Notation. The data are N statistical units scored on k ordinal binary variables. Each unit is
# described by its binary profile (a 0/1 vector of length k); a profile is identified by its decimal
# value, the first variable being the most significant digit (e.g. with k = 6, profile 000101 has
# identifier 5). The optimal bidimensional embedding (OBE) is computed by the function
# OptimalBidimensionalEmbedding of the R package poseticDataAnalysis.
#
# STRUCTURE OF THE FILE
#   A. Installation (only if needed) and loading of the packages
#   B. Frequency distributions over the binary profiles
#   C. Order structure of the input lattice and of the optimal embedding
#   D. Relative contingencies (Figure 4 of the paper)
#   E. Plots: OBE map (Figures 3 and 5) and contingency map (Figure 4)
############################################################################
############################################################################


############################################################################
# A. Installation (only if needed) and loading of the packages
############################################################################
# - poseticDataAnalysis (CRAN): Optimal Bidimensional Embedding and poset tools;
# - viridis (CRAN): color scale of the approximation errors in the OBE map.
# Nothing is installed if the packages are already available.

for (pkg in c("poseticDataAnalysis", "viridis")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(poseticDataAnalysis)
library(viridis)


############################################################################
# B. Frequency distributions over the binary profiles
############################################################################

# Frequency distribution over the OBSERVED binary profiles.
#   matrix  : data frame/matrix of the binary responses (one row per unit, one column per variable)
#   weights : vector of the sample weights of the units
# Output: a list with
#   [[1]] matrix of the observed profiles (one row per profile, sorted by their decimal identifier);
#   [[2]] sample counts of the observed profiles;
#   [[3]] frequencies of the observed profiles obtained by sample weights (sum of the weights of the
#         units having that profile; they sum to 1 if the weights do).
profiles_freq_dist_weighted <- function(matrix, weights){
  k <- ncol(matrix)
  #decimal identifier of the profile of each unit
  dec_profiles <- as.matrix(matrix) %*% (2^((k:1)-1))
  observed_profiles_list <- sort(unique(dec_profiles))
  freq <- freq_W <- vector(length = length(observed_profiles_list))
  for(i in 1:length(observed_profiles_list)){
    freq_W[i] <- sum(weights*(dec_profiles == observed_profiles_list[i]))
    freq[i] <- sum(dec_profiles == observed_profiles_list[i])
  }
  #binary profiles corresponding to the observed identifiers
  obs_Profiles <- t(sapply(observed_profiles_list, function(x){ as.integer(intToBits(x)) })[k:1,])
  return(list(obs_Profiles, freq, freq_W))
}

# Frequency distribution over the observed binary profiles, split according to the classes of a
# grouping variable (e.g. the education degree, Figure 5 of the paper).
#   matrix, weights    : as in profiles_freq_dist_weighted
#   grouping_variable  : vector with the class of each unit
#   grouping_lables    : labels of the classes (if NULL, the distinct values of grouping_variable)
# Output: a list with
#   [[1]], [[2]], [[3]] as in profiles_freq_dist_weighted (whole sample);
#   [[4]] matrix of the sample counts: one row per observed profile, one column per class;
#   [[5]] matrix of the frequencies obtained by sample weights: one row per observed profile, one column
#         per class (the frequency of a profile is split among the classes of the grouping variable).
profiles_freq_dist_weighted_classified <- function(matrix, weights, grouping_variable, grouping_lables = NULL){
  if(is.null(grouping_lables)){
    grouping_lables <- unique(grouping_variable)
  }
  k <- ncol(matrix)
  dec_profiles <- as.matrix(matrix) %*% (2^((k:1)-1))
  observed_profiles_list <- sort(unique(dec_profiles))
  freq <- freq_W <- vector(length = length(observed_profiles_list))
  group_freq <- group_freq_W <- matrix(0, length(observed_profiles_list), length(grouping_lables))
  colnames(group_freq) <- colnames(group_freq_W) <- grouping_lables
  for(i in 1:length(observed_profiles_list)){
    freq_W[i] <- sum(weights*(dec_profiles == observed_profiles_list[i]))
    freq[i] <- sum(dec_profiles == observed_profiles_list[i])
    for(j in 1:length(grouping_lables)){
      group_freq[i,j] <- sum((dec_profiles == observed_profiles_list[i])*(grouping_variable == grouping_lables[j]))
      group_freq_W[i,j] <- sum(weights*(dec_profiles == observed_profiles_list[i])*(grouping_variable == grouping_lables[j]))
    }
  }
  obs_Profiles <- t(sapply(observed_profiles_list, function(x){ as.integer(intToBits(x)) })[k:1,])
  return(list(obs_Profiles, freq, freq_W, group_freq, group_freq_W))
}

# Frequency distribution over ALL of the 2^k binary profiles (unobserved profiles have frequency 0).
#   matrix, weights : as in profiles_freq_dist_weighted
# Output: a list with [[1]] matrix of the 2^k profiles (row i is the profile with identifier i-1) and
#   [[2]] their frequencies obtained by sample weights.
complete_profiles_freq_dist_weighted <- function(matrix, weights){
  k <- ncol(matrix)
  profiles_list <- 0:(2^k-1)
  dec_profiles <- as.matrix(matrix) %*% (2^((k:1)-1))
  freq_W <- vector(length = length(profiles_list))
  for(i in 0:(2^k-1)){
    freq_W[i+1] <- sum(weights*(dec_profiles == i))
  }
  obs_Profiles <- t(sapply(profiles_list, function(x){ as.integer(intToBits(x)) })[k:1,])
  return(list(obs_Profiles, freq_W))
}


############################################################################
# C. Order structure of the input lattice and of the optimal embedding
############################################################################

# Coordinates of ALL of the 2^k binary profiles in the bidimensional embedding generated by a variable
# index permutation (e.g. the optimal one, bestVariablesPriority in the output of
# OptimalBidimensionalEmbedding). The x coordinate is the height of the profile in the lexicographic
# linear extension in which the variables have priority perm[1], perm[2], ..., perm[k]; the y coordinate
# is its height in the lexicographic linear extension with the reversed priority perm[k], ..., perm[1].
#   profiles : matrix of binary profiles (one row per profile)
#   perm     : variable index permutation
# Output: a matrix with columns x and y (one row per profile), with values in 1,...,2^k.
embedding_coordinates <- function(profiles, perm){
  k <- length(perm)
  pow <- pow_r <- numeric(k)
  pow[perm] <- 2^((k-1):0)          #weight of each variable in the x coordinate
  pow_r[rev(perm)] <- 2^((k-1):0)   #weight of each variable in the y coordinate
  cbind(x = as.vector(as.matrix(profiles) %*% pow) + 1, y = as.vector(as.matrix(profiles) %*% pow_r) + 1)
}

# Numbers of comparable and incomparable pairs of distinct elements of a poset, given its elements as
# points of R^d ordered component-wise (profiles of the Boolean lattice, or coordinates of the
# bidimensional embedding): two elements are comparable if one is component-wise lower than or equal to
# the other.
#   points : matrix with one row per element
# Output: a vector with the numbers of comparable and incomparable pairs.
count_comparabilities <- function(points){
  points <- as.matrix(points)
  n <- nrow(points)
  #le[i, j] = TRUE if element i is component-wise lower than or equal to element j
  le <- matrix(FALSE, n, n)
  for(i in 1:n) le[i, ] <- apply(points, 1, function(p) all(points[i, ] <= p))
  comp <- (le | t(le))[upper.tri(le)]
  c(comparable = sum(comp), incomparable = sum(!comp))
}


############################################################################
# D. Relative contingencies (Figure 4 of the paper)
############################################################################

# Matrix of the 2^k binary profiles with their observed frequency, the frequency expected under
# independence of the k variables, the relative contingency and the cell of the profile in the OBE map.
#   matrix          : data frame/matrix of the binary responses (one row per unit)
#   weights         : sample weights of the units (internally rescaled to sum 1)
#   variablepriority: variable index permutation of the embedding (e.g. bestVariablesPriority)
#   varnames        : names of the k variables (used as column names)
# Output: a 2^k x (k+5) matrix (row i = profile with identifier i-1) with columns
#   1..k     : the binary profile;
#   indep    : frequency expected under independence (product of the marginal frequencies of its scores);
#   observed : observed frequency (obtained by sample weights);
#   relative_contingencies : (observed - indep)/indep;
#   x, y     : coordinates of the profile in the embedding (as in embedding_coordinates).
complete_matrix_profiles_function <- function(matrix, weights, variablepriority, varnames){
  weights <- weights/(sum(weights))
  k <- ncol(matrix)
  MIS_w <- complete_profiles_freq_dist_weighted(matrix, weights)
  #1) marginal frequencies of the scores 0 and 1 of each variable
  freqM <- matrix(0, k, 2)
  colnames(freqM) <- c("0_weigth", "1_weigth")
  rownames(freqM) <- varnames
  for(i in 1:k){
    freqM[i,1] <- sum((MIS_w[[1]][,i] == 0)*MIS_w[[2]])
    freqM[i,2] <- sum((MIS_w[[1]][,i] == 1)*MIS_w[[2]])
  }
  #2) complete list of profiles, with their observed frequencies
  complete_matrix_profile <- matrix(0, 2^k, k+5)
  complete_matrix_profile[,1:k] <- MIS_w[[1]]
  complete_matrix_profile[,k+2] <- MIS_w[[2]]
  colnames(complete_matrix_profile) <- c(varnames, "indep", "observed", "relative_contingencies", "x", "y")
  #3) frequencies expected under independence and coordinates in the embedding
  for(i in 1:(2^k)){
    complete_matrix_profile[i,k+1] <- prod((complete_matrix_profile[i,1:k] == 0)*freqM[,1] + (complete_matrix_profile[i,1:k] == 1)*freqM[,2])
  }
  complete_matrix_profile[, c(k+4, k+5)] <- embedding_coordinates(complete_matrix_profile[,1:k], variablepriority)
  #4) relative contingencies
  complete_matrix_profile[,k+3] <- (complete_matrix_profile[,k+2] - complete_matrix_profile[,k+1])/complete_matrix_profile[,k+1]
  return(complete_matrix_profile)
}


############################################################################
# E. Plots
############################################################################

# OBE map: draws the optimal bidimensional embedding in the 2^k x 2^k grid, tessellated according to the
# scores of the variables (in the order given by the optimal variable index permutation).
# The units/profiles can be represented in different ways (argument type):
#   "point"         : one point per profile;
#   "jitter_square" / "jitter_circle": one point per unit, jittered around the profile position
#                     (in a square or in a circle); the profile is marked by a point at the center;
#   "bubbles"       : one circle per profile, with area proportional to its frequency;
#   "bubbles_and_approximation": as "bubbles", with the circles colored according to the profile
#                     approximation error (formula (2) of the paper), and a color legend on the right.
POsetDR.plot<-function(result, #output of OptimalBidimensionalEmbedding
                       N, #number of points drawn in the jittered clouds ("jitter_square", "jitter_circle"):
                       #profile i is represented by round(N*frequency_i) points
                       k,#number of binary variables
                       gV_x=ceiling(k/2),#number of binary variables considered in building the grid on the x-axis
                       gV_y=k-ceiling(k/2),#number of binary variables considered in building the grid on the y-axis
                       Grid = TRUE, #if TRUE the grid of the tessellation is drawn
                       Grid_lty = c(1,5,2,4,3)[1:gV_x], #line type of the grid segments associated to each variable
                       Grid_lwd = c(3,2.5,2,1.5,1)[1:gV_x], #line width of the grid segments associated to each variable
                       Grid_col = "blue", #color of the grid
                       Profile = TRUE, #if TRUE the variable scores are reported on both axes (profile legend)
                       Profiles_col = "black", #color of the profile legend
                       varNames = NULL, #string/expression vector with the variables' names (in the original order)
                       Profiles_expand=c(1,1), #spacing of the rows of the profile legend (x axis, y axis):
                       #if greater (lower) than one, the space between the rows is expanded (reduced)
                       title = "Optimal profile map", #main title of the plot
                       type=c("point", "jitter_square", "jitter_circle", "bubbles", "bubbles_and_approximation"), #type of plot
                       expandX = c(0,0), #extra space on the x-axis: c[1] = added to the left, c[2] = added to the right
                       expandY = c(0,0), #extra space on the y-axis: c[1] = added to the bottom, c[2] = added to the top
                       expand_radius=1, #multiplicative factor of the radius of the bubbles
                       points_col="black", #color of the points ("point", "jitter_square", "jitter_circle") or of the
                       #border of the bubbles ("bubbles", "bubbles_and_approximation")
                       center_col= "red", #color of the points drawn at the center of the clouds of points/bubbles
                       optimize_radius=FALSE, #if TRUE the radii are rescaled by the minimum distance between profiles
                       bubbles_col = "gray", #color of the bubbles (type = "bubbles")
                       points_pch = 16, #symbol of the points ("point", "jitter_square", "jitter_circle")
                       center_pch = 16, #symbol of the points drawn at the center of the clouds of points/bubbles
                       points_cex = 0.5, #size of the points ("point", "jitter_square", "jitter_circle")
                       center_cex = 0.5, #size of the points drawn at the center of the clouds of points/bubbles
                       jitter_amount = 1, #dispersion of the clouds of points around the profiles ("jitter_square", "jitter_circle")
                       rampPalette = "viridis", #color scale of the errors: "grey"/"gray" or a viridis option
                       bubble_sizes = NULL, #not used
                       cex.profile=1, #size of the text of the profile legend (variable scores and names on the axes)
                       cex.legend=cex.profile, #size of the numbers of the color legend ("bubbles_and_approximation")
                       min_radius=0, #minimum radius of the bubbles
                       lwd_circle=1, #line width of the border of the bubbles ("bubbles_and_approximation")
                       error_max=NULL, #upper bound of the color scale ("bubbles_and_approximation"; if NULL, max of the errors)
                       cex.main=par("cex.main"), #size of the title
                       asp=1, #aspect ratio of the plot: with asp = 1 the grid is drawn as a square, whatever the
                       #size of the graphical device (the plot area is enlarged in one direction if needed)
                       xlim=NULL, #optional limits of the plot area (they override the default ones, computed from
                       ylim=NULL, #expandX and expandY); e.g. xlim = ylim = c(0, 2^k) makes the plot area coincide with
                       #the grid, the profile legend and the color legend being drawn in the margins (with par(xpd = NA))
                       Profiles_lines=NULL){ #length (in grid units) of the dotted lines separating the columns/rows of
                       #the profile legend; if NULL, 2^(k-gV)*(gV+1)*Profiles_expand

  #names of the components of the output of OptimalBidimensionalEmbedding (older versions returned 3 components)
  if(length(result)==3){
    names(result)[c(2,3)]<-c("bestVariablesPriority","bestRepresentation")}
  #default variables' names: V[i] is the variable in the i-th column of the original dataset
  if(is.null(varNames)){
    varNames<-vector()
    for(i in 1:k){
      varNames<-append(varNames, parse(text=(paste0("V[",i,"]"))))
    }
  }
  #variables' names permuted according to the optimal variable priority
  varNames<-varNames[result$bestVariablesPriority]

  #DIMENSIONS OF THE PLOT AREA
  #lower bounds: they depend on the presence of the profile legend
  ylow<-(-2^(k-gV_y)*gV_y)*(Profile == TRUE)+(Profile != TRUE)-expandY[1]
  xlow<-(-2^(k-gV_x)*gV_x)*(Profile == TRUE)+(Profile != TRUE)-expandX[1]
  #upper bounds: they depend on the presence of the profile legend and, for x, of the color legend
  yup<- (2^k+2^(k-gV_y)*3/2)*(Profile == TRUE)+2^k*(Profile != TRUE)+expandY[2]
  xup<- 2^k*(type != "bubbles_and_approximation")+(2^k+2^(k-gV_x)*3)*(type == "bubbles_and_approximation")+expandX[2]
  if(!is.null(xlim)){ xlow <- xlim[1]; xup <- xlim[2] }
  if(!is.null(ylim)){ ylow <- ylim[1]; yup <- ylim[2] }

  #EMPTY PLOT AREA
  plot(x=NULL,xlim=c(xlow,xup), ylim=c(ylow,yup),xlab=" ", ylab=" ", main=title, yaxt="n",xaxt="n",bty="n",
       cex.main=cex.main, asp=asp)

  #POSITION OF THE PROFILES (columns 2 and 3 of bestRepresentation; the profiles on the border of the grid
  #are moved inwards by half a unit)
  xcoord<-result$bestRepresentation[,2]+0.5*(result$bestRepresentation[,2]==1)-0.5*(result$bestRepresentation[,2]==2^k)
  ycoord<-result$bestRepresentation[,3]+0.5*(result$bestRepresentation[,3]==1)-0.5*(result$bestRepresentation[,3]==2^k)
  if(type=="point"){
    points(xcoord, ycoord, pch=points_pch, cex=points_cex, col=points_col)
  }
  if(type=="jitter_square"){
    #number of points of each profile, proportional to its frequency (column 4 of bestRepresentation)
    npoints<-round(N*result$bestRepresentation[,4]/sum(result$bestRepresentation[,4]))
    x_unit<-rep(xcoord,npoints)
    y_unit<-rep(ycoord,npoints)
    #uniform jittering in a square (profiles with a single point are not jittered)
    Coordin1<-jitter(x_unit,amount=jitter_amount)*(rep(npoints,npoints)>1)+x_unit*(rep(npoints,npoints)==1)
    Coordin2<-jitter(y_unit,amount=jitter_amount)*(rep(npoints,npoints)>1)+y_unit*(rep(npoints,npoints)==1)
    points(Coordin1, Coordin2, pch=points_pch,  cex=points_cex, col=points_col)
    points(xcoord, ycoord, pch=center_pch,  cex=center_cex, col=center_col)
  }
  if(type=="jitter_circle"){
    npoints<-round(N*result$bestRepresentation[,4]/sum(result$bestRepresentation[,4]))
    x_unit<-rep(xcoord,npoints)
    y_unit<-rep(ycoord,npoints)
    #random jittering in a circle of radius jitter_amount (profiles with a single point are not jittered)
    angolo<-runif(sum(npoints),max=2*pi)
    radius<-runif(sum(npoints), max=jitter_amount)
    Coordin1<-(x_unit+cos(angolo)*radius)*(rep(npoints,npoints)>1)+x_unit*(rep(npoints,npoints)==1)
    Coordin2<-(y_unit+sin(angolo)*radius)*(rep(npoints,npoints)>1)+y_unit*(rep(npoints,npoints)==1)
    points(Coordin1, Coordin2, pch=points_pch,  cex=points_cex, col=points_col)
    points(xcoord, ycoord, pch=center_pch,  cex=center_cex, col=center_col)
  }
  if(type=="bubbles"){
    if(optimize_radius==TRUE){
      DIStz<-(outer(xcoord,xcoord, function(x,y) (x-y)^2)+outer(ycoord,ycoord, function(x,y) (x-y)^2))^0.5
      diag(DIStz)<-max(DIStz)
      rad<-min(DIStz)}
    else{
      rad<-1
    }
    #radius such that the area of the bubble is proportional to the frequency
    raggi<-(result$bestRepresentation[,4]/pi)^0.5
    raggi<-raggi*rad*expand_radius
    raggi[which(raggi<=min_radius)]<-min_radius
    symbols(xcoord, ycoord, circles = raggi, inches = FALSE, add = TRUE, fg = points_col, bg = bubbles_col)
    points(xcoord, ycoord, pch=center_pch, cex=center_cex,  col=center_col)
  }
  if(type=="bubbles_and_approximation"){
    #profile approximation errors (column 5 of bestRepresentation)
    error<-result$bestRepresentation[,5]
    if(optimize_radius==TRUE){
      DIStz<-(outer(xcoord,xcoord, function(x,y) (x-y)^2)+outer(ycoord,ycoord, function(x,y) (x-y)^2))^0.5
      diag(DIStz)<-max(DIStz)
      rad<-min(DIStz)}
    else{
      rad<-1
    }
    raggi<-(result$bestRepresentation[,4]/pi)^0.5
    raggi<-raggi*rad*expand_radius
    raggi[which(raggi<=min_radius)]<-min_radius
    #color of each bubble, from the color scale (200 colors from 0 to error_max)
    n_col<-200
    if(rampPalette %in% c("grey", "gray")){
      scala_colori <- grey.colors(n_col, start = 0.15, end = 0.75)
    } else {
      scala_colori <- viridis(n_col, option = rampPalette)
    }
    if(is.null(error_max)) error_max <- max(error)
    colori_ordine <- scala_colori[pmin(ceiling(error / error_max * n_col) + (error == 0), n_col)]
    symbols(xcoord, ycoord,circles = raggi,inches = FALSE,add = TRUE,fg = points_col,bg = colori_ordine, lwd=lwd_circle)
    points( xcoord, ycoord,  pch = center_pch,  cex = center_cex,  col = center_col)
    #color legend on the right of the grid
    legend_image <- as.raster(matrix(scala_colori, ncol = 1))[n_col:1]
    text(x=2^k+2^(k-gV_x)*2.5, y = seq(0,2^k,l=5), labels = round(seq(0,error_max,l=5),2), cex=cex.legend)
    rasterImage(legend_image, 2^k+2^(k-gV_x)*0.75+expandX[2], 1, 2^k+2^(k-gV_x)*1.75+expandX[2],2^k)
  }

  #GRID AND PROFILE LEGEND
  #external border of the grid
  segments(c(0,0,2^k,2^k), c(0,0,0,2^k),c(0,2^k,2^k,0), c(2^k,0,2^k,2^k),lwd=Grid_lwd[1], lty=Grid_lty[1], col= Grid_col)
  #x-axis: the i-th variable of the permutation splits each of the 2^(i-1) segments of the previous level in two halves
  for(i in 1:gV_x){
    v<-2^(k-i)*seq(1,2^i,2)+0.5
    l<-rep(0,length(v))
    L<-rep(2^k,length(v))
    if(Profile==TRUE){
      L2<-rep(if(is.null(Profiles_lines)) -2^(k-gV_x)*(gV_x+1)*Profiles_expand[1] else -Profiles_lines, length(v))
      segments(v,l,v,L2,lwd=1, lty=3, col=Profiles_col)
      axis(1, pos=-2^(k-gV_x)*(i-1)*Profiles_expand[1], lwd=0 ,at=c(1-2^(k-gV_x)/k,seq(2^(k-gV_x),2^k-1,2^(k-gV_x))+0.5,2^k)-2^(k-gV_x)/2,
           labels=c(varNames[gV_x-i+1],rep(c(0,1), times=2^(gV_x-i), each=2^(i-1))), col.axis=Profiles_col, cex.axis=cex.profile, gap.axis=0.01)
    }
    if(Grid==TRUE){
      segments(v,l,v,L,lwd=Grid_lwd[i], lty=Grid_lty[i], col=Grid_col)
    }
  }
  #y-axis: the same, with the variables of the reversed permutation
  for(i in 1:gV_y){
    h<-2^(k-i)*seq(1,2^i,2)+0.5
    l<-rep(0,length(h))
    L<-rep(2^k,length(h))
    if(Profile==TRUE){
      L2<-rep(if(is.null(Profiles_lines)) -2^(k-gV_y)*(gV_y+1)*Profiles_expand[2] else -Profiles_lines, length(h))
      segments(l,h,L2,h,lwd=1, lty=3, col=Profiles_col)
      axis(2, pos=-2^(k-gV_y)*(i-1)*Profiles_expand[2], lwd=0 ,at=c(seq(2^(k-gV_y),2^k-1,2^(k-gV_y))+0.5,2^k,2^k+2^(k-gV_y)*(1+1/k))-2^(k-gV_y)/2,
           labels=c(rep(c(0,1), times=2^(gV_y-i), each=2^(i-1)),varNames[gV_x+i]),col.axis=Profiles_col, cex.axis=cex.profile, gap.axis=0.01)
    }
    if(Grid==TRUE){
      segments(l,h,L,h,lwd=Grid_lwd[i], lty=Grid_lty[i], col=Grid_col)
    }
  }
  points(xcoord, ycoord, pch=points_pch, cex=points_cex, col=center_col)
}


# Contingency map (Figure 4 of the paper): the 2^k cells of the OBE grid are colored according to the
# relative contingency of the corresponding profile, (observed - expected)/expected, the expected
# frequency being computed under independence of the variables. Negative values (fewer units than
# expected) are colored with rampPalette1, positive values with rampPalette2; a color legend is drawn
# on the right. NB: in this function the grid has gV variables on each axis.
Contingency.plot<-function(complete_matrix_profile, #output of complete_matrix_profiles_function
                           variable_priority, #variable index permutation of the embedding
                           k,#number of binary variables
                           gV=floor(k/2),#number of binary variables used to build the grid on each axis
                           Grid = TRUE, #if TRUE the grid of the tessellation is drawn
                           Grid_lty = c(1,5,2,4,3)[1:gV], #line type of the grid segments associated to each variable
                           Grid_lwd = c(3,2.5,2,1.5,1)[1:gV], #line width of the grid segments associated to each variable
                           Grid_col = "black", #color of the grid
                           Profile = TRUE, #if TRUE the variable scores are reported on both axes (profile legend)
                           Profiles_col = "black", #color of the profile legend
                           varNames = NULL, #string/expression vector with the variables' names (in the original order)
                           Profiles_expand=1, #spacing of the rows of the profile legend
                           title = "Relative Contingencies colormap", #main title of the plot
                           expandX = c(0,0), #extra space on the x-axis: c[1] = added to the left, c[2] = added to the right
                           expandY = c(0,0), #extra space on the y-axis: c[1] = added to the bottom, c[2] = added to the top
                           rampPalette1 = c("white", "darkred"), #color scale of the negative contingencies
                           rampPalette2 = c("white", "darkgreen"), #color scale of the positive contingencies
                           cex.profile=1, #size of the text of the profile legend, of the color legend and of the signs
                           signs=FALSE, #if TRUE, a + or - sign in each cell marks the sign of its relative contingency,
                           #so that positive and negative values can be told apart also in black and white
                           asp=1){ #aspect ratio of the plot: with asp = 1 the grid is drawn as a square
  if(is.null(varNames)){
    varNames<-vector()
    for(i in 1:k){
      varNames<-append(varNames, parse(text=(paste0("V[",i,"]"))))
    }
  }
  #variables' names permuted according to the variable priority
  varNames<-varNames[variable_priority]
  #DIMENSIONS OF THE PLOT AREA (as in POsetDR.plot, with room for the color legend on the right)
  ylow<-(-2^(k-gV)*gV)*(Profile == TRUE)+(Profile != TRUE)-expandY[1]
  xlow<-(-2^(k-gV)*gV)*(Profile == TRUE)+(Profile != TRUE)-expandX[1]
  yup<- (2^k+2^(k-gV)*3/2)*(Profile == TRUE)+2^k*(Profile != TRUE)+expandY[2]
  xup<- 2^k+(2^(k-gV)*3)+expandX[2]

  #EMPTY PLOT AREA
  plot(x=NULL,xlim=c(xlow,xup), ylim=c(ylow,yup),xlab=" ", ylab=" ", main=title, yaxt="n",xaxt="n",bty="n", asp=asp)
  #relative contingencies (column k+3), split into negative and positive values
  conting<-complete_matrix_profile[,k+3]
  neg_c<-which(conting<=0)
  pos_c<-which(conting>0)
  conting_neg<-conting[neg_c]
  conting_pos<-conting[pos_c]
  #colors: negative values scaled by the minimum, positive values by the maximum
  colfunc1 <- colorRampPalette(rampPalette1)
  colfunc2 <- colorRampPalette(rampPalette2)
  n_col<-200
  scala_colori1<-colfunc1(n_col)
  scala_colori2<-colfunc2(n_col)
  colori_ordine1<-scala_colori1[ceiling(conting_neg/min(conting_neg)*n_col)+(conting_neg==0)]
  colori_ordine2<-scala_colori2[ceiling(conting_pos/max(conting_pos)*n_col)+(conting_pos==0)]
  #cell of each profile: each axis is split into 2^(k-gV) cells of width 2^gV
  xcoord <- ceiling(complete_matrix_profile[,k+4]/(2^gV))
  ycoord <- ceiling(complete_matrix_profile[,k+5]/(2^gV))
  #bottom-left and top-right corners of the cells
  x_nl <- (xcoord[neg_c]-1)*(2^gV)
  y_nb <- (ycoord[neg_c]-1)*(2^gV)
  x_pl <- (xcoord[pos_c]-1)*(2^gV)
  y_pb <- (ycoord[pos_c]-1)*(2^gV)
  x_nr <- (xcoord[neg_c])*(2^gV)
  y_nt <- (ycoord[neg_c])*(2^gV)
  x_pr <- (xcoord[pos_c])*(2^gV)
  y_pt <- (ycoord[pos_c])*(2^gV)
  rect(x_nl, y_nb, x_nr, y_nt, col = colori_ordine1, border = NA)
  rect(x_pl, y_pb, x_pr, y_pt, col = colori_ordine2, border = NA)

  #color legend (minimum, 0 and maximum of the relative contingencies)
  legend_image <- as.raster(matrix(c(colfunc1(n_col/2)[(n_col/2):1],colfunc2(n_col/2)), ncol=1))[n_col:1]
  text(x=2^k+2^(k-gV)*2.7, y = seq(0,2^k,l=3), labels = c(min(conting),0,round(max(conting),2)),cex=cex.profile)
  rasterImage(legend_image, 2^k+2^(k-gV)*0.75+expandX[2], 1, 2^k+2^(k-gV)*1.75+expandX[2],2^k)

  #signs of the relative contingencies, at the center of the cells (white on the dark cells, black on the others)
  if(signs){
    colori_celle<-character(length(conting))
    colori_celle[neg_c]<-colori_ordine1
    colori_celle[pos_c]<-colori_ordine2
    luminosity<-drop(c(0.3,0.59,0.11)%*%(col2rgb(colori_celle)/255))
    text((xcoord-0.5)*2^gV, (ycoord-0.5)*2^gV, labels=ifelse(conting>0,"+","-"), col=ifelse(luminosity<0.45,"white","black"), cex=cex.profile)
  }

  #GRID AND PROFILE LEGEND
  segments(c(0,0,2^k,2^k), c(0,0,0,2^k),c(0,2^k,2^k,0), c(2^k,0,2^k,2^k),lwd=Grid_lwd[1], lty=Grid_lty[1], col= Grid_col)
  for(i in 1:gV){
    h<-v<-2^(k-i)*seq(1,2^i,2)
    l<-rep(0,length(h))
    L<-rep(2^k,length(h))
    if(Grid==TRUE){
      segments(l,h,L,h,lwd=Grid_lwd[i], lty=Grid_lty[i], col=Grid_col)
      segments(v,l,v,L,lwd=Grid_lwd[i], lty=Grid_lty[i], col=Grid_col)
    }
    if(Profile==TRUE){
      L2<-rep(-2^(k-gV)*(gV+1)*Profiles_expand,length(h))
      segments(l,h,L2,h,lwd=1, lty=3, col=Profiles_col)
      segments(v,l,v,L2,lwd=1, lty=3, col=Profiles_col)
      #x-axis: variables perm[1],...,perm[gV]
      axis(1, pos=-2^(k-gV)*(i-1)*Profiles_expand, lwd=0 ,at=c(1-2^(k-gV)/6,seq(2^(k-gV),2^k-1,2^(k-gV))+0.5,2^k)-2^(k-gV)/2,
           labels=c(varNames[gV-i+1],rep(c(0,1), times=2^(gV-i), each=2^(i-1))), col.axis=Profiles_col,cex.axis = cex.profile)
      #y-axis: variables perm[k],...,perm[gV+1] (the y coordinate is given by the reversed permutation, as in POsetDR.plot)
      axis(2, pos=-2^(k-gV)*(i-1)*Profiles_expand, lwd=0 ,at=c(seq(2^(k-gV),2^k-1,2^(k-gV))+0.5,2^k,2^k+2^(k-gV)*(1+1/6))-2^(k-gV)/2,
           labels=c(rep(c(0,1), times=2^(gV-i), each=2^(i-1)),varNames[gV+i]),col.axis=Profiles_col,cex.axis = cex.profile)
    }
  }
}
