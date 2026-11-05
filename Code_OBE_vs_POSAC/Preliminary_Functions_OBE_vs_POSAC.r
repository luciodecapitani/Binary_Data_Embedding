############################################################################
############################################################################
# Preliminary_Functions_OBE_vs_POSAC.r
#
# Packages and functions used by OBE_POSAC_comparison.r (comparison between the
# Optimal Bidimensional Embedding (OBE) and POSAC on the Financial knowledge data).
# This file is sourced at the beginning of OBE_POSAC_comparison.r.
#
# STRUCTURE OF THE FILE
#   A. Installation (only if needed) and loading of the packages
#   B. Data preparation
#   C. Comparison of the two maps (incidence matrices, indicators, classification of pairs)
#   D. Plots (OBE map, POSAC map, side-by-side maps, POSAC panel by variable)
#   E. LaTeX tables
############################################################################
############################################################################


############################################################################
# A. Installation (only if needed) and loading of the packages
############################################################################
# - poseticDataAnalysis (CRAN): Optimal Bidimensional Embedding and poset tools.
# - viridis (CRAN): color scales (used by POsetDR.plot).
# - POSAC is NOT on CRAN: it is installed from the GitHub repository of Cambridge Assessment Research
#   (https://github.com/CambridgeAssessmentResearch/POSAC) through the package remotes.
#   The results of the supplementary material were obtained with POSAC version 1.1.0, commit
#   6aaeecfc76e643bde973c613bcdef14588e0431a: this commit is installed below if POSAC is missing or if a
#   different version of POSAC is installed.
# Nothing is installed if the required packages are already available.

POSAC_commit <- "6aaeecfc76e643bde973c613bcdef14588e0431a"

for (pkg in c("poseticDataAnalysis", "viridis")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

POSAC_ok <- requireNamespace("POSAC", quietly = TRUE) &&
  identical(packageDescription("POSAC")$RemoteSha, POSAC_commit)
if (!POSAC_ok) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github(paste0("CambridgeAssessmentResearch/POSAC@", POSAC_commit))
}

library(poseticDataAnalysis)
library(POSAC)
library(viridis)


############################################################################
# B. Data preparation
############################################################################

# Frequency distribution over the observed binary profiles.
#   matrix  : data frame/matrix of the binary responses (one row per unit, one column per variable)
#   weights : vector of the sample weights of the units
# It returns a list: [[1]] matrix of the observed profiles (one row per profile, sorted by their decimal
# value), [[2]] sample counts of the profiles, [[3]] frequencies obtained by sample weights.
profiles_freq_dist_weighted <- function(matrix, weights){
  k <- ncol(matrix)
  #decimal value of the profiles
  dec_profiles <- as.matrix(matrix)%*%(2^((k:1)-1))
  observed_profiles_list <- sort(unique(dec_profiles))
  freq <- freq_W <- vector(length=length(observed_profiles_list))
  for(i in 1:length(observed_profiles_list)){
    freq_W[i] <- sum(weights*(dec_profiles==observed_profiles_list[i]))
    freq[i] <- sum(dec_profiles==observed_profiles_list[i])
  }
  obs_Profiles <- t(sapply(observed_profiles_list, function(x){ as.integer(intToBits(x))})[k:1,])
  return(list(obs_Profiles, freq, freq_W))
}

# Converts an integer into its binary string (e.g. 5 -> "000101" with n_bits = 6). It is used to turn the
# profile identifiers returned by OptimalBidimensionalEmbedding into profile names, consistent with the
# profile names returned by POSAC.
int2bin <- function(x, n_bits = 8) {
  bits <- as.integer(intToBits(x))[1:n_bits]
  paste(rev(bits), collapse = "")
}


############################################################################
# C. Comparison of the two maps
############################################################################

# Incidence (characteristic) matrix of the bidimensional poset induced by a map, given the profile
# names and their X, Y coordinates (used both for the POSAC and for the OBE map):
# M[i, j] = 1 if and only if X[i] <= X[j] and Y[i] <= Y[j].
intersection_incidence <- function(Structuple, X, Y, reflexive = TRUE) {
  n <- length(Structuple)
  if (length(X) != n || length(Y) != n)
    stop("Structuple, X and Y must have the same length")
  M <- (outer(X, X, "<=") & outer(Y, Y, "<=")) * 1L
  if (!reflexive) diag(M) <- 0L
  dimnames(M) <- list(as.character(Structuple), as.character(Structuple))
  M
}

# Indicators of the quality of a map (Section 3 of the supplementary material): numbers and percentages
# of the pairs of distinct profiles whose relation is correctly reproduced by the map, overall and
# separately for the comparable and the incomparable pairs of the input poset.
#   M_t   : incidence matrix of the input (observed) poset
#   M_a   : incidence matrix of the poset induced by the map
#   freqs : profile frequencies; if NULL each pair has weight 1 ("unweighted"), otherwise the pair {i, j}
#           has weight freqs[i]*freqs[j] ("weighted"). The percentages do not depend on the scale of freqs.
# Convention: M[i, j] = 1 if and only if element j dominates (or coincides with) element i.
compare_incidence <- function(M_t, M_a, freqs = NULL) {
  M_t <- as.matrix(M_t); M_a <- as.matrix(M_a)
  n <- nrow(M_t)
  if (any(dim(M_t) != n) || any(dim(M_a) != n))
    stop("M_t and M_a must be square matrices of the same dimension")
  if (!is.null(rownames(M_t)) && !is.null(rownames(M_a))) {
    if (!setequal(rownames(M_t), rownames(M_a)))
      stop("M_t and M_a do not have the same row/column names")
    M_a <- M_a[rownames(M_t), rownames(M_t)]
  }
  if (is.null(freqs)) freqs <- rep(1, n)
  if (length(freqs) != n) stop("the length of freqs must be equal to nrow(M_t)")

  # unordered pairs of distinct elements (i < j)
  idx <- which(upper.tri(M_t), arr.ind = TRUE)
  i <- idx[, 1]; j <- idx[, 2]
  w <- freqs[i] * freqs[j]

  comp_t   <- (M_t[cbind(i, j)] + M_t[cbind(j, i)]) > 0   # comparable in the input poset
  incomp_t <- !comp_t                                     # incomparable in the input poset
  correct  <- (M_t[cbind(i, j)] == M_a[cbind(i, j)]) &
    (M_t[cbind(j, i)] == M_a[cbind(j, i)])

  data.frame(
    npairs              = sum(w),
    npairscorrect       = sum(w * correct),
    pccorrect           = 100 * sum(w * correct) / sum(w),
    npairscomp          = sum(w * comp_t),
    npairscompcorrect   = sum(w * (comp_t & correct)),
    pccompcorrect       = 100 * sum(w * (comp_t & correct)) / sum(w * comp_t),
    npairsincomp        = sum(w * incomp_t),
    npairsincompcorrect = sum(w * (incomp_t & correct)),
    pcincompcorrect     = 100 * sum(w * (incomp_t & correct)) / sum(w * incomp_t)
  )
}

# Classification of the pairs of distinct profiles according to how a map represents their relation
# in the input (observed) poset (Table S2 of the supplementary material):
#   comparable pairs   -> preserved (same dominance), turned into incomparable, reversed;
#   incomparable pairs -> preserved, turned into comparable.
# The numbers of pairs are unweighted; the last two entries are the percentages of comparable pairs
# turned into incomparable ones and of incomparable pairs turned into comparable ones, each pair being
# weighted by the product of the frequencies of its two profiles.
# M_t, M_a and freqs as in compare_incidence (freqs is required here).
pair_classification <- function(M_t, M_a, freqs) {
  M_t <- as.matrix(M_t); M_a <- as.matrix(M_a)[rownames(M_t), colnames(M_t)]
  idx <- which(upper.tri(M_t), arr.ind = TRUE)
  i <- idx[, 1]; j <- idx[, 2]
  w <- freqs[i] * freqs[j]
  # relation of i and j: 1 if i is below j, -1 if j is below i, 0 if incomparable (distinct profiles)
  rel <- function(M) ifelse(M[cbind(i, j)] == 1, 1, ifelse(M[cbind(j, i)] == 1, -1, 0))
  r_t <- rel(M_t); r_a <- rel(M_a)
  comp <- r_t != 0; inc <- !comp
  c(comp_preserved      = sum(comp & r_a == r_t),
    comp_to_incomp      = sum(comp & r_a == 0),
    comp_reversed       = sum(comp & r_a == -r_t),
    incomp_preserved    = sum(inc & r_a == 0),
    incomp_to_comp      = sum(inc & r_a != 0),
    w_pc_comp_to_incomp = 100 * sum(w[comp & r_a == 0]) / sum(w[comp]),
    w_pc_incomp_to_comp = 100 * sum(w[inc & r_a != 0]) / sum(w[inc]))
}

# Incomparable pairs of the input (observed) poset, grouped by the difference in the number of 1s of their
# two profiles (numbers quoted in the paragraph "No context" of the supplementary material).
#   M_t, freqs : as in pair_classification (the row names of M_t must be the profile names, e.g. "010011")
#   maps       : named list of the incidence matrices of the posets induced by the maps
#   D          : matrix of the degrees of dominance of the profiles, with the same names as M_t
# It returns a data frame with one row per difference: the degree of dominance of the pairs (the function
# stops if it is not the same for all of the pairs with that difference), their number and, for each map,
# the number of pairs made comparable (<map>_to_comp), their percentage counting each pair once (<map>_pc)
# and their percentage weighting each pair by the product of the frequencies of its profiles (<map>_w_pc).
incomparable_by_difference <- function(M_t, maps, freqs, D) {
  M_t <- as.matrix(M_t)
  idx <- which(upper.tri(M_t), arr.ind = TRUE)
  i <- idx[, 1]; j <- idx[, 2]
  inc <- (M_t[cbind(i, j)] + M_t[cbind(j, i)]) == 0   # incomparable in the input poset
  i <- i[inc]; j <- j[inc]
  ones <- sapply(strsplit(rownames(M_t), ""), function(x) sum(x == "1"))
  d <- abs(ones[i] - ones[j])
  w <- freqs[i] * freqs[j]
  degree <- pmax(D[cbind(i, j)], D[cbind(j, i)])
  if (any(tapply(degree, d, function(x) diff(range(x))) > 1e-8))
    stop("the degree of dominance is not determined by the difference in the number of 1s")
  out <- data.frame(difference = sort(unique(d)),
                    degree     = as.vector(tapply(degree, d, mean)),
                    pairs      = as.vector(table(d)))
  for (m in names(maps)) {
    M_a <- as.matrix(maps[[m]])[rownames(M_t), colnames(M_t)]
    to_comp <- (M_a[cbind(i, j)] + M_a[cbind(j, i)]) > 0   # comparable in the map
    out[[paste0(m, "_to_comp")]] <- as.vector(tapply(to_comp, d, sum))
    out[[paste0(m, "_pc")]]      <- 100 * as.vector(tapply(to_comp, d, mean))
    out[[paste0(m, "_w_pc")]]    <- 100 * as.vector(tapply(w * to_comp, d, sum) / tapply(w, d, sum))
  }
  out
}


############################################################################
# D. Plots
############################################################################

# OBE map (function of the paper). In this comparison it is used with type = "bubbles".
POsetDR.plot<-function(result, #output of OptimalBidimensionalEmbedding
                       N, #number of statistical units
                       k,#number of binary variables
                       gV_x=ceiling(k/2),#number of binary variables considered in building the grid on the x-axis
                       gV_y=k-ceiling(k/2),#number of binary variables considered in building the grid on the y-axis
                       Grid = TRUE, #if TRUE a grid based on gV variables is added on both axes
                       Grid_lty = c(1,5,2,4,3)[1:gV_x], #Line type of the grid segments associated to each variable in the grid
                       Grid_lwd = c(3,2.5,2,1.5,1)[1:gV_x], #Line width of the grid segments associated to each variable in the grid
                       Grid_col = "blue", #color of the grid
                       Profile = TRUE, #if TRUE a legend is added on both axes in order to recognize the profiles associated to each square in the grid
                       Profiles_col = "black", #color of the profile legend
                       varNames = NULL, #string /expression vector with the variables' names
                       Profiles_expand=c(1,1), # two dimension vector. The first element regards x axis, the second y axis.
                       #if greater (lower) than one, the space between the lines in the profile legend is expanded (reduced).
                       title = "Optimal profile map", #main title of the plot
                       type=c("point", "jitter_square", "jitter_circle", "bubbles", "bubbles_and_approximation"), #type of plot
                       expandX = c(0,0), #add extra space to the default x-axis: c[1]=space added to the left, c[2]=space added to the right
                       expandY = c(0,0), #add extra space to the default y-axis: c[1]=space added to the bottom, c[2]=space added to the top
                       expand_radius=1, #enlarge the radius of bubbles
                       points_col="black", #when type = "point", "jitter_square", "jitter_circle", it indicates the color of the points drawn, when
                       # type= "bubbles", "bubbles_and_approximation" it indicates the color of the border of the bubbles
                       center_col= "red", #when type = "jitter_square", "jitter_circle","bubbles", "bubbles_and_approximation"
                       #it indicates the color of the points drawn in the center of the clouds of points/bubbles
                       optimize_radius=FALSE,
                       bubbles_col = "gray", #when type ="bubbles" it indicates the bubbles' color
                       points_pch = 16, #when type = "point", "jitter_square", "jitter_circle", it indicates the type of the points drawn
                       center_pch = 16, #when type = "jitter_square", "jitter_circle", "bubbles", "bubbles_and_approximation",
                       #it indicates the type of the points drawn at the center of the clouds of points/bubbles
                       points_cex = 0.5, #when type = "point", "jitter_square", "jitter_circle", it indicates the dimension of the points drawn
                       center_cex = 0.5, #when type = "jitter_square", "jitter_circle", "bubbles", "bubbles_and_approximation",
                       #it indicates the dimension of the points drawn at the center of the clouds of points/bubbles
                       jitter_amount = 1, #when type = "jitter_square", "jitter_circle" it is proportional to the dispersion of the clouds of points
                       #around the centers
                       rampPalette = "viridis",
                       bubble_sizes = NULL,
                       cex.profile=1,
                       min_radius=0,
                       lwd_circle=1,
                       error_max=NULL, #when type = "bubbles_and_approximation", upper bound of the color scale (if NULL, max of the errors)
                       xlim=NULL, #optional limits of the plot area (they override the default ones, computed from expandX and
                       ylim=NULL, #expandY); e.g. xlim = ylim = c(0, 2^k) makes the plot area coincide with the grid
                       Profiles_lines=NULL){ #length (in grid units) of the dotted lines of the profile legend;
                       #if NULL, 2^(k-gV)*(gV+1)*Profiles_expand

  if(length(result)==3){
    names(result)[c(2,3)]<-c("bestVariablesPriority","bestRepresentation")}
  #default variables' names: V[i] is the variable on the i-th column on the original dataset
  if(is.null(varNames)){
    varNames<-vector()
    for(i in 1:k){
      varNames<-append(varNames, parse(text=(paste0("V[",i,"]"))))
    }
  }
  #permuting variables names according to variables priority
  varNames<-varNames[result$bestVariablesPriority]

  #DIMENSIONS OF PLOT AREA
  ylow<-(-2^(k-gV_y)*gV_y)*(Profile == TRUE)+(Profile != TRUE)-expandY[1]
  xlow<-(-2^(k-gV_x)*gV_x)*(Profile == TRUE)+(Profile != TRUE)-expandX[1]
  #NB: type "bubbles and approximation" is a special case for xup
  yup<- (2^k+2^(k-gV_y)*3/2)*(Profile == TRUE)+2^k*(Profile != TRUE)+expandY[2]
  xup<- 2^k*(type != "bubbles_and_approximation")+(2^k+2^(k-gV_x)*3)*(type == "bubbles_and_approximation")+expandX[2]
  if(!is.null(xlim)){ xlow <- xlim[1]; xup <- xlim[2] }
  if(!is.null(ylim)){ ylow <- ylim[1]; yup <- ylim[2] }

  #DRAW THE PLOT AREA
  plot(x=NULL,xlim=c(xlow,xup), ylim=c(ylow,yup),xlab=" ", ylab=" ", main=title, yaxt="n",xaxt="n",bty="n")

  #ADDING THE POINTS (with optional jittering or bubbles) TO THE PLOT AREA
  xcoord<-result$bestRepresentation[,2]+0.5*(result$bestRepresentation[,2]==1)-0.5*(result$bestRepresentation[,2]==2^k)
  ycoord<-result$bestRepresentation[,3]+0.5*(result$bestRepresentation[,3]==1)-0.5*(result$bestRepresentation[,3]==2^k)
  if(type=="point"){
    points(xcoord, ycoord, pch=points_pch, cex=points_cex, col=points_col)
  }
  if(type=="jitter_square"){
    npoints<-round(N*result$bestRepresentation[,4]/sum(result$bestRepresentation[,4]))
    x_unit<-rep(xcoord,npoints)
    y_unit<-rep(ycoord,npoints)
    Coordin1<-jitter(x_unit,amount=jitter_amount)*(rep(npoints,npoints)>1)+x_unit*(rep(npoints,npoints)==1)
    Coordin2<-jitter(y_unit,amount=jitter_amount)*(rep(npoints,npoints)>1)+y_unit*(rep(npoints,npoints)==1)
    points(Coordin1, Coordin2, pch=points_pch,  cex=points_cex, col=points_col)
    points(xcoord, ycoord, pch=center_pch,  cex=center_cex, col=center_col)
  }
  if(type=="jitter_circle"){
    npoints<-round(N*result$bestRepresentation[,4]/sum(result$bestRepresentation[,4]))
    x_unit<-rep(xcoord,npoints)
    y_unit<-rep(ycoord,npoints)
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
    raggi<-(result$bestRepresentation[,4]/pi)^0.5
    raggi<-raggi*rad*expand_radius
    raggi[which(raggi<=min_radius)]<-min_radius
    symbols(xcoord, ycoord, circles = raggi, inches = FALSE, add = TRUE, fg = points_col, bg = bubbles_col, lwd = lwd_circle)
    points(xcoord, ycoord, pch=center_pch, cex=center_cex,  col=center_col)
  }
  if(type=="bubbles_and_approximation"){
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
    legend_image <- as.raster(matrix(scala_colori, ncol = 1))[n_col:1]
    text(x=2^k+2^(k-gV_x)*2.5, y = seq(0,2^k,l=5), labels = round(seq(0,error_max,l=5),2), cex=cex.profile)
    rasterImage(legend_image, 2^k+2^(k-gV_x)*0.75+expandX[2], 1, 2^k+2^(k-gV_x)*1.75+expandX[2],2^k)
  }
  #ADDING THE (optional) GRID AND (optional) PROFILES' LEGEND
  #external borders of the grid
  segments(c(0,0,2^k,2^k), c(0,0,0,2^k),c(0,2^k,2^k,0), c(2^k,0,2^k,2^k),lwd=Grid_lwd[1], lty=Grid_lty[1], col= Grid_col)
  #grid on the x-axis
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
  #grid on the y-axis
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


# POSAC map, with the same graphical logic as POsetDR.plot: each observed profile is a bubble centred at
# its POSAC coordinates (X, Y), with area proportional to the profile frequency. The POSAC axes have no
# variable-based semantics, so there is no tessellation: the map is the square [0, n+1] x [0, n+1],
# n = number of observed profiles.
POSAC.plot<-function(posac_out, #output of POSAC()
                     freqs = posac_out$freqs, #profile frequencies (same order as posac_out$Patterns); internally rescaled to sum 1
                     type = "bubbles", #type of plot (only "bubbles" in this comparison)
                     title = "POSAC map", #main title of the plot
                     expandX = c(0,0), #extra space on the x-axis: c[1]=left, c[2]=right
                     expandY = c(0,0), #extra space on the y-axis: c[1]=bottom, c[2]=top
                     expand_radius = 1, #enlarge the radius of bubbles. To get bubbles comparable with a POsetDR.plot map
                     #drawn in a panel of the same size, use expand_radius_OBE*(n+1)/2^k
                     min_radius = 0, #minimum radius of the bubbles
                     points_col = "black", #color of the border of the bubbles
                     bubbles_col = "gray", #bubbles' color (a single color or one color per profile)
                     center_col = "red", #color of the points drawn at the center of the bubbles
                     center_pch = 16, #type of the points drawn at the center of the bubbles
                     center_cex = 0.5, #dimension of the points drawn at the center of the bubbles
                     Frame = TRUE, #if TRUE the border of the map is drawn
                     Frame_col = "darkgray", #color of the border
                     Frame_lwd = 2, #line width of the border
                     axes = TRUE, #if TRUE the ranks X and Y of POSAC are reported on the axes
                     axes_by = 5, #step of the ticks on the axes
                     labels = FALSE, #if TRUE the profiles' binary strings are written next to the bubbles
                     labels_cex = 0.6, #dimension of the profiles' labels
                     labels_col = "black", #color of the profiles' labels
                     cex.profile = 1, #dimension of the numbers on the axes
                     lwd_circle = 1, #line width of the border of the bubbles
                     xlim = NULL, #optional limits of the plot area (they override the default ones);
                     ylim = NULL, #e.g. to give the map the same size as an OBE map drawn in a panel of the same size
                     asp = 1){ #aspect ratio of the plot (as in plot.default); use NA together with xlim and ylim
  xcoord <- posac_out$X
  ycoord <- posac_out$Y
  n <- length(xcoord)
  if(length(freqs) != n) stop("freqs must have one value for each profile of posac_out")
  freqs <- as.numeric(freqs)/sum(freqs)

  #RADIUS OF THE BUBBLES (as in POsetDR.plot: area proportional to the frequency)
  raggi <- (freqs/pi)^0.5
  raggi <- raggi*expand_radius
  raggi[which(raggi <= min_radius)] <- min_radius

  #DIMENSIONS OF PLOT AREA
  s <- (n+1)/8
  #extra space around the frame, so that bubbles on the border are not cut and do not cover the axes
  pad <- max(0, max(c(raggi-xcoord, raggi-ycoord, xcoord+raggi-(n+1), ycoord+raggi-(n+1))))
  axes_pos <- -pad-s*0.1
  xlow <- -pad-expandX[1]-(axes == TRUE)*s*0.6
  ylow <- -pad-expandY[1]-(axes == TRUE)*s*0.6
  xup <- n+1+pad+expandX[2]
  yup <- n+1+pad+expandY[2]
  if(!is.null(xlim)){ xlow <- xlim[1]; xup <- xlim[2] }
  if(!is.null(ylim)){ ylow <- ylim[1]; yup <- ylim[2] }

  #DRAW THE PLOT AREA AND THE BUBBLES (from the largest to the smallest, so that small bubbles are not hidden)
  plot(x=NULL, xlim=c(xlow,xup), ylim=c(ylow,yup), xlab=" ", ylab=" ", main=title, yaxt="n", xaxt="n", bty="n", asp=asp)
  bubbles_col <- rep(bubbles_col, length.out = n)
  ord <- order(raggi, decreasing = TRUE)
  symbols(xcoord[ord], ycoord[ord], circles = raggi[ord], inches = FALSE, add = TRUE, fg = points_col, bg = bubbles_col[ord], lwd = lwd_circle)
  points(xcoord, ycoord, pch=center_pch, cex=center_cex, col=center_col)

  #FRAME, AXES AND (optional) PROFILES' LABELS
  if(Frame == TRUE){
    segments(c(0,0,n+1,n+1), c(0,0,0,n+1), c(0,n+1,n+1,0), c(n+1,0,n+1,n+1), lwd=Frame_lwd, col=Frame_col)
  }
  if(axes == TRUE){
    tk <- unique(c(1, seq(axes_by, n, by=axes_by), n))
    axis(1, pos=axes_pos, at=tk, labels=tk, cex.axis=cex.profile, lwd=0, lwd.ticks=1, tcl=-0.3, mgp=c(3,0.4,0))
    axis(2, pos=axes_pos, at=tk, labels=tk, cex.axis=cex.profile, lwd=0, lwd.ticks=1, tcl=-0.3, mgp=c(3,0.5,0), las=1)
  }
  if(labels == TRUE){
    text(xcoord, ycoord, labels=posac_out$Patterns, pos=4, cex=labels_cex, col=labels_col, offset=0.3)
  }
}


# Panel of POSAC "bubbles" maps, one for each selected variable, in the spirit of itemdiagram() of the
# POSAC package: in the panel of variable V, each bubble is colored according to the score of V in the
# corresponding profile (col0 if V = 0, col1 if V = 1). Figure S2 of the supplementary material.
POSAC.itemdiagram.bubbles<-function(posac_out, #output of POSAC()
                                    freqs = posac_out$freqs, #profile frequencies (bubble areas), same order as posac_out$Patterns
                                    vars = 1:ncol(posac_out$patmat), #indexes of the variables (columns of posac_out$patmat) to be plotted
                                    var_names = colnames(posac_out$patmat), #names of all the variables (character or expression vector)
                                    col0 = "blue", #bubble color when the variable scores 0
                                    col1 = "red", #bubble color when the variable scores 1
                                    legend01 = TRUE, #if TRUE a common legend (0/1 colors) is added at the bottom of the panel
                                    panel_mai = c(0.38, 0.42, 0.30, 0.16), #margins of each panel, in inches (bottom, left, top, right):
                                    #with xlim = ylim = c(0, n+1), passed through ..., the plot area of each panel coincides with its square
                                    legend_height = 0.28, #space for the legend at the bottom of the figure (inches)
                                    legend_cex = 10/12, #size of the text of the legend
                                    ...){ #further arguments passed to POSAC.plot (e.g. expand_radius, min_radius, center_cex)
  patmat <- as.matrix(posac_out$patmat)
  if(is.null(var_names)) var_names <- paste0("V", 1:ncol(patmat))
  dim1 <- ceiling(sqrt(length(vars)))
  dim2 <- ceiling(length(vars)/dim1)
  op <- par(mfrow = c(dim2, dim1), omi = c(legend_height*legend01, 0, 0, 0), mai = panel_mai, xaxs = "i", yaxs = "i", xpd = NA)
  par(cex = 1, cex.main = 1)   #after mfrow, which would reduce the font size
  for(v in vars){
    POSAC.plot(posac_out, freqs = freqs, title = var_names[v],
               bubbles_col = ifelse(patmat[, v] == 1, col1, col0), ...)
  }
  if(legend01 == TRUE){
    par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    legend("bottom", legend = c("score 0", "score 1"), pt.bg = c(col0, col1), pch = 21, pt.cex = 1.5,
           horiz = TRUE, bty = "n", xpd = NA, cex = legend_cex)
  }
  par(op)
}


# POSAC map (left) and OBE map (right) side by side, with bubbles of area proportional to the profile
# frequencies (Figure S1 of the supplementary material). In each panel the plot area coincides with the map
# (the (n+1)x(n+1) square of POSAC and the 2^k x 2^k grid of OBE), whose side is given in inches, so that the
# two maps have the same physical size; axes and legends are drawn in the margins (xpd = NA), and the bubble
# radii are rescaled so that the two maps can be compared directly.
#   obe, posac : outputs of OptimalBidimensionalEmbedding and POSAC
#   freqs_posac: profile frequencies for the POSAC bubbles (same order as posac$Patterns)
#   N, k, varNames: number of units, number of variables, variable names (for the OBE map)
#   expand_radius: radius scale of the OBE bubbles
#   mai, side  : margins of each panel (bottom, left, top, right) and side of the maps (inches)
#   cex.text   : size of the text on the axes (10/12 = 10 points, with a device font of 12 points)
#   row_spacing: distance between the rows/columns of the profile legend of the OBE map (inches)
OBE_POSAC_maps <- function(obe, posac, freqs_posac, N, k, varNames, expand_radius = 30,
                           mai = c(0.53, 0.55, 0.38, 0.28), side = 2.42, cex.text = 10/12, row_spacing = 0.14) {
  gV <- ceiling(k/2)
  n_obs <- length(posac$X)
  unit <- side/2^k   #inches per grid unit of the OBE map
  op <- par(mfrow = c(1, 2), mai = mai, xaxs = "i", yaxs = "i", xpd = NA, cex = 1, cex.main = 1, mgp = c(3, 0.2, 0))
  POSAC.plot(posac, freqs = freqs_posac, bubbles_col = "gray", title = "POSAC", center_col = "white", center_cex = 0.05,
             expand_radius = expand_radius*(n_obs+1)/2^k, min_radius = 0.4, lwd_circle = 0.5, labels = FALSE,
             cex.profile = cex.text, axes_by = 10, Frame_lwd = 1.1, xlim = c(0, n_obs+1), ylim = c(0, n_obs+1), asp = NA)
  POsetDR.plot(obe, N = N, k = k, gV_x = gV, gV_y = k-gV, Grid = TRUE, Grid_lty = c(1,5,2,4,3)[1:gV], Grid_lwd = rep(1.1, gV),
               Grid_col = "darkgray", Profile = TRUE, Profiles_col = "black", Profiles_expand = rep(row_spacing/unit/2^(k-gV), 2),
               Profiles_lines = (mai[1] - 0.02)/unit, cex.profile = cex.text, varNames = varNames, title = "OBE", type = "bubbles",
               xlim = c(0, 2^k), ylim = c(0, 2^k), points_col = "black", bubbles_col = "gray", lwd_circle = 0.5,
               center_col = "white", center_pch = 16, center_cex = 0.05, points_cex = 0.05, expand_radius = expand_radius, min_radius = 0.4)
  par(op)
}


############################################################################
# E. LaTeX tables
############################################################################

fmt <- function(x, d) formatC(x, format = "f", digits = d, big.mark = ",")

# Table S1: indicators based on characteristic matrices (output: lines of a LaTeX table).
#   tab: data frame with columns method, pairs and the columns returned by compare_incidence
latex_table_S1 <- function(tab, weighted_digits = 4) {
  lines <- c("\\begin{table}[ht]", "\\centering",
             paste0("\\caption{Indicators based on characteristic matrices, application with frequencies obtained by sample weights. ",
                    "$n$: (weighted) number of pairs of distinct observed profiles; $n_c$: (weighted) number of pairs whose relation ",
                    "is correctly reproduced by the map; \\%: $100\\,n_c/n$.}"),
             "\\label{tab:incW}", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{llrrrrrrrrr}", "\\hline",
             " & & \\multicolumn{3}{c}{All pairs} & \\multicolumn{3}{c}{Comparable pairs} & \\multicolumn{3}{c}{Incomparable pairs}\\\\",
             "Pairs & Method & $n$ & $n_{c}$ & \\% & $n$ & $n_{c}$ & \\% & $n$ & $n_{c}$ & \\% \\\\", "\\hline")
  for (r in 1:nrow(tab)) {
    d <- if (tab$pairs[r] == "unweighted") 0 else weighted_digits
    lines <- c(lines, paste0(if (r %in% c(1,3)) tab$pairs[r] else "", " & ", tab$method[r], " & ",
      fmt(tab$npairs[r], d), " & ", fmt(tab$npairscorrect[r], d), " & ", fmt(tab$pccorrect[r], 2), " & ",
      fmt(tab$npairscomp[r], d), " & ", fmt(tab$npairscompcorrect[r], d), " & ", fmt(tab$pccompcorrect[r], 2), " & ",
      fmt(tab$npairsincomp[r], d), " & ", fmt(tab$npairsincompcorrect[r], d), " & ", fmt(tab$pcincompcorrect[r], 2), "\\\\"))
    if (r == 2) lines <- c(lines, "\\hline")
  }
  c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
}

# Table S2: classification of the pairs (output: lines of a LaTeX table).
#   tab: matrix with rows POSAC and OBE and the columns returned by pair_classification
#   n_comp, n_incomp: numbers of comparable and incomparable pairs of the observed poset
latex_table_S2 <- function(tab, n_comp, n_incomp) {
  lines <- c("\\begin{table}[ht]", "\\centering",
             paste0("\\caption{Classification of the pairs of distinct observed profiles according to how the two maps represent ",
                    "their relation in the observed poset (frequencies obtained by sample weights). \\textit{Comparable pairs} of the ",
                    "observed poset: Pres.\\ (preserved), the map shows the same dominance; T.I.\\ (turned into incomparable), the two ",
                    "profiles are incomparable in the map; Rev.\\ (reversed), the map shows the opposite dominance (the dominated profile ",
                    "becomes the dominating one). \\textit{Incomparable pairs} of the observed poset: Pres.\\ (preserved), the two profiles ",
                    "are incomparable also in the map; T.C.\\ (turned into comparable), one of the two profiles dominates the other in the ",
                    "map. The numbers of pairs are unweighted; the last two columns report the percentage of the comparable pairs turned ",
                    "into incomparable ones (C$\\to$I) and of the incomparable pairs turned into comparable ones (I$\\to$C), each pair ",
                    "being weighted by the product of the frequencies of its two profiles.}"),
             "\\label{tab:fate}", "\\small", "\\begin{tabular}{lrrrrrrr}", "\\hline",
             paste0(" & \\multicolumn{3}{c}{Comparable pairs (", fmt(n_comp, 0), ")} & \\multicolumn{2}{c}{Incomparable pairs (",
                    fmt(n_incomp, 0), ")} & \\multicolumn{2}{c}{Weighted \\%}\\\\"),
             "Method & Pres. & T.I. & Rev. & Pres. & T.C. & C$\\to$I & I$\\to$C\\\\", "\\hline")
  for (m in rownames(tab)) {
    lines <- c(lines, paste0(m, " & ", fmt(tab[m, "comp_preserved"], 0), " & ", fmt(tab[m, "comp_to_incomp"], 0), " & ", fmt(tab[m, "comp_reversed"], 0), " & ",
      fmt(tab[m, "incomp_preserved"], 0), " & ", fmt(tab[m, "incomp_to_comp"], 0), " & ", fmt(tab[m, "w_pc_comp_to_incomp"], 2), " & ",
      fmt(tab[m, "w_pc_incomp_to_comp"], 2), "\\\\"))
  }
  c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
}
