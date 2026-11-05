############################################################################
############################################################################
# Financial_literacy_JCGS.r
#
# Reproduces the results of the Financial knowledge example (Section "The Financial knowledge example")
# of the paper "Bidimensional Embedding and Visualization of Multidimensional Ordinal Binary Data",
# on the data of the Bank of Italy 2020 survey on financial literacy (file iacofi.csv):
#   - Table 1  : distribution of the number of entries equal to 1 in the knowledge profiles
#   - Figure 2 : synthesis of the computations (four panels)
#   - Figure 3 : the optimal embedding (jittered clouds of units; bubbles and approximation errors)
#   - Figure 4 : relative contingencies
#   - Figure 5 : maps of the sub-samples with different education degrees
#   - the numbers quoted in the text (printed at the end and saved in numbers_quoted_in_the_text.txt)
# The whole workspace is saved in Financial_literacy_JCGS.RData.
#
# STRUCTURE OF THE PROGRAM
#   0. Packages and functions
#   1. Loading and preparing the data
#   2. Frequency distributions over the profiles (whole sample and by education degree)
#   3. Optimal Bidimensional Embedding
#   4. Table 1 and numbers quoted in the text
#   5. Figure 2: synthesis of the computations
#   6. Figure 3: the optimal embedding
#   7. Figure 4: relative contingencies
#   8. Figure 5: sub-samples by education degree
#   9. Saving of tables, numbers and workspace
#
# Run the script with the working directory set to the folder containing this file, the data file
# iacofi.csv and Preliminary_Functions_Financial_literacy.r, e.g. setwd("path/to/Code_Financial_literacy").
############################################################################
############################################################################


############################################################################
# 0. Packages and functions
############################################################################
# installs (only if needed) and loads poseticDataAnalysis and viridis, and defines the functions used below
source("Preliminary_Functions_Financial_literacy.r")

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

#names of the six knowledge variables in the paper (K1,...,K6 correspond to qk4, qk5, qk6, qk7_1, qk7_2, qk7_3)
varNames <- c("K1","K2","K3","K4","K5","K6")   #plain text, as in the paper: subscripts would be printed below 10 points


############################################################################
# 1. Loading and preparing the data
############################################################################

DT <- read.csv("iacofi.csv")

#selecting the variables of interest
SUBDT <- DT[,c("qk4", "qk5", "qk6", "qk7_1","qk7_1alt", "qk7_2", "qk7_3", "qk7_3alt","qd1","qd7","area5","reddito","edu","capofam")]

#contingency table of the answers to qk5 (K2, simple interest) and qk6 (K3, compound interest), keeping the
#correct answers (102 and 1) and the codes of "don't know" (-97) and "refused/missing" (-99), all the other
#(wrong) answers being coded 0. It is used for the comment on the units who can compute compound, but not
#simple, interest (Section "Interpretation and information content of the visualization").
DDD <- SUBDT[,c(2,3)]
DDD[,1] <- DDD[,1]%in%c(-99,-97,102)*DDD[,1]
DDD[,2] <- DDD[,2]%in%c(-99,-97,1)*DDD[,2]
TT <- table(DDD)
TT

#union of the answers to questions qk7_1 and qk7_1alt, and to qk7_3 and qk7_3alt (the same questions,
#asked in two versions to different units)
SUBDT[which(is.na(SUBDT[,7])),7] <- SUBDT[which(is.na(SUBDT[,7])),8]
SUBDT[which(is.na(SUBDT[,4])),4] <- SUBDT[which(is.na(SUBDT[,4])),5]
SUBDT$qk7_1alt <- SUBDT$qk7_3alt <- NULL
#coding of the correct answers: qk4: 0; qk5: 102; qk6: 1; qk7_1 (merged): 1; qk7_2: 1; qk7_3 (merged): 1
corrette <- c(0,102,1,1,1,1)

#sample weights (normalized to sum 1) and sample size
weights <- DT$wght
weights <- weights/sum(weights)
N <- nrow(SUBDT)

#binary data matrix: 1 = correct answer, 0 = wrong answer, "don't know", refused or missing answer
MIS <- as.data.frame(matrix(-2, nrow=N, ncol=6))   #-2 is a fictitious value to initialize the matrix
colnames(MIS) <- c("qk4", "qk5", "qk6", "qk7_1", "qk7_2", "qk7_3")
for(h in 1:nrow(SUBDT)){
  MIS[h,] <- (SUBDT[h,1:6]==corrette)
}
k <- ncol(MIS)


############################################################################
# 2. Frequency distributions over the profiles (whole sample and by education degree)
############################################################################

#whole sample: B[[1]] observed profiles, B[[2]] sample counts, B[[3]] frequencies obtained by sample weights
B <- profiles_freq_dist_weighted(MIS, weights)

#by education degree (variable edu of the survey: 1 = bachelor/master or higher; 2-3 = high school;
#4-5 = middle school; 6-7 = primary school; 8 = no degree). Primary school and no degree are merged
#into a single class ("Elementary or no degree"), since the class with no degree is very small.
edu_labels <- c("Bachelor/master or higher", "High school", "Middle school", "Elementary or no degree")
edu <- rep(NA_character_, N)
edu[which(SUBDT$edu==1)] <- edu_labels[1]
edu[which((SUBDT$edu>=2) & (SUBDT$edu<=3))] <- edu_labels[2]
edu[which((SUBDT$edu>=4) & (SUBDT$edu<=5))] <- edu_labels[3]
edu[which((SUBDT$edu>=6) & (SUBDT$edu<=8))] <- edu_labels[4]
table(edu, useNA = "ifany")
#B_edu[[5]]: frequencies obtained by sample weights, one column per education degree
B_edu <- profiles_freq_dist_weighted_classified(MIS, weights, edu, edu_labels)


############################################################################
# 3. Optimal Bidimensional Embedding
############################################################################
# exhaustive search over the 6!/2 = 360 admissible embeddings (reversed pairs of lexicographic linear
# extensions), with the frequencies obtained by sample weights

obs_Profiles <- B[[1]]
freq_W <- B[[3]]
result <- OptimalBidimensionalEmbedding(obs_Profiles, freq_W, output_every_sec = 1)

result$bestVariablesPriority   #optimal variable index permutation
result$bestLossValue           #global loss of the optimal embedding
#result$bestRepresentation: one row per observed profile, with columns
#  profiles (decimal identifier), x and y (coordinates in the grid), weights (frequency), error (profile loss)
head(result$bestRepresentation)


############################################################################
# 4. Table 1 and numbers quoted in the text
############################################################################

#Table 1: distribution of the number of entries equal to 1 in the units' profiles (sample weights)
n_ones <- rowSums(MIS)
Table_1 <- tapply(weights, factor(n_ones, levels = 0:k), sum)
round(Table_1, 2)

#comparabilities and incomparabilities among the observed profiles, in the input lattice and in the
#optimal embedding
comp_input_obs  <- count_comparabilities(obs_Profiles)
comp_output_obs <- count_comparabilities(result$bestRepresentation[, c("x", "y")])
#the same, considering all of the 2^k possible profiles (coordinates of all of the profiles in the
#embedding generated by the optimal variable index permutation)
all_Profiles <- complete_profiles_freq_dist_weighted(MIS, weights)[[1]]
all_coordinates <- embedding_coordinates(all_Profiles, result$bestVariablesPriority)
#check: for the observed profiles, these coordinates are those of the optimal embedding
stopifnot(all(all_coordinates[result$bestRepresentation$profiles + 1, ] ==
                as.matrix(result$bestRepresentation[, c("x", "y")])))
comp_input_all  <- count_comparabilities(all_Profiles)
comp_output_all <- count_comparabilities(all_coordinates)

#units in the vertical strip (K2=0) and (K3=1), i.e. units who can compute compound, but not simple, interest:
#from the table TT above, by answer given to the simple interest question (qk5)
strip_K2_0_K3_1 <- TT[, "1"][c("-99", "-97", "0")]
names(strip_K2_0_K3_1) <- c("refused/missing", "don't know", "wrong answer")


############################################################################
# 5. Figure 2: synthesis of the computations
############################################################################
#The panels are drawn at the size they have in the paper (0.495\textwidth, i.e. 2.97 inches) with a font of
#10 points, so that all of the text is printed at 10 points (the size required by the journal).
fig2_side <- 2.97   #side of each panel (inches)
c.p <- 0.5          #size of the points
fig2_pdf <- function(file) {
  pdf(file, width = fig2_side, height = fig2_side, pointsize = 10)
  par(mar = c(2.9, 3.1, 2.6, 0.6), mgp = c(1.9, 0.6, 0), tcl = -0.35)
}

#top-left panel: sorted global losses of the 360 admissible embeddings
fig2_pdf("figures/Figure_2_top_left_Global_loss.pdf")
plot(sort(result$allLoss),
     main = "Global loss", xlab = "2-D map index", ylab = "Loss",
     ylim = c(0, max(result$allLoss)), pch = 20, cex = c.p, cex.main = 1, font.main = 1)
dev.off()

#top-right panel: sorted profile approximation errors at the optimal embedding
fig2_pdf("figures/Figure_2_top_right_Profile_approximation_error.pdf")
plot(sort(result$bestRepresentation[,5]),
     main = "Profile's approximation error", xlab = "Unit index", ylab = "Approximation error",
     ylim = c(0, max(result$bestRepresentation[,5])), pch = 20, cex = c.p, cex.main = 1, font.main = 1)
dev.off()

#bottom-left panel: boxplot of the approximation errors of the statistical units (each profile error is
#repeated as many times as the number of units it represents, ceiling(N*frequency))
npoints <- ceiling(N*result$bestRepresentation[,4])
errorep <- rep(result$bestRepresentation[,5], npoints)
fig2_pdf("figures/Figure_2_bottom_left_Units_approximation_error.pdf")
boxplot(errorep, boxwex = .2, range = 0, col = "white", ylab = "Approximation error")
title(main = "Units' approximation error distribution", font.main = 1, cex.main = 1)
dev.off()

#bottom-right panel: population shares vs approximation errors of the profiles
pop_share <- npoints/sum(npoints)
fig2_pdf("figures/Figure_2_bottom_right_Error_vs_population_share.pdf")
plot(pop_share, result$bestRepresentation[,5],
     main = "Profile's approximation error\nvs Population share", xlab = "Population share", ylab = "Approximation error",
     ylim = c(0, max(result$bestRepresentation[,5])), pch = 20, cex = c.p, cex.main = 1, font.main = 1)
dev.off()


############################################################################
# 6. Figure 3: the optimal embedding
############################################################################
#The figure is drawn at the size it has in the paper (one page, the two maps one above the other), so that
#the sizes of the labels, points and lines are those of the printed figure. The plot area of each panel
#coincides with the 64 x 64 grid (xlim = ylim = c(0, 64), xaxs = yaxs = "i"), whose side is fig3_side
#inches in both panels; the profile legend (variable scores and names) and the color legend are drawn in
#the margins around the grid (xpd = NA), whose sizes are fixed in inches: in this way the two grids are
#squares of exactly the same size and in the same horizontal position.
fig3_total_height <- 7.35    #height of the figure (inches): the largest that fits a page with the caption
fig3_bottom <- 0.53          #margin below each grid (inches): rows of the variable scores on the x-axis
fig3_left   <- 0.53          #margin on the left of each grid (inches): columns of the variable scores on the y-axis
fig3_top    <- 0.35          #margin above each grid (inches): names of the variables of the y-axis
fig3_right  <- 1.00          #margin on the right of each grid (inches): color legend of the bottom panel
fig3_gap    <- 0.20          #extra space between the two panels (inches)
fig3_side   <- (fig3_total_height - 2*(fig3_bottom + fig3_top) - fig3_gap)/2   #side of each grid (inches)
fig3_panel_height <- fig3_bottom + fig3_side + fig3_top              #height of each panel, without the gap (inches)
fig3_width  <- fig3_left + fig3_side + fig3_right                    #width of the figure (inches)
fig3_unit   <- fig3_side/2^k                                         #inches per grid unit
fig3_cex_profile <- 10/12    #size of the variable scores and names on the axes and of the color legend (10 points)
fig3_row_spacing <- 0.14     #distance (inches) between the rows/columns of the profile legend
fig3_profiles_expand <- rep(fig3_row_spacing/fig3_unit/2^(k-3), 2)  #the same spacing, as required by POsetDR.plot
fig3_profiles_lines  <- (fig3_bottom - 0.02)/fig3_unit              #length of the dotted lines of the legend

#top panel: units as jittered clouds of points around the profiles (red points). The clouds are random:
#the seed makes the figure reproducible. N = 3000 is the number of points drawn (profile i is represented
#by round(3000*frequency_i) points); it only affects the density of the clouds.
plot_Figure_3_top <- function(margin_right = fig3_right, gap = 0) {
  set.seed(2026)
  #gap: extra space below the panel (inches), when it is drawn above the bottom panel
  op <- par(mai = c(fig3_bottom + gap, fig3_left, fig3_top, margin_right), xaxs = "i", yaxs = "i", xpd = NA,
            mgp = c(3, 0.2, 0))   #mgp: distance of the variable scores from their rows (in lines)
  POsetDR.plot(result, N = 3000, k = k, gV_x = 3, gV_y = 3,
               Grid = TRUE, Grid_lty = c(1,5,2), Grid_lwd = c(1.1,1.1,1.1), Grid_col = "darkgray",
               Profile = TRUE, Profiles_col = "black", Profiles_expand = fig3_profiles_expand, cex.profile = fig3_cex_profile,
               Profiles_lines = fig3_profiles_lines, varNames = varNames, title = "", type = "jitter_circle",
               xlim = c(0, 2^k), ylim = c(0, 2^k),
               points_col = "black", center_col = "red", points_pch = 16, center_pch = 16, points_cex = 0.22, center_cex = 0.4,
               jitter_amount = 2)
  par(op)
}
#bottom panel: profiles as bubbles (area proportional to frequency) colored by their approximation error,
#with the color legend on the right
plot_Figure_3_bottom <- function() {
  op <- par(mai = c(fig3_bottom, fig3_left, fig3_top, fig3_right), xaxs = "i", yaxs = "i", xpd = NA,
            mgp = c(3, 0.2, 0))
  POsetDR.plot(result, N = N, k = k, gV_x = 3, gV_y = 3,
               Grid = TRUE, Grid_lty = c(1,5,2), Grid_lwd = c(1.1,1.1,1.1), Grid_col = "darkgray",
               Profile = TRUE, Profiles_col = "black", Profiles_expand = fig3_profiles_expand, cex.profile = fig3_cex_profile,
               Profiles_lines = fig3_profiles_lines, varNames = varNames, title = "", type = "bubbles_and_approximation",
               xlim = c(0, 2^k), ylim = c(0, 2^k),
               points_col = "black", center_col = "white", rampPalette = "viridis", center_pch = 16, points_cex = 0.05, center_cex = 0.08,
               expand_radius = 30, min_radius = 0.4, lwd_circle = 0.2, cex.legend = fig3_cex_profile)
  par(op)
}
#single panels
pdf("figures/Figure_3_top_jittered_units.pdf", width = fig3_left + fig3_side + 0.1, height = fig3_panel_height)
plot_Figure_3_top(margin_right = 0.1)
dev.off()
pdf("figures/Figure_3_bottom_bubbles_and_errors.pdf", width = fig3_width, height = fig3_panel_height)
plot_Figure_3_bottom()
dev.off()
#Figure 3 of the paper: the two panels one above the other, separated by fig3_gap inches
pdf("figures/Figure_3.pdf", width = fig3_width, height = fig3_total_height)
layout(matrix(1:2, ncol = 1), heights = c(fig3_panel_height + fig3_gap, fig3_panel_height))
plot_Figure_3_top(gap = fig3_gap)
plot_Figure_3_bottom()
dev.off()


############################################################################
# 7. Figure 4: relative contingencies
############################################################################
#matrix of the 64 profiles with observed and expected (under independence) frequencies, relative
#contingencies and cells in the optimal embedding
complete_matrix_profile <- complete_matrix_profiles_function(MIS, weights, result$bestVariablesPriority,
                                                             c("K_1","K_2","K_3","K_4","K_5","K_6"))
#square grid (asp = 1); device width proportional to the plot area (x range 112 units, y range 100 units).
#The figure is reduced in the paper (0.75\textwidth, i.e. a scale of 0.575): cex.profile = 1.46 gives a text of
#18 points on the device (R rounds the font size to an integer), printed at 10.3 points.
#signs = TRUE: the signs in the cells keep positive and negative values distinct also in black and white, where
#the two ends of the diverging color scale are both dark.
pdf("figures/Figure_4_relative_contingencies.pdf", width = 7*112/100, height = 7)
op <- par(mar = c(0,0,0,0))
Contingency.plot(complete_matrix_profile, result$bestVariablesPriority, k = k, gV = 3, Grid = TRUE,
                 varNames = varNames, Grid_col = "black", Grid_lty = c(1,5,2,4,3)[1:3], Grid_lwd = c(3,2.5,2,1.5,1)[1:3],
                 Profile = TRUE, Profiles_col = "black", Profiles_expand = 1, title = NULL,
                 expandX = c(0,0), expandY = c(0,0), rampPalette1 = c("white", "darkred"),
                 rampPalette2 = c("white", "darkgreen"), cex.profile = 1.46, signs = TRUE)
par(op)
dev.off()


############################################################################
# 8. Figure 5: sub-samples by education degree
############################################################################
#the optimal embedding of the whole sample, with the frequencies of each sub-sample (column of B_edu[[5]];
#the frequencies are not rescaled within the sub-sample, so expand_radius differs among the panels. For the
#class "Elementary or no degree", expand_radius = 49 gives bubbles on the same scale as the panel of the
#units with no degree alone, drawn with expand_radius = 135: 135*sqrt(0.0159/0.1207) = 49, where 0.0159
#and 0.1207 are the population shares of the two classes)
result_edu <- lapply(1:length(edu_labels), function(j) {
  r <- result
  r$bestRepresentation$weights <- B_edu[[5]][, j]
  r
})
names(result_edu) <- edu_labels

#panels of Figure 5
panels_Figure_5 <- data.frame(label = c("High school", "Bachelor/master or higher", "Elementary or no degree", "Middle school"),
                              expand_radius = c(34, 48, 49, 55),
                              position = c("top_left", "top_right", "bottom_left", "bottom_right"))
plot_Figure_5_panel <- function(label, expand_radius) {
  POsetDR.plot(result_edu[[label]], N = N, k = k, gV_x = 3, gV_y = 3,
               Grid = TRUE, Grid_col = "darkgray", Grid_lty = c(1,5,2,4,3)[1:3], Grid_lwd = c(2,2,2,1,1)[1:3],
               Profile = FALSE, Profiles_col = "black", varNames = varNames, Profiles_expand = c(0.5,0.5),
               title = label, type = "bubbles", expandY = c(3,3), expandX = c(3,3),
               points_col = "black", center_col = "black", bubbles_col = "gray92", center_pch = 16, points_cex = 0.3,
               center_cex = 0.4, expand_radius = expand_radius, cex.profile = 1, cex.main = 2.2)
}
#Narrow margins, so that the grids fill the figure: each panel is a square plot area of side fig5_side
#inches, with fig5_mai = (bottom, left, top, right) margins in inches (the top one for the title); the
#bubbles on the border may extend slightly outside the plot area (xpd = NA)
fig5_side <- 5
fig5_mai <- c(0.05, 0.05, 0.55, 0.05)
fig5_width <- fig5_side + fig5_mai[2] + fig5_mai[4]
fig5_height <- fig5_side + fig5_mai[1] + fig5_mai[3]
for (p in 1:nrow(panels_Figure_5)) {
  pdf(paste0("figures/Figure_5_", panels_Figure_5$position[p], "_", gsub("[ /]", "_", panels_Figure_5$label[p]), ".pdf"),
      width = fig5_width, height = fig5_height)
  par(mai = fig5_mai, xpd = NA)
  plot_Figure_5_panel(panels_Figure_5$label[p], panels_Figure_5$expand_radius[p])
  dev.off()
}
#the four panels together
pdf("figures/Figure_5.pdf", width = 2*fig5_width, height = 2*fig5_height)
par(mfrow = c(2, 2), mai = fig5_mai, xpd = NA)
for (p in 1:nrow(panels_Figure_5)) plot_Figure_5_panel(panels_Figure_5$label[p], panels_Figure_5$expand_radius[p])
dev.off()


############################################################################
# 9. Saving of tables, numbers and workspace
############################################################################

#Table 1 (LaTeX and csv)
write.csv(data.frame(number_of_1s = 0:k, relative_frequency = as.numeric(Table_1)), "tables/Table_1.csv", row.names = FALSE)
writeLines(c("\\begin{table}[ht]",
             "\t\\caption{Distribution of the number of entries equal to 1, in the knowledge profiles.}",
             paste0("\t\\centering\\begin{tabular}{r|", strrep("c", k+1), "}"), "\t\t\\hline",
             paste0("\t\tNumber of 1s & ", paste(0:k, collapse = " & "), " \\\\ [0ex]"), "\t\t\\hline",
             paste0("\t\tRelative frequency & ", paste(sprintf("%.2f", Table_1), collapse = " & "), " \\\\ [0ex]"), "\t\t\\hline",
             "\t\\end{tabular}", "\t\\label{Tab.1}", "\\end{table}"), "tables/Table_1.tex")

#numbers quoted in the text
numbers <- c(
  paste("Number of statistical units:", N),
  paste("Number of observed profiles:", nrow(obs_Profiles)),
  paste("Optimal variable index permutation:", paste(result$bestVariablesPriority, collapse = ",")),
  sprintf("Global loss of the optimal embedding: %.4f", result$bestLossValue),
  sprintf("Observed profiles (%d pairs): input lattice %d comparabilities and %d incomparabilities; output poset %d and %d",
          sum(comp_input_obs), comp_input_obs[1], comp_input_obs[2], comp_output_obs[1], comp_output_obs[2]),
  sprintf("All of the %d profiles (%d pairs): input lattice %d comparabilities and %d incomparabilities; output poset %d and %d",
          2^k, sum(comp_input_all), comp_input_all[1], comp_input_all[2], comp_output_all[1], comp_output_all[2]),
  sprintf("Observed profiles: the output poset has %.0f%% comparabilities more and %.0f%% incomparabilities less than the input lattice",
          100*(comp_output_obs[1]/comp_input_obs[1]-1), 100*(1-comp_output_obs[2]/comp_input_obs[2])),
  paste0("Units with K2=0 and K3=1 (compound but not simple interest): ", sum(strip_K2_0_K3_1), " (",
         paste(names(strip_K2_0_K3_1), strip_K2_0_K3_1, sep = ": ", collapse = "; "), " to the simple interest question)"),
  paste("Units by education degree:", paste(names(table(edu)), table(edu), sep = ": ", collapse = "; ")),
  paste("Table 1 (relative frequencies of 0,...,6 entries equal to 1):", paste(sprintf("%.2f", Table_1), collapse = " "))
)
cat(numbers, sep = "\n")
writeLines(numbers, "numbers_quoted_in_the_text.txt")

#R version and packages used
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#workspace (load("Financial_literacy_JCGS.RData") to work on the results without running the script again)
save.image("Financial_literacy_JCGS.RData")
