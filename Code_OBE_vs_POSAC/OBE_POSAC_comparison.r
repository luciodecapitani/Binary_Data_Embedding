############################################################################
############################################################################
# OBE_POSAC_comparison.r
#
# Comparison between the Optimal Bidimensional Embedding (OBE) and POSAC (Partial Order Scalogram
# Analysis with base Coordinates) on the Financial knowledge data of the Bank of Italy 2020 survey
# (file iacofi.csv), with frequencies obtained by sample weights.
# The script reproduces all of the results of the (short) supplementary material to
# "Bidimensional Embedding and Visualization of Multidimensional Ordinal Binary Data":
#   - Table S1  : indicators based on characteristic matrices (tables/Table_S1_indicators.tex/.csv)
#   - Table S2  : classification of the pairs of observed profiles (tables/Table_S2_pairs.tex/.csv)
#   - Figure S1 : OBE and POSAC maps (figures/Figure_S1_OBE_POSAC_maps.pdf)
#   - Figure S2 : POSAC map, one panel per variable (figures/Figure_S2_POSAC_by_variable.pdf)
#   - the numbers quoted in the text (printed at the end of the script)
# The whole workspace is saved in OBE_POSAC_comparison.RData.
#
# STRUCTURE OF THE PROGRAM
#   0. Packages and functions
#   1. Loading and preparing the data
#   2. Incidence matrix of the Boolean lattice B_6
#   3. Optimal Bidimensional Embedding
#   4. POSAC
#   5. Comparison of the two maps (Tables S1 and S2)
#   6. Figures S1 and S2
#   7. Numbers quoted in the text, saving of tables and workspace
#
# Run the script with the working directory set to the folder containing this file, the data file
# iacofi.csv and Preliminary_Functions_OBE_vs_POSAC.r, e.g. setwd("path/to/Code_OBE_vs_POSAC").
# NB: POSAC with the greedy refinement takes several minutes (11-13 minutes on an Intel Core i7-6700).
############################################################################
############################################################################


############################################################################
# 0. Packages and functions
############################################################################
# installs (only if needed) and loads poseticDataAnalysis, POSAC (pinned GitHub commit) and viridis,
# and defines the functions used below
source("Preliminary_Functions_OBE_vs_POSAC.r")

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)


############################################################################
# 1. Loading and preparing the data
############################################################################

DT <- read.csv("iacofi.csv")

#selecting the variables of interest
SUBDT <- DT[,c("qk4", "qk5", "qk6", "qk7_1","qk7_1alt", "qk7_2", "qk7_3", "qk7_3alt","qd1","qd7","area5","reddito","edu","capofam")]

#union of the answers to questions qk7_1 and qk7_1alt, and to qk7_3 and qk7_3alt (the same questions,
#asked in two versions)
SUBDT[which(is.na(SUBDT[,7])),7] <- SUBDT[which(is.na(SUBDT[,7])),8]
SUBDT[which(is.na(SUBDT[,4])),4] <- SUBDT[which(is.na(SUBDT[,4])),5]
SUBDT$qk7_1alt <- SUBDT$qk7_3alt <- NULL
#coding of the correct answers: qk4: 0; qk5: 102; qk6: 1; qk7_1 (merged): 1; qk7_2: 1; qk7_3 (merged): 1
corrette <- c(0,102,1,1,1,1)

#sample weights (normalized to sum 1) and sample size
weights <- DT$wght
weights <- weights/sum(weights)
N <- nrow(SUBDT)

#binary data matrix: 1 = correct answer, 0 = wrong/missing answer
MIS <- as.data.frame(matrix(-2, nrow=N, ncol=6))
colnames(MIS) <- c("qk4", "qk5", "qk6", "qk7_1", "qk7_2", "qk7_3")
for(h in 1:nrow(SUBDT)){
  MIS[h,] <- (SUBDT[h,1:6]==corrette)
}

#number of variables and their names in the paper (K1,...,K6 correspond to qk4,...,qk7_3)
k <- ncol(MIS)
varNames <- c("K1","K2","K3","K4","K5","K6")   #plain text, as in the paper: subscripts would be printed below 10 points

#frequency distribution over the observed profiles:
#B[[1]] observed profiles, B[[2]] sample counts, B[[3]] relative frequencies obtained by sample weights
B <- profiles_freq_dist_weighted(MIS, weights)
obs_Profiles <- B[[1]]
#profile names (binary strings, as in the POSAC output) and named vector of the frequencies
profile <- apply(obs_Profiles, 1, paste, collapse = "")
named_w <- B[[3]]
names(named_w) <- profile


############################################################################
# 2. Incidence matrix of the Boolean lattice B_6
############################################################################

pp <- LinearPOSet(c("0","1"))
PPp <- ProductPOSet(pp,pp,pp,pp,pp,pp)
M_t <- DominanceMatrix(PPp)
#rows and columns renamed consistently with the profile names of POSAC (e.g. "0_1_0_0_1_1" -> "010011")
colnames(M_t) <- gsub("_", "", colnames(M_t))
rownames(M_t) <- gsub("_", "", rownames(M_t))


############################################################################
# 3. Optimal Bidimensional Embedding
############################################################################
# exhaustive search over the k!/2 = 360 admissible embeddings (reversed pairs of lexicographic linear
# extensions), with the relative frequencies obtained by sample weights

time_OBE <- system.time(
  obe <- OptimalBidimensionalEmbedding(obs_Profiles, B[[3]], output_every_sec = 1)
)
obe$bestVariablesPriority    # optimal variable index permutation
#profile names of the OBE output (the profile identifiers are decimal numbers)
obe_names <- sapply(obe$bestRepresentation[, 1], function(x) int2bin(x, n_bits = k))


############################################################################
# 4. POSAC
############################################################################
# POSAC is run on the frequencies obtained by sample weights rescaled to the sample size (summing to N),
# which is the scale the package is designed for (frequencies of units), with the greedy refinement.
# NB: the POSAC solution depends on the scale of the frequencies (its optimizer is not scale invariant),
# while all of the percentages computed below do not.

time_POSAC <- system.time(
  posac <- POSAC(obs_Profiles, B[[3]]*N, greedy = TRUE)
)
#check that the profiles of the two outputs are the same
stopifnot(setequal(obe_names, posac$Patterns))


############################################################################
# 5. Comparison of the two maps (Tables S1 and S2)
############################################################################

#incidence matrices of the posets induced by the two maps on the observed profiles
M_a_posac <- intersection_incidence(posac$Patterns, posac$X, posac$Y)
M_a_obe <- intersection_incidence(obe_names, obe$bestRepresentation[, 2], obe$bestRepresentation[, 3])
M_a_obe <- M_a_obe[rownames(M_a_posac), colnames(M_a_posac)]
#incidence matrix of the observed poset (Boolean lattice restricted to the observed profiles)
M_t_obs <- M_t[rownames(M_a_posac), colnames(M_a_posac)]
w_eval <- named_w[rownames(M_a_posac)]

#Table S1: pairs counted once ("unweighted") or weighted by the product of the profile frequencies ("weighted")
Table_S1 <- rbind(compare_incidence(M_t_obs, M_a_posac, freqs = NULL),
                  compare_incidence(M_t_obs, M_a_obe,   freqs = NULL),
                  compare_incidence(M_t_obs, M_a_posac, freqs = w_eval),
                  compare_incidence(M_t_obs, M_a_obe,   freqs = w_eval))
Table_S1 <- data.frame(pairs = c("unweighted", "unweighted", "weighted", "weighted"),
                       method = c("POSAC", "OBE", "POSAC", "OBE"), Table_S1)
print(Table_S1, digits = 4)

#Table S2: classification of the pairs of distinct observed profiles
Table_S2 <- rbind(POSAC = pair_classification(M_t_obs, M_a_posac, w_eval),
                  OBE   = pair_classification(M_t_obs, M_a_obe,   w_eval))
print(round(Table_S2, 2))


############################################################################
# 6. Figures S1 and S2
############################################################################

#Both figures are drawn at the size they have in the supplementary material (6.5 inches wide), so that the
#text is printed at 10 points or more. In each panel the plot area coincides with the square map, and the
#device height is computed so that the maps are squares (of the same size within each figure).
fig_width <- 6.5
fig_cex   <- 10/12    #size of the text (10 points, with the device font of 12 points)
n_obs     <- length(posac$X)

#Figure S1: POSAC map (left) and OBE map (right), bubbles with area proportional to profile frequencies
S1_mai  <- c(0.53, 0.55, 0.38, 0.28)                 #margins of each panel (inches): bottom, left, top, right
S1_side <- fig_width/2 - S1_mai[2] - S1_mai[4]       #side of the maps (inches)
pdf("figures/Figure_S1_OBE_POSAC_maps.pdf", width = fig_width, height = S1_mai[1] + S1_side + S1_mai[3])
OBE_POSAC_maps(obe, posac, freqs_posac = posac$freqs, N = N, k = k, varNames = varNames,
               mai = S1_mai, side = S1_side, cex.text = fig_cex)
dev.off()

#Figure S2: POSAC map, one panel per variable (blue = score 0, red = score 1), on two rows of three panels
S2_mai    <- c(0.38, 0.42, 0.30, 0.16)               #margins of each panel (inches): bottom, left, top, right
S2_legend <- 0.28                                    #space for the legend at the bottom of the figure (inches)
S2_side   <- fig_width/3 - S2_mai[2] - S2_mai[4]     #side of the maps (inches)
pdf("figures/Figure_S2_POSAC_by_variable.pdf", width = fig_width, height = 2*(S2_mai[1] + S2_side + S2_mai[3]) + S2_legend)
POSAC.itemdiagram.bubbles(posac, freqs = posac$freqs, vars = 1:k, var_names = varNames, col0 = "blue", col1 = "red",
                          panel_mai = S2_mai, legend_height = S2_legend, legend_cex = fig_cex,
                          center_cex = 0.05, center_col = "white", expand_radius = 30*(n_obs+1)/2^k,
                          min_radius = 0.4, lwd_circle = 0.3, cex.profile = fig_cex, axes_by = 20, Frame_lwd = 1,
                          xlim = c(0, n_obs+1), ylim = c(0, n_obs+1), asp = NA)
dev.off()


############################################################################
# 7. Numbers quoted in the text, saving of tables and workspace
############################################################################

n_pairs  <- Table_S1$npairs[1]                       # pairs of distinct observed profiles
n_comp   <- Table_S1$npairscomp[1]                   # comparable pairs of the observed poset
n_incomp <- Table_S1$npairsincomp[1]                 # incomparable pairs of the observed poset

cat("\nNumber of observed profiles:", nrow(obs_Profiles), "\n")
cat("Pairs of distinct observed profiles:", n_pairs, "(comparable:", n_comp, ", incomparable:", n_incomp, ")\n")
cat("OBE optimal variable index permutation:", obe$bestVariablesPriority, "\n")
cat(sprintf("Weighted share of correctly represented pairs: POSAC %.2f%%, OBE %.2f%%\n",
            Table_S1$pccorrect[3], Table_S1$pccorrect[4]))
cat("Input comparabilities broken by POSAC:", Table_S2["POSAC", "comp_to_incomp"], "of", n_comp, "\n")
cat("Incomparabilities turned into comparabilities by OBE:", Table_S2["OBE", "incomp_to_comp"], "\n")
cat(sprintf("Computing time: OBE %.3f s, POSAC (with greedy refinement) %.1f min\n",
            time_OBE[["elapsed"]], time_POSAC[["elapsed"]]/60))

#Numbers quoted in the paragraph "No context".
#Degrees of dominance of the observed profiles (they depend only on the Boolean lattice B_6, not on the data)
D_obs <- BLSDominance(PPp)
dimnames(D_obs) <- lapply(dimnames(D_obs), function(x) gsub("_", "", x))
D_obs <- D_obs[rownames(M_t_obs), colnames(M_t_obs)]
#incomparable pairs by difference in the number of 1s of the two profiles
Table_no_context <- incomparable_by_difference(M_t_obs, list(OBE = M_a_obe, POSAC = M_a_posac), w_eval, D_obs)
cat("\nIncomparable pairs by difference in the number of 1s of the two profiles: degree of dominance, number of\n",
    "pairs and, for each map, pairs made comparable (_to_comp) and their percentage, counting each pair once\n",
    "(_pc) or weighting it by the product of the frequencies of its profiles (_w_pc):\n", sep = "")
print(round(Table_no_context, 3), row.names = FALSE)
#the pair of profiles of the example: frequencies, degree of dominance and position in the two maps
example_profiles <- c("100000", "011111")
i_obe <- match(example_profiles, obe_names)
i_posac <- match(example_profiles, posac$Patterns)
example_pair <- data.frame(profile = example_profiles, freq = as.vector(named_w[example_profiles]),
                           OBE_x = obe$bestRepresentation[i_obe, 2], OBE_y = obe$bestRepresentation[i_obe, 3],
                           POSAC_x = posac$X[i_posac], POSAC_y = posac$Y[i_posac])
cat("\nProfiles of the example (degree of dominance of the pair:",
    round(max(D_obs[example_profiles[1], example_profiles[2]], D_obs[example_profiles[2], example_profiles[1]]), 3), "):\n")
print(example_pair, digits = 3, row.names = FALSE)
#incomparable pairs made comparable by OBE, and comparable pairs in the two maps
cat(sprintf("\nIncomparable pairs made comparable by OBE: %.2f%% of the pairs, %.2f%% of their weight\n",
            100*Table_S2["OBE", "incomp_to_comp"]/n_incomp, Table_S2["OBE", "w_pc_incomp_to_comp"]))
n_comp_map <- Table_S2[, "comp_preserved"] + Table_S2[, "comp_reversed"] + Table_S2[, "incomp_to_comp"]
cat("Comparable pairs in the maps: POSAC", n_comp_map["POSAC"], ", OBE", n_comp_map["OBE"], "\n")

#tables (LaTeX and csv)
writeLines(latex_table_S1(Table_S1), "tables/Table_S1_indicators.tex")
writeLines(latex_table_S2(Table_S2, n_comp, n_incomp), "tables/Table_S2_pairs.tex")
write.csv(Table_S1, "tables/Table_S1_indicators.csv", row.names = FALSE)
write.csv(data.frame(method = rownames(Table_S2), Table_S2), "tables/Table_S2_pairs.csv", row.names = FALSE)

#R version and packages used
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#workspace (load("OBE_POSAC_comparison.RData") to work on the results without running POSAC again)
save.image("OBE_POSAC_comparison.RData")
