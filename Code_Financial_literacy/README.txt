Financial knowledge example: code reproducing the results of the paper
=======================================================================

This folder contains the R code and the results of the Financial knowledge example of the paper

  M. Fattore, A. Avellone, L. De Capitani,
  "Bidimensional Embedding and Visualization of Multidimensional Ordinal Binary Data",

based on the data of the Bank of Italy 2020 survey on financial literacy (D'Alessio et al., 2020).


FILES
-----
  Financial_literacy_JCGS.r                  main script: it reproduces all of the results
  Preliminary_Functions_Financial_literacy.r packages (with their installation, if needed) and
                                             functions used by the main script (data preparation,
                                             order structure, relative contingencies, plots); it is
                                             sourced at the beginning of Financial_literacy_JCGS.r
  iacofi.csv                                 data (Bank of Italy 2020 survey on financial literacy)

  Results, already produced by running the script (see "RESULTS" below):
  Financial_literacy_JCGS.RData              R workspace saved at the end of the script
  figures/                                   Figures 2-5 of the paper (pdf)
  tables/                                    Table 1 of the paper (LaTeX and csv)
  numbers_quoted_in_the_text.txt             numbers quoted in the text of the paper
  sessionInfo.txt                            R version and packages used to produce the results


REQUIREMENTS
------------
  - R (the results were obtained with R 4.6.2 under Windows 11).
  - R packages (both on CRAN):
      poseticDataAnalysis   Optimal Bidimensional Embedding and poset tools
      viridis               color scale of the approximation errors

  Preliminary_Functions_Financial_literacy.r installs the missing packages automatically
  (an internet connection is needed). To install them by hand, run in R:

      install.packages(c("poseticDataAnalysis", "viridis"))


HOW TO RUN THE CODE
-------------------
  1. Download the whole folder and keep the files together.
  2. Start R and set the working directory to this folder, e.g.

         setwd("path/to/Code_Financial_literacy")

  3. Run the main script:

         source("Financial_literacy_JCGS.r", echo = TRUE)

     or, from a terminal, in this folder:

         Rscript Financial_literacy_JCGS.r

  The whole script takes less than a minute (the exhaustive search over the 360 admissible
  embeddings takes a fraction of a second; most of the time is spent preparing the data).
  The left panel of Figure 3 (jittered clouds of points) is random: the script fixes the seed,
  so that the figure is reproducible.

  Running the script overwrites the files in figures/ and tables/, numbers_quoted_in_the_text.txt,
  sessionInfo.txt and Financial_literacy_JCGS.RData.


WORKING ON THE RESULTS WITHOUT RUNNING THE SCRIPT
-------------------------------------------------
  Load the saved workspace (working directory set to this folder):

      source("Preliminary_Functions_Financial_literacy.r")   # packages and functions (e.g. for the plots)
      load("Financial_literacy_JCGS.RData")

  Main objects of the workspace:
      MIS          binary data matrix (one row per unit: 1 = correct answer)
      weights      sample weights of the units (normalized to sum 1)
      B            observed profiles (B[[1]]), sample counts (B[[2]]), frequencies obtained by
                   sample weights (B[[3]])
      B_edu        as B, with the frequencies split by education degree (B_edu[[4]] counts,
                   B_edu[[5]] frequencies obtained by sample weights, one column per degree)
      result       output of OptimalBidimensionalEmbedding: all of the 360 losses (allLoss), the
                   optimal variable index permutation (bestVariablesPriority), its global loss
                   (bestLossValue) and the optimal embedding (bestRepresentation: profile
                   identifier, coordinates x and y, frequency and approximation error of each
                   observed profile)
      complete_matrix_profile   the 64 profiles with observed and expected (under independence)
                   frequencies and relative contingencies (Figure 4)
      Table_1      Table 1 of the paper


RESULTS
-------
  tables/Table_1.tex/.csv     Table 1: distribution of the number of entries equal to 1 in the
                              knowledge profiles (frequencies obtained by sample weights)
  figures/Figure_2_*.pdf      Figure 2: the four panels (global losses of the 360 admissible
                              embeddings; profile approximation errors; boxplot of the units'
                              approximation errors; approximation errors vs population shares),
                              drawn at their printed size (2.97 x 2.97 inches each)
  figures/Figure_3_top_jittered_units.pdf, Figure_3_bottom_bubbles_and_errors.pdf
                              Figure 3: the optimal embedding, with units as jittered clouds
                              (top) and with profiles as circles colored by their approximation
                              error (bottom); Figure_3.pdf (the figure of the paper, drawn at its
                              printed size, 4.23 x 7.35 inches) puts the two panels one above the
                              other, with grids of the same size and aligned
  figures/Figure_4_relative_contingencies.pdf
                              Figure 4: relative contingencies (observed vs expected under
                              independence); the signs in the cells mark positive and negative
                              values
  figures/Figure_5_*.pdf      Figure 5: the optimal embedding with the frequencies of the
                              sub-samples with different education degrees (one file per panel);
                              Figure_5.pdf puts the four panels together
  numbers_quoted_in_the_text.txt
                              optimal variable index permutation, global loss, numbers of
                              comparabilities and incomparabilities of the input lattice and of
                              the output poset, and the other numbers quoted in the text


REFERENCE FOR THE DATA
----------------------
  D'Alessio, G., De Bonis, R., Neri, A., Rampazzi, C. (2020), "Financial literacy in Italy: the
  results of the Bank of Italy's 2020 survey", Bank of Italy, Occasional Papers 588.
