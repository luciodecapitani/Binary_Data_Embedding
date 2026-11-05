OBE vs POSAC: code reproducing the supplementary material
=========================================================

This folder contains the R code and the results of the comparison between the Optimal
Bidimensional Embedding (OBE) and POSAC (Partial Order Scalogram Analysis with base
Coordinates), reported in the supplementary material to

  M. Fattore, A. Avellone, L. De Capitani,
  "Bidimensional Embedding and Visualization of Multidimensional Ordinal Binary Data".

The comparison uses the Financial knowledge data of the Bank of Italy 2020 survey
(D'Alessio et al., 2020), with frequencies obtained by sample weights.


FILES
-----
  OBE_POSAC_comparison.r                main script: it reproduces all of the results
  Preliminary_Functions_OBE_vs_POSAC.r  packages (with their installation, if needed) and
                                        functions used by the main script; it is sourced at
                                        the beginning of OBE_POSAC_comparison.r
  iacofi.csv                            data (Bank of Italy 2020 survey on financial literacy)

  Results, already produced by running the script (see "RESULTS" below):
  OBE_POSAC_comparison.RData            R workspace saved at the end of the script
  figures/                              Figures S1 and S2 (pdf)
  tables/                               Tables S1 and S2 (LaTeX and csv)
  sessionInfo.txt                       R version and packages used to produce the results


REQUIREMENTS
------------
  - R (the results were obtained with R 4.6.2 under Windows 11).
  - R packages:
      poseticDataAnalysis (CRAN)        Optimal Bidimensional Embedding and poset tools
      viridis (CRAN)                    color scales
      POSAC (GitHub, not on CRAN)       https://github.com/CambridgeAssessmentResearch/POSAC
      remotes (CRAN)                    only to install POSAC from GitHub

  Preliminary_Functions_OBE_vs_POSAC.r installs the missing packages automatically. In
  particular, POSAC is installed from GitHub at the commit used for the supplementary material
  (version 1.1.0, commit 6aaeecfc76e643bde973c613bcdef14588e0431a) if it is not installed or if a
  different version is installed. Installing a package from GitHub requires an internet
  connection; on Windows, Rtools may be needed if R asks to compile the package from source.
  To install POSAC by hand, run in R:

      install.packages("remotes")
      remotes::install_github("CambridgeAssessmentResearch/POSAC@6aaeecfc76e643bde973c613bcdef14588e0431a")


HOW TO RUN THE CODE
-------------------
  1. Download the whole folder and keep the files together.
  2. Start R and set the working directory to this folder, e.g.

         setwd("path/to/Code_OBE_vs_POSAC")

  3. Run the main script:

         source("OBE_POSAC_comparison.r", echo = TRUE)

     or, from a terminal, in this folder:

         Rscript OBE_POSAC_comparison.r

  The Optimal Bidimensional Embedding takes about one second; POSAC with the greedy refinement
  takes several minutes (11-13 minutes on an Intel Core i7-6700 CPU at 3.40 GHz).
  POSAC is deterministic (it starts from four fixed configurations), so the results do not depend
  on random numbers.

  Running the script overwrites the files in figures/ and tables/, sessionInfo.txt and
  OBE_POSAC_comparison.RData.


WORKING ON THE RESULTS WITHOUT RUNNING THE SCRIPT
-------------------------------------------------
  To avoid running POSAC again, load the saved workspace (working directory set to this folder):

      source("Preliminary_Functions_OBE_vs_POSAC.r")   # packages and functions (e.g. for the plots)
      load("OBE_POSAC_comparison.RData")

  Main objects of the workspace:
      obe        output of OptimalBidimensionalEmbedding (OBE map)
      posac      output of POSAC (POSAC map)
      obe_names  profile names (binary strings) of the rows of obe$bestRepresentation
      B          observed profiles (B[[1]]), sample counts (B[[2]]), frequencies obtained by
                 sample weights (B[[3]])
      named_w    frequencies obtained by sample weights, named by profile
      M_t_obs    incidence matrix of the observed poset
      M_a_obe, M_a_posac   incidence matrices of the posets induced by the two maps
      Table_S1, Table_S2   Tables S1 and S2 of the supplementary material
      Table_no_context     incomparable pairs by difference in the number of 1s of the two
                           profiles (paragraph "No context" of the supplementary material)
      example_pair         frequencies and map positions of the two profiles of the example in
                           that paragraph
      D_obs      degrees of dominance of the observed profiles


RESULTS
-------
  tables/Table_S1_indicators.tex/.csv   Table S1: indicators based on characteristic matrices
                                        (shares of correctly represented pairs, unweighted and
                                        weighted by the product of the profile frequencies)
  tables/Table_S2_pairs.tex/.csv        Table S2: classification of the pairs of observed
                                        profiles (preserved, turned into incomparable, reversed,
                                        turned into comparable)
  figures/Figure_S1_OBE_POSAC_maps.pdf  Figure S1: POSAC map and OBE map (bubbles with area
                                        proportional to the profile frequencies)
  figures/Figure_S2_POSAC_by_variable.pdf
                                        Figure S2: POSAC map, one panel per variable (blue: the
                                        variable scores 0; red: it scores 1)
  Both figures are drawn at their printed size (6.5 inches wide), with text of 10 points or more.

  The numbers quoted in the text of the supplementary material (number of pairs, optimal
  variable permutation of OBE, shares of correctly represented pairs, comparabilities broken by
  POSAC, incomparabilities turned into comparabilities by OBE) are printed at the end of the
  script. The same holds for the numbers of the paragraph "No context": degrees of dominance and
  shares of incomparable pairs made comparable by each map, by difference in the number of 1s of
  the two profiles; position of the profiles 100000 and 011111 in the two maps; comparable pairs
  in the two maps.


REFERENCE FOR THE DATA
----------------------
  D'Alessio, G., De Bonis, R., Neri, A., Rampazzi, C. (2020), "Financial literacy in Italy: the
  results of the Bank of Italy's 2020 survey", Bank of Italy, Occasional Papers 588.
