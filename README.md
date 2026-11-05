# Replication material

R code and results for

> M. Fattore, A. Avellone, L. De Capitani, *Bidimensional Embedding and Visualization of Multidimensional Ordinal Binary Data* (manuscript submitted to the *Journal of Computational and Graphical Statistics*).

The repository has two independent folders. Each one contains everything it needs, including its own copy of the data, and has its own `README.txt` with detailed instructions.

| Folder | What it reproduces | Main script | Running time |
|---|---|---|---|
| [`Code_Financial_literacy`](Code_Financial_literacy) | The Financial knowledge example of the paper: Table 1, Figures 2–5 and the numbers quoted in the text | `Financial_literacy_JCGS.r` | Less than a minute |
| [`Code_OBE_vs_POSAC`](Code_OBE_vs_POSAC) | The comparison between the Optimal Bidimensional Embedding (OBE) and POSAC in the Supplementary Material: Tables S1–S2, Figures S1–S2 and the numbers quoted in the text | `OBE_POSAC_comparison.r` | About 15 minutes |

Figures 1 and 6–11 and Tables 2–3 of the paper (illustrations of the methodology on small examples and computing times) are not produced by this code.

## How to run the code

1. Download the repository (*Code → Download ZIP*) and extract it.
2. Start R and set the working directory to one of the two folders, for example:

   ```r
   setwd("path/to/Code_Financial_literacy")
   ```

3. Run the main script of that folder:

   ```r
   source("Financial_literacy_JCGS.r", echo = TRUE)
   ```

Running a script overwrites the results stored in its folder.

## What each folder contains

- the main R script and a file with the auxiliary functions, both commented;
- the data (`iacofi.csv`);
- `README.txt`, with the requirements, the instructions and the description of the results;
- the results already produced by the script: `figures/`, `tables/`, the R workspace (`.RData`) and `sessionInfo.txt`.

## Requirements

The results were obtained with R 4.6.2 under Windows 11. The required packages are installed automatically if they are missing, which needs an internet connection:

- `poseticDataAnalysis` and `viridis`, from CRAN, in both folders;
- `POSAC`, from [GitHub](https://github.com/CambridgeAssessmentResearch/POSAC) at a fixed commit, with `remotes`, from CRAN, only in `Code_OBE_vs_POSAC`.

## Licence

The code is released under the [MIT License](LICENSE). The licence does not cover the data: see below.

## Data

The file `iacofi.csv` comes from the 2020 edition of the Bank of Italy survey on the financial literacy of adults in Italy (*Indagine sull'alfabetizzazione finanziaria degli adulti in Italia*):

D'Alessio, G., De Bonis, R., Neri, A., Rampazzi, C. (2020), "Financial literacy in Italy: the results of the Bank of Italy's 2020 survey", Bank of Italy, Occasional Papers 588.

**Conditions of use.** The data belong to the Bank of Italy and are subject to its conditions, not to any licence of the code in this repository:

- the Bank of Italy distributes these anonymized elementary data for research purposes;
- whoever uses them does so under their own full and sole responsibility, without involving the Bank of Italy;
- the source must be cited as "Banca d'Italia", together with the name of the dataset and, if available, the version of the archive used;
- commercial use and modifications of the data published on the Bank of Italy website are not allowed.

Before using the data, please refer to the [survey page](https://www.bancaditalia.it/statistiche/tematiche/indagini-famiglie-imprese/alfabetizzazione/index.html), from which they can be downloaded, the [Research Data Center page](https://www.bancaditalia.it/statistiche/basi-dati/rdc/index.html) and the [copyright page](https://www.bancaditalia.it/footer/copyright/) of the Bank of Italy.
