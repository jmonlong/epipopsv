This repository hosts the source code and instructions to reproduce the analysis and results from our study on epilepsy and structural variation.

To look at the code and resulting graphs/numbers, have a look at the R-markdown reports in the `reports` folder and scripts in the `src` folder.

To rerun the analysis, follow these steps:

# Step 1: Download the relevant data

The necessary data has been deposited on [FigShare](). Depending on the analysis, you might not need to download all the data.

The easiest way to download the data is to use the scripts `download*.sh`. For example, to download all the data (XX Mb):

```sh
sh downloadAll.sh
```

The different packs are:

+ `downloadBenchmark.sh` (XX Mb) for the *in silico* benchmark of PopSV and comparison with existing methods, using the [Twin study]() and [CageKid]() datasets.
+ `downloadEpilepsy.sh` (XX Mb) for the SV analysis in the epilepsy/control cohort.


# Step 2: Install R dependencies

Many different packages are used throughout the analysis. The commands to install them are written in the `installDependencies.R`. To install all the necessary packages open R and run `source("installDependencies.R")`.

# Step 3: Compile the R-markdown reports

The raw R-markdown reports are located in the `src` folder. The reports produced by these scripts are located in the `reports` folder. To recompile them simply run:

```r
library(rmarkdown)
render("XXX.Rmd")
```

You can also compile a bunch of reports using the `compile*.R` scripts. For example:

```sh
Rscript compileEpilepsy.R
```

# Notes

Results might differ slightly from the ones in our paper because several results are based on permutations or sampling procedures.

The number of permutations for a few analysis have been reduced in order to be able to compute them on a laptop in a reasonable amount of time. For the paper we used high performance computers to increase the number of permutations. 

For the most time-consuming steps, we used the *caching* option during the report compilation. It means that it will take the normal time to compile the first time, but will avoid rerunning the long steps on the following compilations.
