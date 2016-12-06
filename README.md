This repository hosts the source code and instructions to reproduce the analysis and results from our study on epilepsy and copy number variation.

If you just want the CNV calls, find them in the [FigShare repository](https://figshare.com/s/20dfdedcc4718e465185), or directly [here](https://ndownloader.figshare.com/files/6994223?private_link=20dfdedcc4718e465185).

To review the code and resulting graphs/numbers, have a look at the R-markdown reports in the [`reports`](https://github.com/jmonlong/epipopsv/tree/master/reports) folder and scripts in the  [`src`](https://github.com/jmonlong/epipopsv/tree/master/src) folder.

To rerun the analysis, follow these steps:

# Step 1: Download the relevant data

The necessary data has been deposited on [FigShare](https://figshare.com/s/20dfdedcc4718e465185). Depending on the analysis, you might not need to download all the data.

Still, the easiest way is to download all the data (1.5 GB) and unzip it in the `data` folder.

Soon, we will prepare different packs that will be downloadable with:

+ `downloadBenchmark.sh` (XX Mb) for the *in silico* benchmark of PopSV and comparison with existing methods, using the [Twin study](https://www.ebi.ac.uk/ena/data/view/PRJEB8308) and [CageKid](https://www.ebi.ac.uk/ega/studies/EGAS00001000083) datasets.
+ `downloadEpilepsy.sh` (XX Mb) for the SV analysis in the epilepsy/control cohort.


# Step 2: Install R dependencies

Many different packages are used throughout the analysis. The commands to install them are written in the `installDependencies.R`. To install all the necessary packages open R and run `source("installDependencies.R")`.

# Step 3: Compile the R-markdown reports

The raw R-markdown reports are located in the `src` folder. To recompile them simply run:

```r
library(rmarkdown)
render("XXX.Rmd")
```

You can also compile a bunch of reports using the `compile*.sh` scripts. For example:

```sh
sh compileEpilepsy.sh
```

You can already see the reports produced by these scripts in the `reports` folder. 

# Notes

The code was tested on fresh dockerized Ubuntu with [R 3.2.5](). Windows is not recommended as the `parallel` package is not available.

Results might differ slightly from the ones in our paper because several results are based on permutations or sampling procedures.

The number of permutations for a few analysis have been reduced in order to be able to compute them on a laptop in a reasonable amount of time. For the paper we used high performance computers to increase the number of permutations.

For the most time-consuming steps, we used the *caching* option during the report compilation. It means that it will take the normal time to compile the first time, but will avoid rerunning the long steps on the following compilations.
