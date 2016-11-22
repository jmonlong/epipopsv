# Other scripts

## Running the methods

Find how PopSV and the other methods were run in their respective folder.

## Annotation files

 The `annotation` folder contains scripts to get public annotations (mainly from UCSC server). 

`downloads.sh` will download the necessary data and `dataFormat.R` will format it for easy use in the analysis. Of note, packages `GenomicRanges` and `dplyr` are required for `dataFormat.R`.

Hence the two commands to prepare the annotation files would be : 

```sh
sh downloads.sh
Rscript dataFormat.R
```

