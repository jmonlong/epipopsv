for rmd in src/PopSV*.Rmd src/*wgs-bias.Rmd
do
    Rscript -e "rmarkdown::render('$rmd')"
    mv ${rmd%Rmd}md reports/
    mv ${rmd%.Rmd}_files reports/
done
