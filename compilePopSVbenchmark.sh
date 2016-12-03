for rmd in src/PopSV*.Rmd src/*wgs-bias.Rmd
do
    Rscript -e "rmarkdown::render('$rmd')"
    mv -f ${rmd%Rmd}md reports/
    mv  -f ${rmd%.Rmd}_files reports/
done
