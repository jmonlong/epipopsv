for rmd in src/epilepsy*.Rmd
do
    Rscript -e "rmarkdown::render('$rmd')"
    mv -f ${rmd%Rmd}md reports/
    mv  -f ${rmd%.Rmd}_files reports/
done
