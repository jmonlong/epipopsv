for rmd in src/epilepsy*.Rmd
do
    Rscript -e "rmarkdown::render('$rmd')"
    mv ${rmd%Rmd}md reports/
    mv ${rmd%.Rmd}_files reports/
done
