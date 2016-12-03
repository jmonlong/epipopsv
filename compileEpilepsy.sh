for rmd in src/epilepsy*.Rmd
do
    Rscript -e "rmarkdown::render('$rmd')"
    mv ${rmd%Rmd}md reports/
    cp -r ${rmd%.Rmd}_files reports/
    rm -rf ${rmd%.Rmd}_files
done
