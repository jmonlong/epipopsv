for rmd in src/epilepsy*.Rmd
do
    Rscript -e "rmarkdown::render('$rmd')"
    mv -f ${rmd%Rmd}md reports/
    cp -r ${rmd%.Rmd}_files reports/
    rm -rf ${rmd%.Rmd}_files
done
