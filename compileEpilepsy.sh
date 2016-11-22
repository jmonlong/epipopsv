for rmd in src/epilepsy*.Rmd
do
    Rscript -e "rmarkdown::render('$rmd', output_dir='reports/')"
done
