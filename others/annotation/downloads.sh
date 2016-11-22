wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz

## Gene annotation
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
zcat gencode.v19.annotation.gtf.gz | awk '{if($3=="gene"){print $0}}' | gzip > gencode.v19.annotation.gene.gtf.gz
zcat gencode.v19.annotation.gtf.gz | awk '{if($3=="exon"){print $0}}' | gzip > gencode.v19.annotation.exon.gtf.gz
