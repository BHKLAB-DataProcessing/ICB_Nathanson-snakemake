# wget -O ~/Documents/ICBCuration/data_source/Nathanson/cohort.csv https://raw.githubusercontent.com/hammerlab/melanoma-reanalysis/master/files/cohort.csv
# wget -O ~/Documents/ICBCuration/data_source/Nathanson/cufflinks.tar.gz https://github.com/hammerlab/melanoma-reanalysis/raw/master/files/cufflinks.tar.gz
# wget -O ~/Documents/ICBCuration/data_source/Nathanson/snvs.csv https://raw.githubusercontent.com/hammerlab/melanoma-reanalysis/master/files/snvs.csv

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

# CLIN.txt
clin <- read.csv(file.path(work_dir, 'cohort.csv'))
write.table(clin, file=file.path(work_dir, 'CLIN.txt'), quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE)

# EXPR.csv.gz
untar(file.path(work_dir, 'cufflinks.tar.gz'), exdir=file.path(work_dir))
expr <- read.csv(file.path(work_dir, 'cufflinks.csv'))
gz <- gzfile(file.path(work_dir, 'EXPR.csv.gz'), "w")
write.table( expr , file=gz , quote=FALSE , sep=',' , col.names=TRUE , row.names=FALSE )
close(gz)
file.remove(file.path(work_dir, 'cufflinks.csv'))

# SNV.txt.gz
snv <- read.csv(file.path(work_dir, 'snvs.csv'))
gz <- gzfile(file.path(work_dir, 'SNV.txt.gz'), "w")
write.table( snv , file=gz , quote=FALSE , sep='\t' , col.names=TRUE , row.names=FALSE )
close(gz)