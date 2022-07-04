library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

snv = as.data.frame( fread(  file.path(input_dir, "SNV.txt.gz") , stringsAsFactors=FALSE , sep="\t" ))
snv = sort( unique( snv[ , "sample" ] ) ) 

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )
rna = as.data.frame( fread( file.path(input_dir, "EXPR.csv.gz") , stringsAsFactors=FALSE  , sep="," ) )
rna = sort( unique( rna$sample))

patient = sort( unique( c( intersect( clin$sample , rna ) ,  intersect( clin$sample , snv ) ) ) )

case = as.data.frame( cbind( patient , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) ) )
colnames(case) = c( "patient" , "snv" , "cna" , "expr" )
rownames(case) = patient

case$snv = as.numeric( as.character( case$snv ) )
case$cna = as.numeric( as.character( case$cna ) )
case$expr = as.numeric( as.character( case$expr ) )

for( i in 1:nrow(case)){
  if( rownames(case)[i] %in% snv ){
    case$snv[i] = 1
  }
  if( rownames(case)[i] %in% rna ){
    case$expr[i] = 1
  }
}

case$patient = paste( "P" , case$patient , sep="" )

write.table( case , file=file.path(output_dir, "cased_sequenced.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )


