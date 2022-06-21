library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

expr = as.data.frame( fread( file.path(input_dir, "EXPR.csv.gz") , stringsAsFactors=FALSE  , sep="," ) )
expr = expr[ expr$FPKM_status %in% "OK" , ]

geneID = sort( unique( expr$gene_short_name ) )
patient = unique(expr$sample)

data = matrix( NA , nrow = length(geneID) , ncol = length(patient) )
colnames(data) = patient
rownames(data) = geneID

for(i in 1:length(patient)){

	e = expr[ expr$sample %in% patient[i] , ]$FPKM
	names(e) = expr[ expr$sample %in% patient[i] , ]$gene_short_name

	data[ names(e) , i ] = e

}

colnames(data) = paste( "P" , colnames(data) , sep="" )

data = as.data.frame( data )
for( i in 1:ncol(data) ){ data[ , i ] = as.numeric( as.character( data[ , i] ) ) }


#################################
#################################

fpkmToTpm <- function(fpkm)
{
    exp( log( fpkm ) - log( sum( fpkm , na.rm = TRUE ) ) + log( 1e6 ) )
}

expr = fpkmToTpm( fpkm = as.matrix( data ) )

#################################
#################################

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
data = log2( data[ , colnames(data) %in% case[ case$expr %in% 1 , ]$patient ] + 1 )


write.table( data , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
