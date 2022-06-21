library(data.table)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

snv = as.data.frame( fread(  file.path(input_dir, "SNV.txt.gz") , stringsAsFactors=FALSE , sep="\t" ))

effect = sapply( snv[ , "effect" ] , function(x){
  output = NULL
  if( x %in% "" ){
    output = "Silent"
  } else{
    z = unlist( strsplit( x , "[0-9]+" , perl=TRUE ))
    ref = unlist( strsplit( z , "." , fixed=TRUE ))[2]
    alt = z[2]
    output = ifelse( ref==alt, "Silent" , "Missense_Mutation" )
  }
})

snv$ref = ifelse( snv$ref %in% "-" , "" , snv$ref )
snv$alt = ifelse( snv$alt %in% "-" , "" , snv$alt )

data = cbind( snv[ , c( "sample" , "pos" , "ref" , "alt" ) ] ,
              sapply( snv[ , "chr" ] , function(x){ paste( "chr" , x , sep="" ) } ) ,
              apply( snv[ , c( "ref", "alt" ) ] , 1 , function(x){ ifelse( nchar(x[1]) != nchar(x[2]) , "INDEL", "SNV") }  ) ,
              effect
)

colnames(data) = c( "Sample" , "Pos" , "Ref" , "Alt" , "Chr" , "MutType" , "Effect"   )
data = data[ , c( "Sample" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" ) ]

## Get Gene from position
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
coords = getBM(attributes=c("hgnc_symbol"  , "chromosome_name" , "start_position", "end_position"), mart=human)
colnames(coords) = c("gene" , "chr","start","end")		
tab_coords = table( coords$chr )
coords = coords[ coords$chr %in% names( tab_coords[ tab_coords >= 400 ] ) | coords$chr %in% "MT" , ]
coords$chr = paste0( 'chr' , coords$chr )

getGene = function(x, coords){
  gene = coords[ as.character( coords$chr ) %in% as.character( x[2] ) & 
                   as.numeric( as.character( coords$start )) <= as.numeric( as.character( x[3] )) & 
                   as.numeric( as.character( coords$end )) >= as.numeric( as.character( x[3] )) , ]
  g = NULL
  if(nrow(gene)){
    if( nrow(gene) == 1 ){
      g = as.character( gene$gene )
    } 
    if( nrow(gene) > 1 ){
      g = paste( unique( as.character( gene$gene ) ) , collapse="," )
    } 
  } else{ g = NA }
  g
}

chr = sort( unique( data$Chr ) )
snv = NULL
for( i in 1:length( chr ) ){
  dat = data[ data$Chr %in% chr[i] , ]
  coor = coords[ coords$chr %in% chr[i] , ]
  print( chr[i] )
  genes = apply(dat,1, function(x){ getGene( x=x , coords=coor ) } )
  snv = rbind( snv , cbind( dat , genes ) )
}

data = snv
colnames(data) = c( "Sample" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" , "Gene"  )

data$Sample = paste( "P" , data$Sample , sep="" )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
data = data[ as.character( data$Sample ) %in% as.character( case[ case$snv %in% 1 , ]$patient ) , c( "Sample" , "Gene" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" ) ]

output = NULL
for( i in 1:nrow(data) ){
  gene = data$Gene[i]
  if( length( grep( "," , gene ) ) ){
    gene = unlist( strsplit( gene , "," , fixed=TRUE ) )
    dat = cbind( data[ i , ]$Sample , gene , data[ i , c( -1 , -2 ) ] )
    colnames( dat ) = colnames( data )
    output = rbind( output , dat )
  } else{
    output = rbind( output , data[ i , ] )
  }
}

write.table( output , file=file.path(output_dir, "SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
