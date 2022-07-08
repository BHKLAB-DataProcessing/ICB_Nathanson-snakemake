library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")

clin_original = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )
selected_cols <- c( "sample","Age","Gender","Primary.Melanoma.type","M.stage","OS","Alive","Response.duration..weeks.")
clin = cbind( clin_original[ , selected_cols ] , "Melanoma" , "CTLA4", NA, NA , NA, NA , NA , NA )
colnames(clin) = c( "patient" , "age" , "sex" , "histo" , "stage" , "t.os"  ,"os","response.other.info", "primary" , "drug_type" , "recist" , "t.pfs" , "pfs" , "dna" , "rna" , "response" )

clin$os = ifelse( clin$os %in% 1 , 0 , 1 )
clin$t.os = clin$t.os * 12
clin$response.other.info = ifelse( clin$response.other.info <= 24 , "NR" , 
							ifelse( clin$response.other.info > 24 , "R" , NA ))


clin$stage = ifelse( clin$stage %in% 3 , "III" , 
				ifelse( clin$stage %in% 4 , "IV" , NA )) 

clin$response = Get_Response( data=clin )
clin$patient = paste( "P" , clin$patient , sep="" )
clin_original$sample = paste( "P" , clin_original$sample , sep="" )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
clin$rna[ clin$patient %in% case[ case$expr %in% 1 , ]$patient ] = "fpkm"
clin$dna[ clin$patient %in% case[ case$snv %in% 1 , ]$patient ] = "wes"

rownames(clin) = clin$patient
clin = clin[ case$patient , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]
rownames(clin_original) <- clin_original$sample
clin_original <- clin_original[case$patient, ]

clin <- format_clin_data(clin_original, 'sample', selected_cols, clin)
clin$Ipi.dosing..mg.kg.x... <- str_replace_all(clin$Ipi.dosing..mg.kg.x..., ';', ',')
clin$BRAF.NRASstatus <- str_replace_all(clin$BRAF.NRASstatus, ';', ',')

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

