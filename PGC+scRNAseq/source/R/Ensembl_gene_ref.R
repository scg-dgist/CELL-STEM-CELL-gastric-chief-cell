##### Code for Mouse Stomach Pgc+ Single-cell RNA-seq data
##### Author: Somi Kim (ksomi47@dgist.ac.kr)
##### Last Update: 2021/02/17


### Make Mmusculus ensembl gene reference set
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset="mmusculus_gene_ensembl", 
                   host="www.ensembl.org")
ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'chromosome_name', 'gene_biotype'), 
                      mart=ensembl)
rownames(ensemblGenes) <- ensemblGenes[,1]
mtGenes = ensemblGenes[ensemblGenes[,3]=="MT",]
