library(data.table)
library(magrittr)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#------------------------------------------ Define some useful functions -----------------------------------------------------------------------
linkage <- function(gff) { # creates a 2-column table with children->parent linkages. Takes original gff annotation as its argument.
  output <- apply(gff[,9], 1, function(x) {
    id <- substr(x, 4, regexpr(';', x) - 1)
    if(regexpr('Parent=', x)[[1]] > 0) { parent <- substr(x, regexpr('Parent=', x) + 7, gregexpr(';', x)[[1]][2] -1); return(c(id, parent))}
    else {return(c(id, "Primary"))}
  }) %>% t() %>% as.data.frame(., stringsAsFactors = F) %>% setNames(., c("ID", "Parent1"))
  return(output)
}

chain <- function(x) { # reconstitutes a full chain of parents for every child record. Takes the output of the linkage function as its argument.
  temp <- x[match(x[, ncol(x)], x$ID), "Parent1"]
  temp[is.na(temp)] <- "Primary"
  temp <- as.data.frame(temp, stringsAsFactors = F) %>% setNames(., paste0("Parent", ncol(x)))
  print(length(unique(temp[,1])))
  if(length(unique(temp[,1])) == 1) {return(x)}
  else{ return(chain(cbind(x, temp))) }
}

remove.features <- function(gff, type, parents){ # removes unwanted feature types from the gff file together with their children and parents
              id <- gff[unlist(gff[,3]) %in% type, 9] %>% apply(., 1, function(x) {substr(x, 4, regexpr(';', x) - 1)})
              #parents <- chain(linkage(gff))
              primary_id <- apply(parents[match(id, parents$ID), ], 1, function(x) {
                record <- x[x != "Primary"] 
                record <- record[length(record)]
                return(record)
              }) %>% unname() %>% unique()
              
              discard_lines <- sapply(parents, function(x) {  x %in% primary_id    }) %>% apply(., 1, function(x) { any(x) }) 
              output <- gff[!discard_lines, ]
              return(output)
}

extract.id <- function(gff, type, level = "Primary", parents){ # simplified version of remove.features() that only reports top-level (Primary) or their own local IDs for given feature types
              id <- gff[unlist(gff[,3]) %in% type, 9] %>% apply(., 1, function(x) {substr(x, 4, regexpr(';', x) - 1)})
              if(level != "Primary") {return(id)}
              else { #parents <- chain(linkage(gff))
                     primary_id <- apply(parents[match(id, parents$ID), ], 1, function(x) {
                        record <- x[x != "Primary"] 
                        record <- record[length(record)]
                        return(record)
                      }) %>% unname() %>% unique()
              return(primary_id)
              }
}

#----------------------------------------------------------------------------------------------------------------------------------------
# Load GFF annotation file
gff    <- fread(file="GRCh38.p12.Refseq.gff", skip = 8, stringsAsFactors = F, header = F, fill = T, na.strings = c("", "NA"), sep="\t") %>% na.omit() #deals with unwanted #comment lines in the gff
gff$ID <- apply(gff[,9], 1, function(x) { id <- substr(x, 4, regexpr(';', x) - 1) })

# create a table of parents-children relations
parents.table <- chain(linkage(gff))

# Create a list of top parents with all childrens listed (Caution: takes ~1 hour) # Save the R object for future use to avoid re-calculations.
parents.tree <- lapply(unique(parents.table[parents.table$Parent1 == 'Primary','ID']), function(x) { return(NULL) }) %>% setNames(unique(parents.table[parents.table$Parent1 == 'Primary', 'ID']))
for(i in 1:nrow(parents.table)) {
          if (!isTRUE(parents.table[i,'ID'] %in% names(parents.tree)) ) {
            children <- parents.table[i, -c(1)] %>% unname() %>% unlist()
            parent   <- children[children %in% names(parents.tree)]
            parents.tree[[parent]] <- c(parents.tree[[parent]], parents.table[i, 'ID']) %>% unique()   }
} 

#saveRDS(parents.tree, file = "parents_tree.rds")
#parents.tree <- readRDS('parents_tree.rds')

# fix gene boundaries to be the same as the span of children features. The discrepancy occured when I removed Gnomon records. Some gene names were shared between BestRefseq and Gnomon as 'BestRefSeq%2CGnomon'. Their boundaries are wider than corresponding BestRefSeq childs.
setindex(gff, ID)
for(name in names(parents.tree)) {
    slice <- gff[name, on = 'ID']

    if(nrow(slice) == 1){
      gene_start <- slice[, V4] 
      gene_end   <- slice[, V5] 

      if(!is.null(parents.tree[[name]]))    {
        children_start <- gff[parents.tree[[name]], V4, on = 'ID'] %>% min()
        children_end   <- gff[parents.tree[[name]], V5, on = 'ID'] %>% max()
        
        if(children_start > gene_start) { gene_start = children_start }
        if(children_end   < gene_end  ) { gene_end   = children_end   }
      } 
      gff[name, on = 'ID', V4 := gene_start]
      gff[name, on = 'ID', V5 := gene_end]
    }
} 
        # con <- file("GRCh38.p12.Refseq.gff", "r")
        # header <- readLines(con, n = 8)
        # write.table(header, file = "GRCh38.p12.Refseq.fixed.gff", col.names = F, row.names = F, quote = F)
        # write.table(gff[,1:9], file = "GRCh38.p12.Refseq.fixed.gff", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
        # close(con); rm(con)

# Remove non-coding features
# check all present top level features: table(gff[,3])
# Tips: NCBI RefSeq. Difference between 'mRNA' and 'transcript' features is the former have CDS and the latter do not.
# Gene --> mRNA --> CDS and exons.
# Gene --> transcript --> exons
# Therefore, do not remove 'transcripts' otherwise you remove their gene parents.

gff2 <- remove.features(gff, c('antisense_RNA','biological_region','cDNA_match','centromere','guide_RNA','lnc_RNA','match','miRNA','primary_transcript','pseudogene','region','RNase_MRP_RNA','RNase_P_RNA','rRNA','scRNA','snoRNA','snRNA','telomerase_RNA','tRNA','vault_RNA','Y_RNA'), parents.table)
        con <- file("GRCh38.p12.Refseq.gff", "r")
        header <- readLines(con, n = 8)
        write.table(header, file = "GRCh38.p12.Refseq.coding.gff", col.names = F, row.names = F, quote = F)
        write.table(gff2[,1:9], file = "GRCh38.p12.Refseq.coding.gff", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
        close(con); rm(con)




















#---------------------------------------------------------------------------------------------------------------------------------------
# #gff2 <- copy(gff)
# gff <- copy(gff2)
# setindex(gff, ID)
# 
# coding_gene.id     <- intersect(extract.id(gff, "CDS", level = "Primary", parents.table), extract.id(gff, "mRNA", level = "Primary", parents.table) ) # foolproof method because some coding genes might have only one type (CDS or mRNA) listed in annotation by mistake
# coding_gene.lines  <- match(coding_gene.id, gff$ID)  #line indexing for instant search
# 
# primary.id    <- names(parents.tree)                 # speeds up execution of next snippet
# primary.lines <- match(primary.id, gff$ID)           # speeds up execution of next snippet
# 
# 
# 
# 
# for(name in names(parents.tree)) {
#   print(name)
#   gene_model <- gff[c(name, parents.tree[[name]]), on = 'ID']
# #  rna.id <- parents.tree[[name]][grep('rna', parents.tree[[name]])]
#   
# 
# #  rna <- gene_model[gene_model$ID == rna.id & gene_model$V3 == 'mRNA', c('V4', 'V5')]
#   
#   
# 
#   #find lower neighbour
#   lower_neighbor_model <- data.table()
#   flag <- 0
#   limit <- 1
#   while(flag == 0) {
#     neighbor <- names(parents.tree)[match(name, names(parents.tree)) + limit]
#     if(!is.na(neighbor)) { 
#        lower_neighbor_model <-  gff[c(neighbor, parents.tree[[neighbor]]), on = 'ID'] 
#        if(lower_neighbor_model[1, V7] == gene_model[1, V7] & lower_neighbor_model[1, V1] == gene_model[1, V1]) { flag <- 1  }
#        else {  limit = limit + 1; if(limit == 100) {flag <- 2}  }
#     }
#     else{flag <- 2}
#   }
#   if(flag == 2) {lower_neighbor_model <- 'abscent'} 
#   
#   #find upper neighbour
#   upper_neighbor_model <- data.table()
#   flag <- 0
#   limit <- 1
#   while(flag == 0) {
#     neighbor <- names(parents.tree)[match(name, names(parents.tree)) - limit]
#     if(length(neighbor) > 0) { 
#       upper_neighbor_model <-  gff[c(neighbor, parents.tree[[neighbor]]), on = 'ID'] 
#       if(upper_neighbor_model[1, V7] == gene_model[1, V7] & upper_neighbor_model[1, V1] == gene_model[1, V1]) { flag <- 1  }
#       else {  limit = limit + 1; if(limit == 100) {flag <- 2}  }
#     }
#     else{flag <- 2}
#   }
#   if(flag == 2) {upper_neighbor_model <- 'abscent'} 
#     
# 
#   
#       
#         #boundaries of exons, mRNAs and CDSs
#         if(length(gene_model[V3 == 'exon', V4]) != 0) {exon.left  <- gene_model[V3 == 'exon', c('V4', 'ID')]}  else {exon.left  <- 'abscent'}
#         if(length(gene_model[V3 == 'exon', V5]) != 0) {exon.right <- gene_model[V3 == 'exon', c('V5', 'ID')]}  else {exon.right <- 'abscent'}
#         if(length(gene_model[V3 == 'mRNA', V4]) != 0) {mRNA.left  <- gene_model[V3 == 'mRNA', c('V4', 'ID')]}  else {mRNA.left  <- 'abscent'}
#         if(length(gene_model[V3 == 'mRNA', V5]) != 0) {mRNA.right <- gene_model[V3 == 'mRNA', c('V5', 'ID')]}  else {mRNA.right <- 'abscent'}
#         if(length(gene_model[V3 == 'CDS',  V4]) != 0) {CDS.left   <- gene_model[V3 == 'CDS',  c('V4', 'ID')]}  else {CDS.left   <- 'abscent'}
#         if(length(gene_model[V3 == 'CDS',  V5]) != 0) {CDS.right  <- gene_model[V3 == 'CDS',  c('V5', 'ID')]}  else {CDS.right  <- 'abscent'}  
#       
#         if(typeof(CDS.left)  != 'character' & typeof(mRNA.left)  != 'character') {UTR5 <- min(CDS.left[,1])   - min(mRNA.left[,1])}  else {UTR5 <- 'abscent'}
#         if(typeof(CDS.right) != 'character' & typeof(mRNA.right) != 'character') {UTR3 <- max(mRNA.right[,1]) - max(CDS.right[,1])}  else {UTR3 <- 'abscent'}
#         
#         if(UTR5 != 'abscent' & typeof(upper_neighbor_model) != 'character') {
#            if(UTR5 < 100 & gene_model[1, V4] - upper_neighbor_model[1,V5] > 1000) {
#                gff[name, on = 'ID', V4 := V4 - 100 + UTR5]
#                gff[mRNA.left[mRNA.left$V4 == min(mRNA.left[,1]), ID], on = 'ID', V4 := V4 - 100 + UTR5]
#                gff[exon.left[exon.left$V4 == min(exon.left[,1]), ID], on = 'ID', V4 := V4 - 100 + UTR5]     }
#         }
#         
#         if(UTR3 != 'abscent' & typeof(lower_neighbor_model) != 'character')  {
#            if(UTR3 < 100 & lower_neighbor_model[1,V4] - gene_model[1, V5] > 1000) {
#               gff[name, on = 'ID', V5 := V5 + 100 - UTR3]
#               gff[mRNA.right[mRNA.right$V5 ==  min(mRNA.right[,1]), ID], on = 'ID', V5 := V5 + 100 - UTR3]
#               gff[exon.right[exon.right$V5 ==  min(exon.right[,1]), ID], on = 'ID', V5 := V5 + 100 - UTR3]  }
#         }
#   
# 
#     
# }
#   
#   
#   
#   
# 















