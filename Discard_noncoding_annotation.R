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
gff    <- fread(file="GRCh38.p12.Refseq.gff", skip = 8, stringsAsFactors = F, header = F, fill = T, na.strings = c("", "NA"), sep="\t")
gff    <- na.omit(gff) # deals with unwanted #comment lines in the gff
gff$ID <- apply(gff[,9], 1, function(x) { id <- substr(x, 4, regexpr(';', x) - 1) })

# create a table of parents-children relations
parents.table <- chain(linkage(gff))

# Create a list of top parents with all childrens listed (Caution: takes ~1 hour) # Save the R object for future use to avoid re-calculations.
parents.tree <- lapply(parents.table[parents.table$Parent1 == 'Primary','ID'], function(x) { return(NULL) }) %>% setNames(parents.table[parents.table$Parent1 == 'Primary',] %>% .[,'ID'])
for(i in 1:nrow(parents.table)) {
          if (!isTRUE(parents.table[i,'ID'] %in% names(parents.tree)) ) {
            children <- parents.table[i, -c(1)] %>% unname() %>% unlist()
            parent   <- children[children %in% names(parents.tree)]
            parents.tree[[parent]] <- c(parents.tree[[parent]], parents.table[i, 'ID']) %>% unique()   }
} #saveRDS(parents.tree, file = "parents_tree.rds")
#parents.tree <- readRDS('parents_tree.rds')

# fix gene boundaries to be the same as the span of children features. The discrepancy occured when I removed Gnomon records. Some gene names were shared between BestRefseq and Gnomon as 'BestRefSeq%2CGnomon'. Their boundaries are wider than corresponding BestRefSeq childs.
setindex(gff, ID)
for(name in names(parents.tree)) {
    slice <- gff[name, on = 'ID']

    if(nrow(slice) == 1){
      gene_start <- slice[, V4] 
      gene_end   <- slice[, V5] 

      if(!is.null(parents.tree[[name]]))    {
        children_start <- sapply(parents.tree[[ name ]], function(x) {  return(gff[x, V4, on = 'ID'] %>% min()     )     }) %>% min()
        children_end   <- sapply(parents.tree[[ name ]], function(x) {  return(gff[x, V5, on = 'ID'] %>% max()     )     }) %>% max()

        if(children_start > gene_start) { gene_start = children_start }
        if(children_end   < gene_end  ) { gene_end   = children_end   }
      } 
      gff[name, on = 'ID', V4 := gene_start]
      gff[name, on = 'ID', V5 := gene_end]
    }
} 
        #con <- file("GRCh38.p12.Refseq.gff", "r")
        #header <- readLines(con, n = 8)
        #write.table(header, file = "GRCh38.p12.Refseq.fixed.gff", col.names = F, row.names = F, quote = F)
        #write.table(gff2[,1:9], file = "GRCh38.p12.Refseq.fixed.gff", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
        #close(con); rm(con)

# Remove non-coding features
# check all present top level features: table(gff[,3])
# Tips: NCBI RefSeq. Difference between 'mRNA' and 'transcript' features is the former have CDS and the latter do not.
# Gene --> mRNA --> CDS and exons.
# Gene --> transcript --> exons
# Therefore, do not remove 'transcripts' otherwise you remove their gene parents.

gff2 <- remove.features(gff, c('antisense_RNA','biological_region','cDNA_match','centromere','guide_RNA','lnc_RNA','match','miRNA','primary_transcript','pseudogene','region','RNase_MRP_RNA','RNase_P_RNA','rRNA','scRNA','snoRNA','snRNA','telomerase_RNA','tRNA','vault_RNA','Y_RNA'), parents.table)
        # con <- file("GRCh38.p12.Refseq.coding.gff", "r")
        # header <- readLines(con, n = 8)
        # write.table(header, file = "GRCh38.p12.Refseq.coding.gff", col.names = F, row.names = F, quote = F)
        # write.table(gff2[,1:9], file = "GRCh38.p12.Refseq.coding.gff", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
        # close(con); rm(con)











# # Extract coding genes ID that have both CDS and mRNA children
# gene.id   <- intersect(extract.id(gff, "CDS", level = "Primary", parents.table), extract.id(gff, "mRNA", level = "Primary", parents.table) ) # foolproof method because some coding genes might have only one type (CDS or mRNA) listed in annotation
# gene.lines    <- match(gene.id, gff$ID)
# 
#         # check if these genes are proper genes and not some weird IDs
#         # if (all(grepl('gene', gene.id))) {print('All found gene IDs are in the form of "geneX" ')} else { print('Some gene IDs are not the form of "geneX" ') }
#  
# primary.id    <- unique(parents.table[parents.table[,2] == "Primary", 1]) # speeds up execution of next snippet
# primary.lines <- match(primary.id, gff$ID)                                # speeds up execution of next snippet
# 
# gene.record <- lapply(1:2, function(x) {
#   gene <- gff[gene.lines[x], ] %>% unlist() %>% unname()
#   
#   #find upper neighbour on the same strand
#   upper <- function(a){
#     temp <- gff[primary.lines[match(gene.lines[x], primary.lines) - a], ] %>% unlist() %>% unname()
#     if(temp[1] == gene[1]) {
#       if (temp[7] == gene[7]) { return(temp) }   else { return (upper(a + 1)) }
#     }
#     else {return(NULL)}
#   }
#   
#   #find lower neighbour on the same strand
#   lower <- function(a){
#     temp <- gff[primary.lines[match(gene.lines[x], primary.lines) + a], ] %>% unlist() %>% unname()
#     if(temp[1] == gene[1]) {
#       if (temp[7] == gene[7]) { return(temp) }   else { return (lower(a + 1)) }
#     }
#     else {return(NULL)}
#   }
#   
#   upper.neighbour <- upper(1)
#   lower.neighbour <- lower(1)
#   
#   #extract gff lines for the gene by taking a range between gff line with gene name and the lower neighbour
#   if (length(lower.neighbour) > 0) { lines <- gff[gene.lines[x]:(match(lower.neighbour[10], gff$ID)-1), ]  }
#   else { print(gene.id[x])
#     if(!is.na(unique(unname(unlist(gff[,1])))[match(gene[1], unique(unname(unlist(gff[,1])))) + 1])) { lines <- gff[gene.lines[x]:match((unique(unname(unlist(gff[,1]))))[match(gene[1], unique(unname(unlist(gff[,1])))) + 1], unname(unlist(gff[,1]))), ] }
#     else { lines <- gff[gene.lines[x]:nrow(gff), ] }
#   }
#   #discard lines for features on the opposite strand
#   lines <- as.data.frame(lines)
#   lines <- lines[lines$V7 == gene[7], ]
#   
#   #boundaries of exons, mRNAs and CDSs
#   exon.left.ends <- c()
#   exon.right.ends <- c()
#   mRNA.left.ends <- c()
#   mRNA.right.ends <- c()
#   CDS.left.ends <- c()
#   CDS.right.ends <- c()
#   for (i in 1:nrow(lines)) { 
#     if(lines[i,3] == 'exon') { 
#       exon.left.ends = c(exon.left.ends, lines[i,4]) 
#       exon.right.ends = c(exon.right.ends, lines[i,5])   } 
#     else if (lines[i,3] == 'mRNA') {
#       mRNA.left.ends = c(mRNA.left.ends, lines[i,4])
#       mRNA.right.ends = c(mRNA.right.ends, lines[i,5])   }
#     else if (lines[i,3] == 'CDS') {
#       CDS.left.ends <- c(CDS.left.ends, lines[i,4])
#       CDS.right.ends <- c(CDS.right.ends, lines[i,5])
#     }
#   }
#   exon.limits <- c(min(exon.left.ends), max(exon.right.ends))
#   mRNA.limits <- c(min(mRNA.left.ends), max(mRNA.right.ends))
#   CDS.limits <- c(min(CDS.left.ends), max(CDS.right.ends))
#   
#   utr.left <- exon.limits[1] - CDS.limits[1]
#   utr.right <- exon.limits[2] - CDS.limits[2]
#   
#   gap.left <- as.integer(gene[4]) - as.integer(upper.neighbour[5]) - 1
#   gap.right <- as.integer(lower.neighbour[4]) - as.integer(gene[5]) - 1
#   
#   new.gene.left = as.integer(gene[4])
#   new.gene.right = as.integer(gene[5])
#   if(utr.left  < 100 & gap.left  > 1000)  { 
#     new.gene.left  <- as.integer(gene[4]) - 100 + utr.left
#     new.mRNA.left  <- new.gene.left
#     new.exon.left  <- new.gene.left
#   }
#   if(utr.right < 100 & gap.right > 1000)  { 
#     new.gene.right <- as.integer(gene[5]) + 100 - utr.right 
#     new.mRNA.right <- new.gene.right
#     new.exon.right <- new.gene.right
#   }
#   
#   return(list( "gene"  = list('id' = gene.id[x], 
#                               'coordinates' = gene[4:5],
#                               'lines' = lines,
#                               'exon_limits' = exon.limits, 
#                               'mRNA_limits' = mRNA.limits, 
#                               'CDS_limits' = CDS.limits,
#                               'utr_left' = utr.left,
#                               'utr_right' = utr.right,
#                               'gap_left' = gap.left,
#                               'gap_right' = gap.right,
#                               'new_gene_left' = new.gene.left,
#                               'new_gene_right' = new.gene.right), 
#                "upper" = list("id" = upper.neighbour[10], "coordinates" = upper.neighbour[4:5]),
#                "lower" = list("id" = lower.neighbour[10], "coordinates" = lower.neighbour[4:5]) ))
#   
#   
# })  %>% setNames(., gene.id[1:2])    
# 








