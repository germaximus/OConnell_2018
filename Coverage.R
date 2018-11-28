library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#############################################################################################################
###########  Prepare mRNA coverage files for analysis  #######################################################
#############################################################################################################
      #Instuction how to load Riso-seq coverage data for every sample, separated by individual reads. Takes a lot of time, therefore I saved the result as Rdata file.
      # library("R.utils")
      # coverage_samples <- vector("list", length=length(list.files(path="./sample_coverage", pattern='.coverage')))
      # names(coverage_samples) <- list.files(path="./sample_coverage", pattern='.coverage')
      # for(i in 1:length(coverage_samples)) {
      #   handle <- file(paste0("./sample_coverage/", names(coverage_samples)[i]), open = "rb"); linecount <- countLines(handle);  close(handle);
      #   handle <- file(paste0("./sample_coverage/", names(coverage_samples)[i]), open = "r");
      #   file <- vector("list", length=linecount)
      #   for(a in 1:linecount) {
      #     line <- readLines(handle, n=1)
      #     line <- unlist(strsplit(line, "\t"))
      #     names(file)[a] <- line[1]
      #     file[[a]] <- as.integer(line[-1])
      #   }
      #   close(handle)
      #   # in early version of the Coverage.pl genes were written in a random order (as unsorted hash keys). To fix this, sort file names:
      #   file <- file[sort(names(file))]
      #   coverage_samples[[i]] <- file
      # }
      #  dir.create("Rdata")
      #  saveRDS(coverage_samples, file = "./Rdata/coverage_samples.Rdata")

readRDS(coverage_samples, file="./Rdata/coverage_samples.Rdata")
matching.table <- read.table("./sample_coverage/sample_table.txt", header = TRUE, stringsAsFactors = FALSE)

#accepts coverage_samples Rdata as the input. length - mRNA length cut-off (includes 5`-UTR), save - should the file be saved
coverage_processor <- function(input, length = 2000, align_by = 'start', save = FALSE, report = FALSE) {
  if(align_by == 'start'){ #align transcripts by start codon
    coverage_filtered <- lapply(input, function(x) { lapply(x, function(y)  {   if(length(y) >= length) {return(y[1:length])}})}) #leaves only mRNA longer than specified by length
  } 
  else if(align_by == 'stop') { #align transcripts by stop codon
    coverage_filtered <- lapply(input, function(x) { lapply(x, function(y)  {   if(length(y) >= length) {return(tail(y, length))}})}) #leaves only mRNA longer than specified by length
  }
  
  coverage_filtered <- lapply(coverage_filtered, function(x) {x[!sapply(x, is.null)] }) 
  coverage_normalized <- lapply(coverage_filtered, function(x) {
    lapply(x, function(y) {
      average <- mean(y)  
      if(average >= 0.25){  return(y/average)   }
      else {return(NULL)}
    })
  }); rm(coverage_filtered)
  coverage_normalized <- lapply(coverage_normalized, function(x) {x[!sapply(x, is.null)] })
  coverage_matrix     <- lapply(coverage_normalized, function(x) { matrix(unlist(x), byrow=TRUE, nrow=length(x))  })
  densities_sample    <- lapply(coverage_matrix, function(x) {colSums(x) / nrow(x)})
  if(isTRUE(save)) { saveRDS(densities_sample, file = "./Rdata/density_samples.Rdata")  }
  if(isTRUE(report)) { print(paste0("processing is completed for cut-off ",length))  }
  return(densities_sample)
} #example: coverage_processor(input = coverage_samples, length = 2000, save = FALSE, report = TRUE)

density_start <- coverage_processor(input = coverage_samples, length = 1500, align_by = 'start', save = FALSE, report = FALSE)
density_stop  <- coverage_processor(input = coverage_samples, length = 1500, align_by = 'stop', save = FALSE, report = FALSE)

#plot stop site
x11(height=3, width=6) 
par(cex=0.6, cex.lab=2, cex.axis=2, cex.main=2, lend = 1, lwd = 3)
plot(NULL, xlab="", ylab="", yaxt='n', xaxt='n', xlim= c(1,1500), ylim= c(0,3))
abline(v = 100, lty = 2, col = 'grey25')
lines(density_start[["mRNA_0128-01_p6_1.coverage"]], type="l", col="grey50")
lines(density_start[["ribo_0128-01_p6_1.coverage"]], type="l", col="#00a087")
lines(density_start[["ribo_2127_p9_2.coverage"]], type="l", col="#e64b35")
lines(density_start[["ribo_HDF_p5_2.coverage"]], type="l", col="#3c5488")
axis(1, at=c(0, 500, 1000, 1500), labels=c(-100, 400, 900, 1400), lwd = 3)
axis(2, las=1, at=c(0,1,2,3), labels=c(0,1,2,3), lwd = 3)
dev.copy2pdf(file="Start_codon_coverage.pdf")
dev.off()

#plot stop site
x11(height=3, width=6)
par(cex=0.6, cex.lab=2, cex.axis=2, cex.main=2, lend = 1, lwd = 3)
plot(NULL, xlab="", ylab="", yaxt='n', xaxt='n', xlim= c(1,1500), ylim= c(0,3))
abline(v = 1400, lty = 2, col = 'grey25')
lines(density_stop[["mRNA_0128-01_p6_1.coverage"]], type="l", col="grey50")
lines(density_stop[["ribo_0128-01_p6_1.coverage"]], type="l", col="#00a087")
lines(density_stop[["ribo_2127_p9_2.coverage"]], type="l", col="#e64b35")
lines(density_stop[["ribo_HDF_p5_2.coverage"]], type="l", col="#3c5488")
axis(1, at=c(0, 500, 1000, 1500), labels=c(-1400, -900, -400, 100), lwd = 3)
axis(2, las=1, at=c(0,1,2,3), labels=c(0,1,2,3), lwd = 3)
dev.copy2pdf(file="Stop_codon_coverage.pdf")
dev.off()







