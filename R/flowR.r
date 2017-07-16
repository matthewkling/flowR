

library(diveRsity)
library(dplyr)
library(tidyr)
select <- dplyr::select


#' Convert MIGRATE-N infiles to divMigrate format.
#'
#' This function takes a MIGRATE-N infile (genepop format) and converts it into
#' a divMigrate infile format.
#'
#' @param infile Path to genepop text file
#' @param outfile Path to output file
#' @export
mn2dm <- function(infile, outfile=NULL){

      require(stringr)

      lines <- readLines(infile)
      n_loci <- as.integer(strsplit(lines[1], " ")[[1]][2])
      loci_names <- paste0("locus", 1:n_loci)
      loci <- paste(loci_names, collapse=", ")
      lines <- c(lines[1], loci, lines[2:length(lines)]) # insert locus labels on second line
      lines <- sub("_ ", "_ , ", lines) # comma after each individual's name
      lines <- gsub("\\?", "000", lines) # missing data is encoded as zeroes in popgen
      lines[nchar(lines)<25] <- "pop" # these have to be called pop
      lines[1] <- "title" # title can't start with a number


      # add leading zeroes to all alleles so that they're 3 digits long
      pad <- function(x) str_pad(x, 3, "left", 0)
      expand <- function(x){
            if(all(nchar(x)==max(nchar(x)))) return(x)
            xd <- data.frame(a=x) %>%
                  separate(a, c("a1", "a2")) %>%
                  mutate_each(funs(pad), a1:a2) %>%
                  unite(a, a1:a2, sep=".")
            return(xd$a)
      }
      ar <- setdiff(which(nchar(lines)>10), 2) # rows with allele data
      a <- data.frame(raw=lines[ar]) %>%
            separate(raw, c("ind", "alleles"), sep=",", remove=F) %>%
            mutate(alleles=str_trim(alleles)) %>%
            separate(alleles, loci_names, sep=" ") %>%
            mutate_each_(funs(expand), loci_names) %>%
            select(-raw) %>%
            mutate(ind=paste0(ind, ", ")) %>%
            apply(1, paste, collapse=" ")

      lines[ar] <- a

      lines <- gsub("\\.", "", lines) # no separator between alleles

      if(is.null(outfile)) outfile <- paste0(dirname(infile), "/infile_divmigrate.txt")
      con <- file(outfile)
      writeLines(lines, con)
      close(con)
      return(outfile)
}


#' Convert MIGRATE-N infiles to BayesAss format.
#'
#' This function takes a MIGRATE-N infile (genepop format) and converts it into
#' a BayesAss infile format.
#'
#' @param infile Path to genepop text file
#' @param outfile Path to output file
#' @export
mn2ba <- function(infile, outfile=NULL){

      lines <- readLines(infile)
      lines <- gsub("\\?", "0", lines)

      # Bring the genotype data into a data frame
      # I assume that the individual lines are longer than the first line
      # That may not be a good assumption
      dframe <- data.frame(matrix(unlist(strsplit(lines[nchar(lines)>nchar(lines[1])], " ")), length(lines[nchar(lines)>nchar(lines[1])]), byrow=T))

      # Fill in the population names as the second column of the data frame
      # I assume that the population lines are shorter than the first line
      # That may not be a good assumption either
      dframe[[2]] <- rep(sapply(strsplit(lines[nchar(lines)<nchar(lines[1])], " "),function(x) x[2]), as.integer(sapply(strsplit(lines[nchar(lines)<nchar(lines[1])], " "),function(x) x[1])))

      # Unpack individuals into separate lines for each locus
      # Loci are arbitrarily named "locus1" to "locusn"
      # Also split the locus genotypes into separate columns for each allele
      dframe2 <- data.frame()
      for(i in 3:ncol(dframe)) dframe2 <- rbind(dframe2, data.frame(c(dframe[c(1:2)], L=paste0("locus", i), data.frame(do.call('rbind', strsplit(as.character(dframe[[i]]), '.', fixed=TRUE))))))

      # Save a text file with given name or "infile_bayesass.txt" (default)
      if(is.null(outfile)) outfile <- paste0(dirname(infile), "/infile_bayesass.txt")
      write.table(dframe2, outfile, sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)

      return(outfile)
}




#' Convert divMigrate infiles to BayesAss format.
#'
#' This function takes a divMigrate infile and converts it into
#' a BayesAss infile format. Currently does not support haplotype.
#'
#' @param infile Path to genepop text file
#' @param outfile Path to output file
#' @export
dm2ba <- function(infile, outfile=NULL){

      require(stringr)
      require(dplyr)
      require(tidyr)

      lines <- readLines(infile)
      header <- lines[1:2]
      lines <- lines[3:length(lines)]
      loci <- unlist(str_split(header[2], ", "))

      # remove individual IDs
      lines <- substr(lines, regexpr(",", lines)+2, nchar(lines))

      # split into populations
      i1 <- which(grepl("pop", lines))
      i2 <- c((i1 - 1)[2:length(i1)], length(lines))
      p <- lapply(1:length(i1), function(i){ lines[(i1[i]+1):i2[i]] })

      # convert to long format
      p <- lapply(1:length(p), function(i){
            x <- p[[i]]
            x <- data.frame(ind = paste0("ind",
                                         str_pad(i, nchar(length(p)), "right", 0),
                                         str_pad(1:length(x), 3, "right", 0)),
                            pop = paste0("pop", i),
                            geno = str_trim(x)) %>%
                  separate(geno, into=loci, sep=" ", remove=T) %>%
                  gather(locus, geno, -ind, -pop)
            x <- separate(x, geno, into=c("allele1", "allele2"), sep=nchar(x$geno[1])/2)
      })

      p <- do.call("rbind", p)

      # export
      if(is.null(outfile)) outfile <- paste0(dirname(infile), "/infile_BayesAss.txt")
      write.table(p, outfile, sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
      return(outfile)
}





#' Load MIGRATE-N outfile.
#'
#' This function loads MIGRATE-N results as a data frame.
#'
#' @param path Path to genepop text file
#' @return a data frame of migration parameters
#' @export
load_mn_results <- function(path){
      require(dplyr)
      require(tidyr)
      require(data.table)
      d <- readLines(path)
      d <- d[which(d=="Bayesian estimates")+3:length(d)]
      d <- d[1:(which(d=="")[1]-2)]
      d <- d[!grepl("---------------", d)]
      writeLines(d, "temp.txt")
      d <- fread("temp.txt")
      file.remove("temp.txt")

      gf <- d %>%
            filter(Locus=="All", grepl("M_", Parameter)) %>%
            dplyr::select(Parameter, median) %>% # using mode or mean rather than median yields very similar results
            separate(Parameter, c("from", "to"), sep="->") %>%
            mutate(from=as.integer(sub("M_", "", from)),
                   to=as.integer(to))
      gen_flow <- spread(gf, to, median) %>%
            dplyr::select(-from) %>%
            as.matrix()
      gen_flow[is.na(gen_flow)] <- 0

      return(gen_flow)
}



#' Load BayesAss outfile.
#'
#' This function loads BayesAss results as a data frame. Currently it retrieves
#' only the pairwise population stats, but it may later be extended to retrieve
#' individual ancestry estimates.
#'
#' @param path Path to genepop text file
#' @return a list of two data frames, one with population id numbers and one
#'   with migration parameters
#' @export
load_ba_results <- function(path){
      require(dplyr)
      require(tidyr)
      require(stringr)

      # load file
      d <- readLines(path)

      # grab and munge population dictionary
      p <- d[which(d==" Population Index -> Population Label:")+2]
      p <- str_trim(p)
      p <- data.frame(raw=str_split(p, " ")[[1]]) %>%
            separate(raw, c("id", "label"), sep="->") %>%
            mutate(id=as.integer(id))

      # grab pairwise migration data
      d <- d[(which(d==" Migration Rates:")+2) :
                   (which(d==" Inbreeding Coefficients:")-2)]

      # get data into clean tabular format
      strip <- function(x) gsub("\\[|\\]|\\(|\\)", "", x)
      d <- data.frame(raw=d) %>%
            separate(raw, paste0("e", 0:length(d)), sep=" m") %>%
            select(-e0) %>%
            gather(edge, value) %>%
            separate(value, c("pair", "value"), sep=": ") %>%
            separate(pair, c("p1", "p2"), sep="\\]\\[") %>%
            separate(value, c("value", "sd"), sep="\\(") %>%
            mutate_each(funs(strip)) %>%
            select(-edge) %>%
            mutate_each(funs(as.numeric))

      return(list(populations=p, data=d))
}



#' Run BayesAss.
#'
#' @param infile Path to input gene data file
#' @param outfile Output filename including .txt
#' @param iter Number of iterations (integer)
#' @param burn Burn-in iterations (integer)
#' @param interval Sampling interval (integer)
#' @param seed Random starting seed (integer)
#' @param other Any other args to BA3.exe
#' @param exe Path to BayesAss executable file
#'
#' @export
BayesAss <- function(infile="infile_bayesass.txt",
                     outfile="outfile_bayesass.txt",
                     iter=10000000,
                     burn=1000000,
                     interval=1000,
                     deltaM=0.10,
                     deltaA=0.10,
                     deltaF=0.10,
                     seed=10,
                     other=NULL,
                     exe="E:/flow/BA3Windows64/BA3.exe"){

      y <- system(paste0(exe,
                    " -u -v", #" -g -t",
                    #" -m0.10 -a0.45 -f0.45",
                    " -i", iter,
                    " -b", burn,
                    " -s", seed,
                    " -n", interval,
                    " -a", deltaA,
                    " -f", deltaF,
                    " -m", deltaM,
                    " ", other,
                    " -o", outfile,
                    " ", infile),
             wait=TRUE, invisible=FALSE, intern=TRUE)
      return(y)
}



#' Get acceptance rates from BayesAss output.
#'
#' @param x An object returned by BayesAss
#'
#' @export
acceptance <- function(x){
      require(dplyr)
      require(stringr)
      x <- x[grepl("accepted", x)] %>%
            str_split("accepted") %>%
            unlist()
      x <- x[substr(x, 1, 1)==":"]
      x <- substr(x, 4, nchar(x))
      i <- gregexpr(")", x) %>% sapply(function(x) x[1])
      x <- substr(x, 1, i-1) %>%
            str_split(", ") %>%
            lapply(as.numeric) %>%
            do.call("cbind", .) %>%
            t()
      colnames(x) <- c("migrate", "indiv", "allele", "inbreed", "missing")
      return(x)
}




#' Calculate Bayesian deviance from a BayesAss trace file.
#'
#' Code modified from supplementary files for: Meirmans, P. G. Nonconvergence in
#' Bayesian estimation of migration rates. 726–733 (2014).
#' doi:10.1111/1755-0998.12216
#'
#' @param file Path to BA3trace.txt
#' @param burn Number of burn iterations used
#' @param interval Sampling intrval used
#'
#' @return Deviance value
#'
#' @export
BA_deviance <- function(file, burn, interval){

      # Read the data from the trace file
      trace=read.table(file, header=TRUE)

      # Calculate the deviance
      range = (trace$State > burn & trace$State %% interval == 0)
      D = -2*mean(trace$LogProb[range])

}



