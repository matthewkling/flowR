

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
            separate(value, c("value", "ci"), sep="\\(") %>%
            mutate_each(funs(strip)) %>%
            select(-edge) %>%
            mutate_each(funs(as.numeric))

      return(list(populations=p, data=d))
}



#' Run BayesAss.
#'
#' @param iter Number of iterations (integer)
#' @param burn Burn-in iterations (integer)
#' @param interval Sampling interval (integer)
#' @param seed Random starting seed (integer)
#' @param other Any other args to BA3.exe
#' @param exe Path to BayesAss executable file
#' @param infile Path to input gene data file
#' @param wd Path to working directory for file export
#' @param tag Text identifier to append to outfiles
#'
#' @export
bayesass <- function(iter=10000000,
                     burn=1000000,
                     interval=1000,
                     seed=10,
                     exe="E:/flow/BA3Windows64/BA3.exe",
                     infile="infile_bayesass.txt",
                     other=NULL,
                     wd=getwd(),
                     tag=NULL){
      setwd(wd)
      system(paste0(exe,
                    " -u -v -g -t -m0.10 -a0.45 -f0.45",
                    " -i", iter,
                    " -b", burn,
                    " -s", seed,
                    " -n", interval,
                    " ", other,
                    " ",
                    infile),
             wait=FALSE, invisible=FALSE)

      if(!is.null(tag)){
            files <- c("BA3out.txt", "BA3trace.txt", "BA3indiv.txt")
            file.rename(files, gsub(".txt", paste0("_", tag, ".txt"), files))
      }
}

