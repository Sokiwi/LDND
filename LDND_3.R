LDND <- function(A, B) {
	# A and B are two ASJP language names
	# read the data file, remove % and trim whitespaces, 
	# but only if it is not already in memory
	if ( length(grep("in.memory",ls(envir=.GlobalEnv)))==0 ) {
		in.memory <- "nothing"
	} else {
		in.memory <- get("in.memory", envir=.GlobalEnv)
	}
	if ( in.memory!="asjp" ) {
		cat("The data is not in memory, it will take a few minutes to process it.\n")
		cat("But next time you do a distance it will only take a split second.\n\n")
		cat("reading data...\n")
		data <- read.table(file="listss17_formatted.tab", header=TRUE, sep="\t", quote="", na.strings="", strip.white=TRUE, comment.char="", stringsAsFactors=FALSE, colClasses="character",numerals="no.loss")
		data <- data[,c(1,11,12,13,21,22,28,29,31,32,33,35,38,40,41,44,49,50,51,53,54,57,58,61,63,64,67,68,71,76,82,84,85,87,92,95,96,102,105,106,110)]
		cat("getting rid of the loanword symbol...\n")
		# get rid of the loanword symbol
		eliminate.loans <- function(x) {
			where.loan <- grep("%", x)
			if(length(where.loan) > 0) {
				for (i in 1:length(where.loan)) {
					x[where.loan[i]] <- strsplit(x[where.loan[i]], "%")[[1]][2]
				}
			}
			return(x)
		}
		data[,2:41] <- apply(data[,2:41], 2, eliminate.loans)
		# trim whitespace
		cat("trimming whitespace...\n")
		trim <- function(x) {
			whites <- which(strsplit(x, "")[[1]]==" ")
			if( length(whites) > 0 ) {
				x <- paste(strsplit(x, " ")[[1]][1:(length(whites)+1)],collapse="")
			} else {
				x <- x
			}
		return(x)
		}
		data <- apply(data, c(1,2), trim)

		# get rid of more than two synonyms
		cat("getting rid of more than two synonyms...\n")
		synred <- function(x) {
			if( length(strsplit(x,",")[[1]]) > 2 ) {
				x <- paste(strsplit(x,",")[[1]][1:2],collapse=",")
			} else {
				x <- x
			}
		return(x)
		}
		data <- apply(data, c(1,2), synred)

		# get a vector of unique symbols, including strings
		# to be treated as one symbol
		# collecting all unique sound representations
		cat("getting a list of all unique sound representations...\n")
		# find all unique sound representations (comma included)
		# first a function for splitting words, sw
		sw <- function(x) {
			ws <- strsplit(x,"")[[1]]
			# QUOTES
			wq <- grep("\"",ws)
			if ( length(wq) > 0 ) {
				for (i in 1:length(wq)) {
					ws[wq[i]-1] <- paste(ws[wq[i]-1],ws[wq[i]],sep="")
				}
			ws <- ws[-wq]
			}
			# STARS
			wst <- grep("\\*",ws)
			if ( length(wst) > 0 ) {
				for (i in 1:length(wst)) {
					ws[wst[i]-1] <- paste(ws[wst[i]-1],ws[wst[i]],sep="")
				}
			ws <- ws[-wst]
			}
			# TILDE
			wt <- grep("\\~",ws)
			if ( length(wt) > 0 ) {
				for (i in 1:length(wt)) {
					ws[wt[i]-2] <- paste(ws[wt[i]-2],ws[wt[i]-1],sep="")
				}
			ws <- ws[-c(wt,wt-1)]
			}
			# DOLLAR
			wd <- grep("\\$",ws)
			if ( length(wd) > 0 ) {
				for (i in 1:length(wd)) {
					ws[wd[i]-3] <- paste(ws[wd[i]-3],ws[wd[i]-2],ws[wd[i]-1],sep="")
				}
			ws <- ws[-c(wd,wd-1,wd-2)]
			}
		return(unique(ws))
		}
		l <- lapply(data[,2:41],sw)
		uni <- unique(unlist(l))
		rm(l)
		# the following is a trick to get the comma in position 44
		wcomma <- which(uni==",")
		uni <- uni[-wcomma]
		uni <- c(uni[1:43],",",uni[44:length(uni)])
		assign("uni", uni, envir=.GlobalEnv)

		# translate sound-encoding sequences to numbers
		# then to utf8 and then to words
		cat("transforming sound representations to single utf8 symbols...\n")
		transf <- function(x) {
			ws <- strsplit(x,"")[[1]]
			# QUOTES
			wq <- grep("\"",ws)
			if ( length(wq) > 0 ) {
				for (i in 1:length(wq)) {
					ws[wq[i]-1] <- paste(ws[wq[i]-1],ws[wq[i]],sep="")
				}
			ws <- ws[-wq]
			}
			# STARS
			wst <- grep("\\*",ws)
			if ( length(wst) > 0 ) {
				for (i in 1:length(wst)) {
					ws[wst[i]-1] <- paste(ws[wst[i]-1],ws[wst[i]],sep="")
				}
			ws <- ws[-wst]
			}
			# TILDE
			wt <- grep("\\~",ws)
			if ( length(wt) > 0 ) {
				for (i in 1:length(wt)) {
					ws[wt[i]-2] <- paste(ws[wt[i]-2],ws[wt[i]-1],sep="")
				}
			ws <- ws[-c(wt,wt-1)]
			}
			# DOLLAR
			wd <- grep("\\$",ws)
			if ( length(wd) > 0 ) {
				for (i in 1:length(wd)) {
					ws[wd[i]-3] <- paste(ws[wd[i]-3],ws[wd[i]-2],ws[wd[i]-1],sep="")
				}
			ws <- ws[-c(wd,wd-1,wd-2)]
			}
			for (i in 1:length(ws)) {
				ws[i] <- which(uni==ws[i])
			}
			out <- intToUtf8(ws)
		}
		data[,2:41] <- apply(data[,2:41], c(1,2), transf)
		assign("in.memory", "asjp", envir=.GlobalEnv)
		assign("data", data, envir=.GlobalEnv)
	}
	# reduce the matrix to the two relevant languages
	whereA <- which(data[,1]==A)
	whereB <- which(data[,1]==B)
	data.red <- data[c(whereA, whereB),2:41]
	# do LDN for individual word pairs in data.red
	LDN <- function(x,y) {
		LDN <- as.vector(adist(x,y)/max(nchar(x),nchar(y)))
		return(LDN)
	}
	# reduce further when something is missing for the purpose
	# of LDN between words referring to the same concept
	missing <- union(which(data.red[1,]=="\004\004\004"),which(data.red[2,]=="\004\004\004"))
	if(length(missing) > 0) {
		data.red1 <- data.red[,-missing]
	} else {
		data.red1 <- data.red
	}
	# do average LDN for all pairs of words referring to the same concept,
	# take into account synonyms
	LDNs <- rep(0,length(data.red1[1,]))
	for (i in 1:length(data.red1[1,])) {
		L1 <- length(grep(",",data.red1[1,i]))
		L2 <- length(grep(",",data.red1[2,i]))
		if ( L1==0 & L2==0 ) {
			LDNs[i] <- LDN(data.red1[1,i],data.red1[2,i])
		} else if ( L1 > 0 & L2==0 ) {
			w1 <- strsplit(data.red1[1,i],",")[[1]][1]
			w2 <- strsplit(data.red1[1,i],",")[[1]][2]
			w3 <- data.red1[2,i]
			LDNs[i] <- mean(c(LDN(w1,w3),LDN(w2,w3)))
		} else if ( L1==0 & L2 > 0 ) {
			w1 <- data.red1[1,i]
			w2 <- strsplit(data.red1[2,i],",")[[1]][1]
			w3 <- strsplit(data.red1[2,i],",")[[1]][2]
			LDNs[i] <- mean(c(LDN(w1,w2),LDN(w1,w3)))
		} else {
			w1 <- strsplit(data.red1[1,i],",")[[1]][1]
			w2 <- strsplit(data.red1[1,i],",")[[1]][2]
			w3 <- strsplit(data.red1[2,i],",")[[1]][1]
			w4 <- strsplit(data.red1[2,i],",")[[1]][2]
			LDNs[i] <- mean(c(LDN(w1,w3),LDN(w1,w4),LDN(w2,w3),LDN(w2,w4)))
		}
	}
	LDN.lg <- mean(LDNs)
	list1 <- data.red[1,][data.red[1,]!="\004\004\004"]
	list2 <- data.red[2,][data.red[2,]!="\004\004\004"]
	names.grid <- expand.grid(names(list1),names(list2))
	diagonal <- which(as.vector(names.grid[,1])==as.vector(names.grid[,2]))
	w <- expand.grid(list1,list2,stringsAsFactors=FALSE)[-diagonal,] # w is grid of word pairs
	L <- length(w[,1])
	gamma <- rep(0,L)
	for (i in 1:L) {
		c1 <- length(grep(",",w[i,1]))
		c2 <- length(grep(",",w[i,2]))
		if ( c1==0 & c2==0 ) {
			gamma[i] <- LDN(w[i,1],w[i,2])
		} else if ( c1 > 0 & c2==0 ) {
			w1 <- strsplit(w[i,1],",")[[1]][1]
			w2 <- strsplit(w[i,1],",")[[1]][2]
			w3 <- w[i,2]
			gamma[i] <- mean(c(LDN(w1,w3),LDN(w2,w3)))
		} else if ( c1==0 & c2 > 0 ) {
			w1 <- w[i,1]
			w2 <- strsplit(w[i,2],",")[[1]][1]
			w3 <- strsplit(w[i,2],",")[[1]][2]
			gamma[i] <- mean(c(LDN(w1,w2),LDN(w1,w3)))
		} else {
			w1 <- strsplit(w[i,1],",")[[1]][1]
			w2 <- strsplit(w[i,1],",")[[1]][2]
			w3 <- strsplit(w[i,2],",")[[1]][1]
			w4 <- strsplit(w[i,2],",")[[1]][2]
			gamma[i] <- mean(c(LDN(w1,w3),LDN(w1,w4),LDN(w2,w3),LDN(w2,w4)))
		}
	}
	LDND_out <- round(LDN.lg/mean(gamma),4)
	return(LDND_out)
}

