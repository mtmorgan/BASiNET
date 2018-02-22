#'@title Performs the classification methodology using complex network theory
#'@name classification
#'
#'@description Given two distinct data sets, one of mnRNA and one of lncRNA. 
#'The classification of the data is done from the structure of the networks formed by the sequences. 
#'After this is done classifying with the J48 classifier and randomForest. 
#'It is also created in the current directory a file of type arff called' result 'with all values so that it can be used later. 
#There is also the graphic parameter that when TRUE generates graphs based on the results of each measure.
#'
#'@param mRNA Directory where the file .FASTA lies with the mRNA sequences
#'@param lncRNA Directory where the file is located fasta with lncRNA sequences
#'@param word Integer that defines the size of the word to parse. By default the word parameter is set to 3
#'@param step Integer that determines the distance that will be traversed in the sequences for creating a new connection. By default the step parameter is set to 1
#'@param sncRNA Directory where the file is located fasta with the sncRNA sequences (OPTIONAL)
#'@param graphic Parameter of the logical type, TRUE or FALSE for graphics generation. As default graphic gets FALSE
#'@param graphic3D Parameter of the logical type, TRUE or FALSE for 3D graphics generation. As default graphic3D gets FALSE
#'@param classifier Character Parameter. By default the classifier is J48, but the user can choose to use randomForest by configuring as classifier = "RF"
#'
#'
#' @return Data.frame with the results of measures
#'
#' @author Eric Augusto Ito
#'
#'
#' @examples
#' arqSeqMRNA <- system.file("extdata", "sequences2.fasta", package = "BASiNET")
#' arqSeqLNCRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
#' classification(mRNA=arqSeqMRNA,lncRNA=arqSeqLNCRNA,word=3,step=3,graphic=FALSE,graphic3D=FALSE)
#'
#' @importFrom Biostrings readBStringSet
#' @import igraph
#' @import RWeka
#' @import randomForest
#' @export

classification <- function(mRNA, lncRNA, word=3, step=1, sncRNA, graphic, graphic3D, classifier = c('J48', 'RF')){


	# if(missing(classifier)){
	# 	classifier<-"J48"
	# }

	# if((classifier!="RF")&&(classifier!="J48")){
	# 	print("INVALID CLASSIFIER! Classifier configured for J48")
	# 	classifier<-"J48"
	# }
	classifier <- match.arg(classifier)

	seqMRNA<-readBStringSet(mRNA)
	seqLNCRNA<-readBStringSet(lncRNA)

	if(!missing(sncRNA)){
		seqSNCRNA<-readBStringSet(sncRNA)
		numClass<-3
	}else{
		seqSNCRNA<-NULL
		numClass<-2
	}

	numSeq<-(length(seqMRNA)+length(seqLNCRNA)+length(seqSNCRNA))
	averageShortestPathLengths <- matrix(nrow=numSeq,ncol=1)
	clusteringCoefficient <- matrix(nrow=numSeq,ncol=1)
	standardDeviation <- matrix(nrow=numSeq,ncol=1)
	maximum <- matrix(nrow=numSeq,ncol=1)
	assortativity<- matrix(nrow=numSeq,ncol=1)
	betweenness <- matrix(nrow=numSeq,ncol=1)
	degree <- matrix(nrow=numSeq,ncol=1)
	minimum <- matrix(nrow=numSeq,ncol=1)
	motifs3 <- matrix(nrow=numSeq,ncol=1)
	motifs4 <- matrix(nrow=numSeq,ncol=1)

	for(k in seq_len(numClass)){
		if (k == 1) {
			message("Analyzing mRNA from number: ")
			seq<-seqMRNA
		} else if (k == 2) {
			message("Analyzing lncRNA from number: ")
			seq<-seqLNCRNA
		} else {			# k == 3
			message("Analyzing sncRNA from number: ")
			seq<-seqSNCRNA
		}

		for(x in seq_along(seq)){
			message(x)
			sequence<-strsplit(toString(seq[x]),split='')
			sequence<-sequence[[1]]
			net<-createNet(word, step, sequence)
			limitThreshold<-max(net[])

			n_new_columns <- limitThreshold - ncol(averageShortestPathLengths)
			m_new <- matrix(nrow = numSeq, ncol = n_new_columns)
			averageShortestPathLengths <- cbind(averageShortestPathLengths, m_new)
			clusteringCoefficient <- cbind(clusteringCoefficient, m_new)
			...
			
			vector <- sapply(seq_len(limitThreshold), function(t) {
			    		net<<-threshold(t, net)
			    		measures(net)
			})
			
			cidx <- seq_len(nrow(vector))
			averageShortestPathLengths[aux + x, cidx] <- vector[1,]
			clusteringCoefficient[aux+x, cidx] <- vector[2,]
			...
		}
	}
	message("Rescaling values")
	numCol<-ncol(averageShortestPathLengths)
	
	listMatrix<-list(
	    averageShortestPathLengths,clusteringCoefficient,
	    standardDeviation,maximum,assortativity,betweenness,degree,
	    minimum,motifs3,motifs4)
	listMatrix <- lapply(
	    listMatrix, reschedule,
	    length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA)
	)
	
	namesMeasure<-c(
	    "Average shortest path length", "Cluster Coefficient",
	    "Standard deviation", "Maximum", "Assortativity",
	    "Betweenness", "Degree", "Minimal", "Motifs 3", "Motifs 4"
	)
 

	
	if(missing(graphic)||graphic==FALSE){
	}else{
		if(graphic==TRUE){
			for(i in seq_len(10)){
				createGraph2D(listMatrix[[i]],length(seqMRNA),length(seqLNCRNA), namesMeasure[i])
			}
		}
	}	
	if(missing(graphic3D)||graphic3D==FALSE){
	}else{
		if(graphic3D==TRUE){
			for(x in seq_len(10)){
				for(y in seq_len(10)){
					if(y>x){
						createGraph3D(listMatrix[[x]],listMatrix[[y]],length(seqMRNA),length(seqLNCRNA),namesMeasure[x],namesMeasure[y])
					}
				}
			}
		}
	}

	
	#cria uma janela em branco para exibicao da arvore do J48
	# class1<-matrix(nrow=length(seqMRNA),ncol=1)
	# class2<-matrix(nrow=length(seqLNCRNA),ncol=1)
	# class1[,1]<-"mRNA"
	# class2[,1]<-"lncRNA"
	# classes<-rbind(class1,class2)
	# colnames(classes) <- paste('CLASS') 

	message("Creating data frame")
	data<-cbind(assortativity,betweenness,averageShortestPathLengths, clusteringCoefficient,degree,minimum,maximum,standardDeviation,motifs3,motifs4)
	data<-data.frame(data)
	if(numClass==2){
		data["CLASS"]<-factor(c("lncRNA"), levels = c("mRNA","lncRNA"));
		for(i in seq_along(seqMRNA)){
			data$CLASS[i] <- "mRNA"
		}
	}else{
		if(numClass==3){
			data["CLASS"]<-factor(c("lncRNA"), levels = c("mRNA","lncRNA","sncRNA"));
			data$CLASS[seq_along(seqMRNA)] <- "mRNA"
			data$CLASS[(length(seqMRNA)+length(seqLNCRNA)+1):(length(seqMRNA)+length(seqLNCRNA)+length(seqSNCRNA))] <- "sncRNA"
		}
	}
	data[is.na(data)] <- 0
	rmcfs::write.arff(data, file = "Result.arff")
	message("Result.arff file generated in the current R directory")

	if(classifier=="J48"){
		message("Sorting the data with the J48")
		obj <- J48(CLASS ~ ., data = data)
		result <- evaluate_Weka_classifier(obj, numFolds = 10, complexity = TRUE, seed = 1, class = TRUE)
		print(obj)
		print(result)
		# plot(obj)
	}

	if(classifier=="RF"){
		message("Sorting the data with the Random Forest")
	    set.seed(1)
		rf <- randomForest(data[,seq_along(data[1,])], data[,"CLASS"])
		print(rf)
		print(getTree(randomForest(data[,-40], data[,5], ntree=10), 3, labelVar=TRUE))
	}

	return(invisible(data))
}
