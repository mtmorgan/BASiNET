#'@title Performs the classification methodology using complex network theory
#'@name classification
#'
#'@description Given two distinct data sets, one of mnRNA and one of lncRNA.
#'The classification of the data is done from the structure of the networks formed by the sequences.
#'After this is done classifying with the J48 classifier and randomForest. It is also created in the current directory
#'a file of type arff called' result 'with all values ​​so that it can be used later. There is
#'also the graphic parameter that when TRUE generates graphs based on the results of each measure.
#'
#'@param word Size of the word to parse
#'@param step # to be given to each call
#'@param mRNA Directory where the file lies with the mRNA sequences
#'@param lncRNA Directory where the file is located fasta with lncRNA sequences
#'@param sncRNA Directory where the file is located fasta with the sncRNA sequences (OPTIONAL)
#'@param graphic TRUE or FALSE for graphics generation. As default graphic gets FALSE
#'@param graphic3D TRUE or FALSE for 3D graphics generation. As default graphic3D gets FALSE
#'
#' @details
#'
#' @return Void
#'
#' @author Eric Augusto Ito
#'
#' @seealso
#'
#' @examples
#' arqSeqMRNA <- system.file("extdata", "sequences2.fasta", package = "classificador")
#' arqSeqLNCRNA <- system.file("extdata", "sequences.fasta", package = "classificador")
#' classification(3, 3, arqSeqMRNA,arqSeqLNCRNA)
#'
#' @import seqinr 
#' @import igraph
#' @import RWeka
#' @import randomForest
#' @export

classification <- function(word, step, mRNA, lncRNA, sncRNA, graphic, graphic3D){

	seqMRNA<-read.fasta(file=mRNA,forceDNAtolower=FALSE)
	seqLNCRNA<-read.fasta(file=lncRNA,forceDNAtolower=FALSE)

	if(!missing(sncRNA)){
		seqSNCRNA<-read.fasta(file=sncRNA,forceDNAtolower=FALSE)
		numClass<-3
	}else{
		seqSNCRNA<-NULL
		numClass<-2
	}

	numSeq<-(length(seqMRNA)+length(seqLNCRNA)+length(seqSNCRNA))
	caminhoMinimoMedio <- matrix(nrow=numSeq,ncol=1)    
	coeficienteDeCluster <- matrix(nrow=numSeq,ncol=1) 
	desvioPadrao <- matrix(nrow=numSeq,ncol=1) 
	maximo <- matrix(nrow=numSeq,ncol=1) 
	assortativity<- matrix(nrow=numSeq,ncol=1) 
	betweenness <- matrix(nrow=numSeq,ncol=1) 
	degree <- matrix(nrow=numSeq,ncol=1) 
	minimo <- matrix(nrow=numSeq,ncol=1) 
	motifs3 <- matrix(nrow=numSeq,ncol=1) 
	motifs4 <- matrix(nrow=numSeq,ncol=1)

	for(k in 1:numClass){
		
		if(k==1){
			print("Analyzing mRNA from number: ")
			aux<-0
			seq<-seqMRNA
			type<-"mRNA"
		}else{
			if(k==2){
				print("Analyzing lncRNA from number: ")
				seq<-seqLNCRNA
				aux<-length(seqMRNA)
				type<-"lncRNA"
			}else{
				if(k==3){
					print("Analyzing sncRNA from number: ")
					seq<-seqSNCRNA
					aux<-(length(seqMRNA)+length(seqLNCRNA))
				}
			}
		}

		for(x in 1:length(seq)){
			print(x)
			limitThreshold<-0
			sequence<-getSequence(seq[[x]])
			net<-createNet(word, step, sequence)
			limitThreshold<-max(net[])

			for(t in 1:limitThreshold){
				if(t>length(caminhoMinimoMedio[1,])){
					caminhoMinimoMedio <- cbind(caminhoMinimoMedio,matrix(nrow=numSeq,ncol=1))
					coeficienteDeCluster <- cbind(coeficienteDeCluster,matrix(nrow=numSeq,ncol=1))
					desvioPadrao <- cbind(desvioPadrao,matrix(nrow=numSeq,ncol=1))
					maximo <- cbind(maximo,matrix(nrow=numSeq,ncol=1))
					assortativity <- cbind(assortativity,matrix(nrow=numSeq,ncol=1))
					betweenness <- cbind(betweenness,matrix(nrow=numSeq,ncol=1))
					degree <-cbind(degree,matrix(nrow=numSeq,ncol=1))
					minimo <- cbind(minimo,matrix(nrow=numSeq,ncol=1))
					motifs3 <- cbind(motifs3,matrix(nrow=numSeq,ncol=1))
					motifs4 <- cbind(motifs4,matrix(nrow=numSeq,ncol=1))
				}

				net<-threshold(t, net)
				vector<-measures(net)
				caminhoMinimoMedio[aux+x,t] <- vector[1]    
				coeficienteDeCluster[aux+x,t] <- vector[2] 
				degree[aux+x,t] <- vector[3] 
				assortativity[aux+x,t] <- vector[4] 
				betweenness[aux+x,t] <- vector[5] 
				desvioPadrao[aux+x,t] <- vector[6] 
				maximo[aux+x,t] <- vector[7] 
				minimo[aux+x,t] <- vector[8] 
				motifs3[aux+x,t] <- vector[9] 
				motifs4[aux+x,t] <- vector[10] 
				
			}
		}
	}
	print("Rescaling values")
	numCol<-length(caminhoMinimoMedio[1,])
	caminhoMinimoMedio<-reescalar(caminhoMinimoMedio,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(caminhoMinimoMedio) <- paste('CMM',1:numCol)
	coeficienteDeCluster<-reescalar(coeficienteDeCluster,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(coeficienteDeCluster) <- paste('CdC', 1:numCol)
	desvioPadrao<-reescalar(desvioPadrao,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(desvioPadrao) <- paste('DP', 1:numCol)
	maximo<-reescalar(maximo,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(maximo) <- paste('MAX', 1:numCol)
	assortativity<-reescalar(assortativity,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(assortativity) <- paste('ASS', 1:numCol)
	betweenness<-reescalar(betweenness,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(betweenness) <- paste('BET', 1:numCol)
	degree<-reescalar(degree,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(degree) <- paste('DEG', 1:numCol)
	minimo<-reescalar(minimo,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(minimo) <- paste('MIN', 1:numCol)
	motifs3<-reescalar(motifs3,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(motifs3) <- paste('MT3', 1:numCol)
	motifs4<-reescalar(motifs4,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(motifs4) <- paste('MT4', 1:numCol)

	listMatrix<-list(caminhoMinimoMedio,coeficienteDeCluster,desvioPadrao,maximo,assortativity,betweenness,degree,minimo,motifs3,motifs4)
	namesMeasure<-c("Average path length", "Cluster Coefficient", "Standard deviation", "Maximum", "Assortativity", "Betweenness", "Degree", "Minimal", "Motifs 3", "Motifs 4")
 

	
	if(missing(graphic)||graphic==FALSE){
	}else{
		if(graphic==TRUE){
			for(i in 1:10){
				createGraph2D(listMatrix[[i]],length(seqMRNA), namesMeasure[i])
			}
		}
	}	
	if(missing(graphic3D)||graphic==FALSE){
	}else{
		if(graphic3D==TRUE){
			for(x in 1:10){
				for(y in 1:10){
					if(y>x){
						createGraph3D(listMatrix[[x]],listMatrix[[y]],length(seqMRNA),namesMeasure[x],namesMeasure[y])
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

	print("Creating data frame")
	data<-cbind(assortativity,betweenness,caminhoMinimoMedio, coeficienteDeCluster,degree,minimo,maximo,desvioPadrao,motifs3,motifs4)
	data<-data.frame(data)
	if(numClass==2){
		data["CLASS"]<-factor(c("lncRNA"), levels = c("mRNA","lncRNA"));
		for(i in 1:length(seqMRNA)){
			data$CLASS[i] <- "mRNA"
		}
	}else{
		if(numClass==3){
			data["CLASS"]<-factor(c("lncRNA"), levels = c("mRNA","lncRNA","sncRNA"));
			data$CLASS[1:length(seqMRNA)] <- "mRNA"
			data$CLASS[(length(seqMRNA)+length(seqLNCRNA)+1):(length(seqMRNA)+length(seqLNCRNA)+length(seqSNCRNA))] <- "sncRNA"
		}
	}
	data[is.na(data)] <- 0
	rmcfs::write.arff(data, file = "Result.arff")
	print("Result.arff file generated in the current R directory")
	print("Sorting the data with the J48")
	obj <- J48(CLASS ~ ., data = data)
	result <- evaluate_Weka_classifier(obj, numFolds = 10, complexity = TRUE, seed = 1, class = TRUE)
	print(obj)
	print(result)
	# plot(obj)
	print("Sorting the data with the Random Forest")
    set.seed(1)
	rf <- randomForest(data[,1:length(data[1,])], data[,"CLASS"])
	print(rf)
	print(getTree(randomForest(data[,-40], data[,5], ntree=10), 3, labelVar=TRUE))

	return()
}