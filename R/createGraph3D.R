#'@title Creates a three-dimensional chart between two measurements and the threshold
#'@name createGraph3D
#'
#'@description For an analysis of each measure, the createGraph2D () function was created in order to visualize the behavior of each measurement in relation to the threshold. This function creates a graph (Measure x Measure x Threshold) from two arrays, mRNA sequences are given the blue color, the lncRNA sequences are given a red color. In cases where there is a third class this will be given the green color
#'
#'@param matrix1 matrix of the first measure for the creation of the three-dimensional chart
#'@param matrix2 matrix of the second measure for the creation of the three-dimensional chart
#'@param numSeqMRNA number of mRNA sequences
#'@param numSeqLNCRNA number of lncRNA sequences
#'@param nameMeasure1 name of the first measure to put in the title of the graph
#'@param nameMeasure2 name of the second measure to put in the title of the graph
#'
#'@return Show a graphic three-dimensional
#'@author Eric Augusto Ito
#'
#'@import rgl



createGraph3D <- function(matrix1, matrix2, numSeqMRNA, numSeqLNCRNA, nameMeasure1, nameMeasure2){

	somador<-(1/(length(matrix1[1,])-1))
	threshold<-seq(0,1,somador)
	open3d()
	plot3d(matrix1[1,],matrix2[1,], threshold, xlab=nameMeasure1, ylab=nameMeasure2, zlab="Threshold", col=c("blue"), type="l2")
	for(i in 2:length(matrix1[,1])){
		if(i<=numSeqMRNA){
			lines3d(matrix1[i,], matrix2[i,], threshold,color=c("blue"))	
		}else{
			if(i<=(numSeqLNCRNA+numSeqMRNA)){
				lines3d(matrix1[i,], matrix2[i,], threshold,color=c("red"))	
			}else{
				lines3d(matrix1[i,], matrix2[i,], threshold,color=c("green"))	
			}
		}
	}

	return()
}