library(preprocessCore)

quantNorm <- function(filename)
{
	rawData<-read.table(filename,header=T,row.names=1)
	mat <- as.matrix(rawData)
	norm <- normalize.quantiles(mat,copy=TRUE)
	normData<-data.frame(rownames(rawData),norm)
	colnames(normData)<-c("ID",colnames(rawData))
	write.table(normData,paste(filename,'.norm',sep=""),sep='\t',quote=F, row.names=F)
}
