test_that("getIndex works", {
	me<-list(cis=list(eqtls=list(FDR=c(0.1,0.05,0.01,0.2))))
	expect_equal(getIndex(TRUE,me,0.1),c(1,2,3))
	#TODO: Will need to fix
	expect_equal(getIndex(TRUE,me,0),integer(0))
	me<-list(trans=list(eqtls=list(FDR=c(0.1,0.05,0.01,0.2))))
	expect_equal(getIndex(FALSE,me,0.1),c(1,2,3))
})

test_that("getEQTLS works", {
	me<-list(cis=list(eqtls<-data.frame(snps=c('snp1','snp2'),gene=c('gene2','gene2'),statistic=c(6.09,5.97),pvalue=c(0.04e-01,4.23e-1),FDR=c(0.02,0.012),beta=c(3.78,4.49)) ))

	index <- c(1)

	res <- me$cis$eqtls[1,]

	expect_equal(getEQTLS(TRUE,me,index),res)

    me<-list()
    me$trans<-list() 
    me$trans$eqtls<-data.frame(snps=c('snp1','snp2'),gene=c('gene2','gene2'),statistic=c(6.09,5.97),pvalue=c(0.04e-01,4.23e-1),FDR=c(0.02,0.012),beta=c(3.78,4.49)) 

	index <-c(2)

	res <- me$trans$eqtls[2,]

	expect_equal(getEQTLS(FALSE,me,index),res)
})

#TODO: Add non-trivial correlation example
#TODO: Add test with NAs
test_that("getCorr works", {
	genotype <- c(1,1,1,1,1,1,2,2,2,2,2,2)
	phenotype <- c(1.5,1.5,1.5,1.5,1.5,1.5,2.5,2.5,2.5,2.5,2.5,2.5)	

	expect_equal(getCorr(phenotype, genotype), c(1))

	phenotype <- c(0,1,0,1,0,1,1,0,1,0,1,0)

	expect_equal(getCorr(phenotype, genotype), c(0))
})

test_that("getFileData works", {
	readData <- getFileData('data/snps')
	expectedData <- data.frame(p1=c(1,0,0),p2=c(2,2,2),p3=c(0,2,1),p4=c(1,0,2),p5=c(2,2,2))
	rownames(expectedData) <- c('snp1','snp2','snp3')
	expect_equal(readData,expectedData)

	readData <- getFileData('data/expr')
	expectedData <- data.frame(p1=c(2.86511,7.171009,4.036252),p2=c(2.574775,8.350821,1.955288),p3=c(5.352287,7.714878,2.051268),p4=c(3.849516,7.249339,2.697464),p5=c(3.019739,5.440955,1.892165))
	rownames(expectedData) <- c('gene1','gene2','gene3')
	expect_equal(readData,expectedData)
})
