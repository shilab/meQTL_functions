test_that("getIndex works", {
	me$cis$eqtls<-list()
	me$cis$eqtls$FDR<-c(0.1,0.05,0.01,0.2)
	expect_equal(getIndex(TRUE,me,0.1),c(1,2,3))
	#TODO: Will need to fix
	expect_equal(getIndex(TRUE,me,0),integer(0))
	me$trans$eqtls<-list()
	me$trans$eqtls$FDR<-c(0.1,0.05,0.01,0.2)
	expect_equal(getIndex(FALSE,me,0.1),c(1,2,3))
})

test_that("getEQTLS works", {
	me<-list()
	me$cis<-list()
	me$cis$eqtls<-data.frame(snps=c('snp1','snp2'),gene=c('gene2','gene2'),statistic=c(6.09,5.97),pvalue=c(0.04e-01,4.23e-1),FDR=c(0.02,0.012),beta=c(3.78,4.49))

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
