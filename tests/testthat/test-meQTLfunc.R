test_that("setModel function works", {
	expect_error(setModel(''))
	expect_equal(setModel('linear'),modelLINEAR)
	expect_equal(setModel('anova'),modelANOVA)
	expect_equal(setModel('linear_cross'),modelLINEAR_CROSS)
})

test_that("getGenotypes works", {
	snps<-getGenotypes("\t","NA",TRUE, TRUE, 'data/snps')
	expect_equal(snps$fileDelimiter, "\t")
	expect_equal(snps$fileOmitCharacters, "NA")
	expect_equal(snps$GetAllRowNames(), c('snp1','snp2','snp3'))
	expect_equal(snps$columnNames, c('p1', 'p2', 'p3', 'p4', 'p5'))
	mat <- as.matrix(snps$getSlice(1))
	expected_mat = matrix(c(1,2,0,1,2,0,2,2,0,2,0,2,1,2,2),nrow=3,byrow=T)
	expect_equal(mat, expected_mat)
})

test_that("setSNPOptions function works", {
	fileOptions<-setFileOptions("\t","NA",TRUE,TRUE)
	gene<-fileOptions$gene;
	cvrt<-fileOptions$cvrt;
    expect_equal(gene$fileDelimiter, "\t")
    expect_equal(gene$fileOmitCharacters, "NA")
    expect_equal(cvrt$fileDelimiter, "\t")
    expect_equal(cvrt$fileOmitCharacters, "NA")

	fileOptions<-setFileOptions(" ","NaN",TRUE,TRUE)
    gene<-fileOptions$gene;
    cvrt<-fileOptions$cvrt;
    expect_equal(gene$fileDelimiter, " ")
    expect_equal(gene$fileOmitCharacters, "NaN")
    expect_equal(cvrt$fileDelimiter, " ")
    expect_equal(cvrt$fileOmitCharacters, "NaN")
})

test_that("skip rows with header", {
	header<-TRUE
	expect_equal(rowsToSkip(header),1)
	header<-FALSE
	expect_equal(rowsToSkip(header),0)
})

test_that("skip columns with rownames", {
	rownames<-TRUE
	expect_equal(colsToSkip(rownames),1)
	rownames<-FALSE
	expect_equal(colsToSkip(rownames),0)
})

test_that("covariates filename", {
	expect_equal(getCovariates(''),character())
	expect_equal(getCovariates('cov_file'),'cov_file')
})
