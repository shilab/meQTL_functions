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

	snps<-getGenotypes(" ", "NaN", TRUE, TRUE, 'data/snps2')
	expect_equal(snps$fileDelimiter, " ")
	expect_equal(snps$fileOmitCharacters, "NaN")
})

test_that("getExpression works", {
	gene<-getExpression("\t","NA",TRUE,TRUE, 'data/expr')
	expect_equal(gene$fileDelimiter, "\t")
	expect_equal(gene$fileOmitCharacters, "NA")
	expect_equal(gene$GetAllRowNames(), c('gene1', 'gene2', 'gene3'))
	expect_equal(gene$columnNames, c('p1', 'p2', 'p3', 'p4', 'p5'))
	mat <- as.matrix(gene$getSlice(1))
	expected_mat = matrix(c(2.86511,2.574775,5.352287,3.849516,3.019739,7.171009,8.350821,7.714878,7.249339,5.440955,4.036252,1.955288,2.051268,2.697464,1.892165),nrow=3, byrow=T)
	expect_equal(mat, expected_mat)

	gene<-getExpression(" ","NA",TRUE,TRUE, 'data/expr2')
	expect_equal(gene$fileDelimiter, " ")
	mat <- as.matrix(gene$getSlice(1)) 
	expect_equal(mat, expected_mat)
})

test_that("getCovariates function works", {
	cvrt<-getCovariates("\t","NA",TRUE,TRUE, 'data/cvrt')
    expect_equal(cvrt$fileDelimiter, "\t")
    expect_equal(cvrt$fileOmitCharacters, "NA")
	expect_equal(cvrt$GetAllRowNames(), c('cov1', 'cov2'))
	expect_equal(cvrt$columnNames, c('p1', 'p2', 'p3', 'p4', 'p5'))
	mat <- as.matrix(cvrt$getSlice(1))
	expected_mat = matrix(c(0,0,1,0,1,1,0,1,1,0),nrow=2,byrow=T)
	expect_equal(mat, expected_mat)

	cvrt<-getCovariates(" ","NaN",TRUE,TRUE, 'data/cvrt2')
    expect_equal(cvrt$fileDelimiter, " ")
	mat <- as.matrix(cvrt$getSlice(1))
	expect_equal(mat, expected_mat)
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
	expect_equal(getCovFileName(''),character())
	expect_equal(getCovFileName('cov_file'),'cov_file')
})
