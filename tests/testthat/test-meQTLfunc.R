test_that("setModel function works", {
	expect_error(setModel(''))
	expect_equal(setModel('linear'),modelLINEAR)
	expect_equal(setModel('anova'),modelANOVA)
	expect_equal(setModel('linear_cross'),modelLINEAR_CROSS)
})

test_that("setSNPOptions function works", {
	fileOptions<-setFileOptions("\t","NA",TRUE,TRUE)
	snps<-fileOptions$snps;
	gene<-fileOptions$gene;
	cvrt<-fileOptions$cvrt;
	expect_equal(snps$fileDelimiter, "\t")
	expect_equal(snps$fileOmitCharacters, "NA")
    expect_equal(gene$fileDelimiter, "\t")
    expect_equal(gene$fileOmitCharacters, "NA")
    expect_equal(cvrt$fileDelimiter, "\t")
    expect_equal(cvrt$fileOmitCharacters, "NA")

	fileOptions<-setFileOptions(" ","NaN",TRUE,TRUE)
	snps<-fileOptions$snps;
    gene<-fileOptions$gene;
    cvrt<-fileOptions$cvrt;
    expect_equal(snps$fileDelimiter, " ")
    expect_equal(snps$fileOmitCharacters, "NaN")
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
