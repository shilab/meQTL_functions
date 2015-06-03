test_that("setModel function works", {
	expect_error(setModel(''))
	expect_equal(setModel('linear'),modelLINEAR)
	expect_equal(setModel('anova'),modelANOVA)
	expect_equal(setModel('linear_cross'),modelLINEAR_CROSS)
})

test_that("setSNPOptions function works", {
	fileOptions<-setFileOptions("\t","NA")
	snps<-fileOptions$snps;
	gene<-fileOptions$gene;
	cvrt<-fileOptions$cvrt;
	expect_equal(snps$fileDelimiter, "\t")
	expect_equal(snps$fileOmitCharacters, "NA")
    expect_equal(gene$fileDelimiter, "\t")
    expect_equal(gene$fileOmitCharacters, "NA")
    expect_equal(cvrt$fileDelimiter, "\t")
    expect_equal(cvrt$fileOmitCharacters, "NA")

	fileOptions<-setFileOptions(" ","NaN")
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
