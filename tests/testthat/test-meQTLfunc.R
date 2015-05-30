test_that("setModel function works", {
	expect_error(setModel(''))
	expect_equal(setModel('linear'),modelLINEAR)
	expect_equal(setModel('anova'),modelANOVA)
	expect_equal(setModel('linear_cross'),modelLINEAR_CROSS)
})

test_that("setSNPOptions function works", {
	snps<-setSNPOptions("\t", "NA")
	expect_equal(snps$fileDelimiter, "\t")
	expect_equal(snps$fileOmitCharacters, "NA")
	snps<-setSNPOptions(" ", "NaN")
    expect_equal(snps$fileDelimiter, " ")
    expect_equal(snps$fileOmitCharacters, "NaN")
})
