test_that("get_genotype_class categorises correctly", {
	expect_equal(get_genotype_class("0/0"), "HomoRef")
	expect_equal(get_genotype_class("1/1"), "HomoAlt")
	expect_equal(get_genotype_class("0/1"), "HeteroRef")
	expect_equal(get_genotype_class("2/1"), "HeteroAlt")
	expect_equal(get_genotype_class(as.character(NA)), as.character(NA))
	expect_error(get_genotype_class(c("0/0", "2/1")))
})

