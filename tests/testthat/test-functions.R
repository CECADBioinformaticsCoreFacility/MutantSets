test_that("get_genotype_class categorises correctly", {
	expect_equal(get_genotype_class("0/0"), "HomoRef")
	expect_equal(get_genotype_class("1/1"), "HomoAlt")
	expect_equal(get_genotype_class("0/1"), "HeteroRef")
	expect_equal(get_genotype_class("2/1"), "HeteroAlt")
	expect_equal(get_genotype_class(as.character(NA)), as.character(NA))
	expect_error(get_genotype_class(c("0/0", "2/1")))
})

test_that("EFF strings parsed", {
	ex <- tibble(EFF = "downstream_gene_variant(MODIFIER||283|c.*341_*342insG|210|C10A4.7|protein_coding|CODING|C10A4.7||TACCCCCCCCCCCCCCAAACACAAT)")
	exParsed <- tibble(
		Effect = "downstream_gene_variant", Effect_Impact = "MODIFIER", 
		Functional_Class = "", Codon_Change = "283",
		Amino_Acid_Change = "c.*341_*342insG", 
		Amino_Acid_Length = "210", Gene_Name = "C10A4.7",
		Transcript_BioType = "protein_coding", 
		Gene_Coding = "CODING", Transcript_ID = "C10A4.7", Exon_Rank = "", 
		Genotype_Number = "TACCCCCCCCCCCCCCAAACACAAT",
		ERRORS = "", WARNINGS = ""
	)
	expect_equal(eff_parser(ex), exParsed)
	exraw <- "downstream_gene_variant(MODIFIER||283|c.*341_*342insG|210|C10A4.7|protein_coding|CODING|C10A4.7||TACCCCCCCCCCCCCCAAACACAAT),intergenic_region(MODIFIER|||n.7413327_7413328insC|||||||TACCCCCCCCCCCCCCAAACACAAT)"
	expect_equal(eff_extractor(exraw)[1,], exParsed)

	expect_equal(
		split_geno(tibble(a = "a,b", EFF = exraw, ALT = "A,A")),
		tibble(a = c("a", "b"), EFF = strsplit(exraw,",")[[1]], ALT = c("A", "A"))
	)
	
	expect_equal(
		split_geno(tibble(a = "a", EFF = "a,b,c,d", ALT = "A")),
		tibble(a = "a", EFF = "a,b,c,d", ALT = "A")
	)
	expect_equal(
		split_geno(tibble(a = "a,b,c", EFF = "a,b,c,d", ALT = "A,A,A")),
		tibble(a = c("a", "b", "c"), EFF = c("a,d","b","c"), ALT = "A")
	)
	expect_equal(
		split_geno(tibble(a = "a,b", EFF = "a,b,c,d", ALT = "A,A")),
		tibble(a = c("a", "b"), EFF = c("a,c","b,d"), ALT = "A")
	)
})

test_that("Embedded wormbase", {
	expect_equal(
		wormbase_view("III", 100000, 105000, 1002500, 1002501),
		"<embed width = '100%' height = '500px' src='https://wormbase.org/tools/genome/jbrowse-simple/?data=data%2Fc_elegans_PRJNA13758&loc=III%3A1e+05..105000&tracks=Curated_Genes%2CClassical_alleles&highlight=III%3A1002500..1002501'>"
	)
	expect_equal(
		wormbase_view("III", 100000, 105000),
		"<embed width = '100%' height = '500px' src='https://wormbase.org/tools/genome/jbrowse-simple/?data=data%2Fc_elegans_PRJNA13758&loc=III%3A1e+05..105000&tracks=Curated_Genes%2CClassical_alleles'>"
	)
	expect_error(
		wormbase_view("chrIII", 100000, 105000)
	)
})
