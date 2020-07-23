
# Global data -----------------------------------------------------------------
# 
gtfnames <- c(
	"sequence", "source", "feature", "start", "end",
	"score", "strand","phase","attributes"
)
# seq <- ape::as.DNAbin(getSeq(BSgenome.Celegans.UCSC.ce11))
# vcfin <- vcfR::read.vcfR(
# "~/Downloads/Galaxy171-[SnpSift_Filter_on_data_158].vcf", verbose = FALSE
# )

# ann <- as.data.frame( # Package this data!!
# 	vroom::vroom(
# 		"/media/richardjamesacton/ext/data/genomes/ce11.ensGene.gtf.gz",
# 		col_names = gtfnames
# 	)
# )

# var type colours
var_type_colours <- c(
	"complex" = "#e41a1c",
	"del" = "#377eb8",
	"snp" = "#4daf4a",
	"mnp" = "#984ea3",
	"ins" = "#ff7f00"
)

# Functions -------------------------------------------------------------------
## |Effect Annotation ---------------------------------------------------------

#' effParser
#' parsing snpeff strings EFF field 
#' see: http://snpeff.sourceforge.net/SnpEff_manual.html#eff
#' @param df a dataframe with an EFF column containing and EFF string
#' @return a dataframe with columns for 
effParser <- function(df) {
	tidyr::extract(
		df, EFF,
		into = c(
			"Effect", "Effect_Impact", "Functional_Class", "Codon_Change",
			"Amino_Acid_Change", "Amino_Acid_Length", "Gene_Name",
			"Transcript_BioType", "Gene_Coding", "Transcript_ID",
			"Exon_Rank", "Genotype_Number", "ERRORS", "WARNINGS"
		),
		regex = paste0(
			"(.+)\\(", # Effect
			"(.*)\\|(.*)\\|(\\d*)\\|(.*)\\|(\\d*)\\|(.*)",
			"\\|(.*)\\|(.*)\\|(.*)\\|(.*)\\|(.*)",
			"(?:\\|(.*))?(?:\\|(.*))?", # Errors and Warnings (optional)
			"\\)"
		)#,
		#convert = TRUE
	)
}

#' effExtractor
#' split multiple EFF strings into separate individual EFF strings
#'  and convert each to a row in a table
#' @param eff an EFF string or multiple comma separated EFF strings
#' @return a tibble with an EFF column where each row is a single EFF string
effExtractor <- function(eff) {
	effParser(
		tibble::tibble(
			EFF = strsplit(eff, ",")[[1]]
		)
	)
}


#' gtformfun
#' @param dt a DT data table 
#' @param id the name of the column to colour by genotype
gtformfun <- function(dt, id) {
	dt %>%
		DT::formatStyle(
			id, id,
			backgroundColor = DT::styleEqual(
				c("0/0", "0/1", "1/1"),
				c("#458B00","#FF7F00", "#CD3333")
			),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		)
}

## | Quality ------------------------------------------------------------------
highDepth <- function(df) {
	df %>%
		group_by(CHROM) %>%
		mutate(
			meanDepth = mean(DP),
			if_else(
				(mean(DP)*3) > DP, paste0(FILTER, "highDepth"), FILTER
			)
		) %>%
		ungroup()
}

# strand biases?
# mapping quality

## | Multi genotype locus handeling -------------------------------------------
splitGeno <- function(dfr) {
	spltcls <- dfr %>% 
		as.list() %>% 
		purrr::map(~{ # only split character vectors - keeps types consistent later
			if(is.character(.x)){
				strsplit(.x, ",")[[1]]
			} else {
				.x
			}
		})
	EFF <- character(length = length(spltcls$ALT))
	# genotypes alternate in the eff string so assigning every Nth value to genotype N
	for (i in seq_along(spltcls$ALT)) {
		EFF[i] <- paste0(
			spltcls$EFF[seq(i, length(spltcls$EFF), by = length(spltcls$ALT))],
			collapse = ","
		)
	}
	spltcls$EFF <- EFF
	spltcls %>% tibble::as_tibble()
}

lociByGenotype <- function(df) {
	nms <- colnames(df)
	gts <- nms[grepl("gt_", nms)]
	gts <- gts[gts != "gt_GT"]
	
	vcftidy_wide <- df %>%
		#mutate(POS = as.integer(POS)) %>%
		tidyr::unite(pos, CHROM, POS, sep = "_", remove = FALSE) %>%
		dplyr::select(-all_of(gts)) %>%
		tidyr::pivot_wider(
			names_from = Indiv, values_from = gt_GT
		)
	
	vcftidy_wide %>%
		dplyr::filter(NUMALT > 1) %>%
		dplyr::group_by(pos) %>%
		dplyr::group_split() %>%
		purrr::map_dfr(splitGeno) %>%
		#readr::type_convert() %>%
		mutate(
			POS = as.integer(POS),
			QUAL = as.double(QUAL),
			AC = as.character(AC)
		) %>%
		dplyr::bind_rows(vcftidy_wide %>% dplyr::filter(NUMALT <= 1)) %>%
		dplyr::mutate(AF = as.numeric(AF)) %>%
		dplyr::arrange(desc(AF), desc(QUAL), CHROM, POS)
}

## | Deletion filtering -------------------------------------------------------
categoriseSample <- function(x, nm) {
	nm <- gsub("sample_(.*)", "\\1", nm)
	tb <- NULL
	if (any(c("1/1", "0/1") %in% x)) {
		tb <- tibble::tibble_row(sample = nm, keep = TRUE)
	} else if("0/0" %in% x) {
		tb <- tibble::tibble_row(sample = nm, keep = FALSE)
	}
	return(tb)
}

delFeatFun <- function(df, sampids, lst) {
	#lst <- input[paste0("sample_", sampids)]
	if(all(sapply(lst, is.null))) { return(df) }
	samples2keep <- purrr::map2_dfr(lst, names(lst), ~categoriseSample(.x, .y))
	df <- df %>% filter(
		dplyr::across(
			dplyr::matches(
				paste0(
					"score_",
					dplyr::filter(samples2keep, keep == TRUE)$sample
				)
			),
			~ !is.na(.x)
		),
		dplyr::across(
			dplyr::matches(
				paste0(
					"score_",
					dplyr::filter(samples2keep, keep == FALSE)$sample
				)
			),
			~is.na(.x)
		)
	)
	return(df)
}

## | Plots --------------------------------------------------------------------
lociPlot <- function(df) { ## var_type_colours !! global
	ggplot2::ggplot(df, ggplot2::aes(POS, AF)) + 
		ggplot2::geom_point(ggplot2::aes(colour = TYPE, alpha = QUAL)) + 
		#geom_point_interactive(aes(colour = TYPE, alpha = QUAL)) + 
		ggplot2::facet_wrap(~CHROM, nrow = 1, scales = "free_x") + 
		ggplot2::scale_x_continuous(labels = scales::comma) +
		ggplot2::scale_colour_manual(values = var_type_colours) + 
		ggplot2::theme_light() +
		ggplot2::theme(
			axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
		) + 
		ggplot2::labs(
			x = "Position (bp)",
			y = "Allele Frequency"
		)
}


# mutTypeFreqPlot
mutTypeFreqPlot <- function(df) { ## var_type_colours !! global
	df %>%
		dplyr::count(TYPE) %>% 
		dplyr::mutate(
			pc = (n/sum(n)) * 100
		) %>% { 
			ggplot2::ggplot(., aes(TYPE, n)) + 
				ggplot2::geom_col(
					ggplot2::aes(
						fill = TYPE,
						text = sprintf(
							"%.1f%% %s (N = %s)",
							pc, TYPE, format(n, big.mark = ",")
						)
					)
				) + 
				ggplot2::scale_fill_manual(values = var_type_colours) + 
				ggplot2::theme_minimal() + 
				ggplot2::theme(legend.position = "bottom") + 
				ggplot2::labs(x = "Type", y = "Number of Mutations")
		} %>%
		plotly::ggplotly(tooltip = "text")
}

## | UI generators ------------------------------------------------------------
alias_samples <- function(samples) {
	renamers <- lapply(samples, function(sampid) {
		textInput(
			inputId = paste0("alias_", sampid),
			label = "Sample Name",
			placeholder = sampid
		)
	})
	tagList(renamers)
}
