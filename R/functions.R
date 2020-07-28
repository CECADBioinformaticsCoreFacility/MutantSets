
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

#' eff_parser
#' parsing snpeff strings EFF field 
#' see: http://snpeff.sourceforge.net/SnpEff_manual.html#eff
#' @param df a dataframe with an EFF column containing and EFF string
#' @return a dataframe with columns for 
eff_parser <- function(df) {
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

#' eff_extractor
#' split multiple EFF strings into separate individual EFF strings
#'  and convert each to a row in a table
#' @param eff an EFF string or multiple comma separated EFF strings
#' @return a tibble with an EFF column where each row is a single EFF string
eff_extractor <- function(eff) {
	eff_parser(
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

gen_var_eff_DT <- function(loci, row_clicked) {
	df <- loci %>% 
		dplyr::slice(row_clicked) %>%
		dplyr::pull(EFF) %>%
		eff_extractor() %>%
		dplyr::select(-Amino_Acid_Length, -ERRORS, -WARNINGS, -Genotype_Number)
	
	DT::datatable(
		df, rownames = FALSE, options = list(scrollX = TRUE),
	) %>%
		DT::formatStyle(
			"Effect", "Effect_Impact", 
			backgroundColor = DT::styleEqual(
				toupper(c("High", "Moderate", "Low", "Modifier")),
				c("#CD3333", "#FF7F00", "#458B00", "#969696")
			),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		)
}

## | Quality ------------------------------------------------------------------
# high_depth <- function(df) {
# 	df %>%
# 		group_by(CHROM) %>%
# 		mutate(
# 			meanDepth = mean(DP),
# 			if_else(
# 				(mean(DP)*3) > DP, paste0(FILTER, "highDepth"), FILTER
# 			)
# 		) %>%
# 		ungroup()
# }

DP_slider <- function(allloci) {
	rng <- range(allloci$DP, na.rm = TRUE)
	sliderInput(
		"DP_filter", "Depth (DP)",
		min = rng[1], max = rng[2],
		value = c(rng[1], mean(allloci$DP) * 3),
		step = 1
	)
}

QUAL_slider <- function(allloci) {
	rng <- range(allloci$QUAL, na.rm = TRUE)
	sliderInput(
		"QUAL_filter", "Variant Quality (QUAL)",
		min = rng[1], max = rng[2],
		value = 20, step = 1
	)
}

QR_slider <- function(allloci) {
	rng <- range(allloci$QR, na.rm = TRUE)
	sliderInput(
		"QR_filter", "Reference Allele Quality (QR)",
		min = rng[1], max = rng[2],
		value = 20, step = 1
	)
}

QA_slider <- function(allloci) {
	rng <- range(allloci$QA, na.rm = TRUE)
	sliderInput(
		"QA_filter", "Alternate Allele Quality (QA)",
		min = rng[1], max = rng[2],
		value = 20, step = 1
	)
}

quality_sliders <- function(allloci) {
	tagList(
		DP_slider(allloci),
		QUAL_slider(allloci),
		QR_slider(allloci),
		QA_slider(allloci)
	)
}

quality_filters <- function(df, input) {
	df %>% dplyr::filter(
		DP >= input$DP_filter[1], DP <= input$DP_filter[2],
		QUAL >= input$QUAL_filter,
		QR >= input$QR_filter,
		QA >= input$QA_filter,
	)
}

# strand biases?
# mapping quality

## | Multi genotype locus handeling -------------------------------------------
split_geno <- function(dfr) {
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

loci_by_genotype <- function(df) {
	nms <- colnames(df)
	gts <- nms[grepl("gt_", nms)]
	gts <- gts[gts != "gt_GT"]
	
	vcftidy_wide <- df %>%
		mutate(QA = as.integer(QA)) %>%
		tidyr::unite(pos, CHROM, POS, sep = "_", remove = FALSE) %>%
		dplyr::select(-all_of(gts)) %>%
		tidyr::pivot_wider(
			names_from = Indiv, values_from = gt_GT
		)
	
	vcftidy_wide %>%
		dplyr::filter(NUMALT > 1) %>%
		dplyr::group_by(pos) %>%
		dplyr::group_split() %>%
		purrr::map_dfr(split_geno) %>%
		#readr::type_convert() %>%
		mutate(
			POS = as.integer(POS),
			QA = as.integer(QA),
			QUAL = as.double(QUAL),
			AC = as.character(AC)
		) %>%
		dplyr::bind_rows(vcftidy_wide %>% dplyr::filter(NUMALT <= 1)) %>%
		dplyr::mutate(AF = as.numeric(AF)) %>%
		dplyr::arrange(desc(AF), desc(QUAL), CHROM, POS)
}

## | Deletion filtering -------------------------------------------------------
import_deletions <- function(filenm) {
	vroom::vroom(filenm, col_names = gtfnames, delim = "\t") %>%
		as.data.frame() %>% 
		tidyr::extract(
			attributes,
			into = c("sample","p_value"),
			regex = "sample=(\\d+);p_value=(-?[\\d.]+(?:e-?\\d+)?)",
			convert = TRUE
		) %>%
		dplyr::select(-source, -strand, -phase)
}

categorise_sample <- function(x, nm) {
	nm <- gsub("sample_(.*)", "\\1", nm)
	tb <- NULL
	if (any(c("1/1", "0/1") %in% x)) {
		tb <- tibble::tibble_row(sample = nm, keep = TRUE)
	} else if("0/0" %in% x) {
		tb <- tibble::tibble_row(sample = nm, keep = FALSE)
	}
	return(tb)
}

del_feat_fun <- function(df, sampids, lst) {
	#lst <- input[paste0("sample_", sampids)]
	if(all(sapply(lst, is.null))) { return(df) }
	samples2keep <- purrr::map2_dfr(lst, names(lst), ~categorise_sample(.x, .y))
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
loci_plot <- function(df) { ## var_type_colours !! global
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


# mut_type_freq_plot
mut_type_freq_plot <- function(df) { ## var_type_colours !! global
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

genotype_selector <- function(samples, input) {
	selector <- lapply(samples, function(sampid) {
		label <- paste0("Genotypes to Include for ", sampid)
		alias <- input[[paste0("alias_", sampid)]]
		if(!is.null(alias)) {
			if(alias != "") {
				label <- paste0(
					"Genotypes to Include for ", sampid,
					" (", alias, ")"
				)
			}
		}
		shinyWidgets::checkboxGroupButtons(
			inputId = paste0("sample_", sampid),
			label = label, 
			choices = c("0/0", "0/1", "1/1"), 
			justified = TRUE, status = "primary",
			checkIcon = list(
				yes = icon("ok", lib = "glyphicon"),
				no = icon("remove", lib = "glyphicon")
			)
		)
	})
	tagList(selector)
}

variant_column_selector <- function(meta){
	opts <- meta %>%
		dplyr::filter(!grepl("gt_", ID)) %>%
		dplyr::filter(ID != "EFF")
	
	optlabeled <- opts$ID
	
	names(optlabeled) <- paste0("(", opts$ID, ") ", opts$Description)
	
	shinyWidgets::pickerInput(
		inputId = "pick_cols", 
		label = "Select/deselect all Columns", 
		choices = optlabeled, 
		selected = c( # factor out as default col variable?
			"CHROM", "POS", "REF", "ALT", "QUAL", "DP", "QR", "QA", "ODDS",
			"TYPE", "NUMALT", "EFF", "AF"
		),
		options = shinyWidgets::pickerOptions(
			actionsBox = TRUE,
			selectedTextFormat = "count > 3",
			liveSearch = TRUE,
			size = 10
		),
		multiple = TRUE
	)
}

## | Main Variant Table ----------------------------

gtformatter <- function(dt, samples) {
	Reduce(gtformfun, samples, init = dt)
}

gen_var_tab <- function(loci, samples, aliases) {
	# don't display the long EFF column or the unique position identifier
	df <- loci %>% dplyr::select(-EFF, -pos)
	
	# colour scale for allele frequency column
	brks <- quantile(df$AF, probs = seq(.05, .95, .05), na.rm = TRUE)
	clrs <- scales::colour_ramp(
		#c("#ffeda0", "#feb24c", "#f03b20")
		c("#efedf5", "#bcbddc", "#756bb1")
	)(c(0, brks))
	
	
	aliases[aliases == ""] <- samples[aliases == ""]
	cnms <- colnames(df)
	cnms[cnms %in% samples] <- aliases
	
	dt <- df %>% 
		DT::datatable(
			rownames = FALSE,
			selection = "single",
			colnames = cnms,
			options = list(
				#autoWidth = TRUE,
				scrollX = TRUE,
				columnDefs = list(
					list(
						render = htmlwidgets::JS(
							"function(data, type, row, meta) {",
							"return type === 'display' && data.length > 10 ?",
							"'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
							"}"),
						targets = list(2, 3)
					)
				)#"_all"
			)
		) %>%
		## barplots ~ reactive to column selection ~ errors when deselecting default columns
		DT::formatStyle(
			'DP','DP',
			background = DT::styleColorBar(
				df[["DP"]],
				'lightblue'
			),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		) %>%
		DT::formatStyle(
			'QUAL','QUAL',
			background = DT::styleColorBar(
				df[["QUAL"]],
				'lightblue'
			),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		) %>% DT::formatStyle(
			"TYPE", "TYPE",
			backgroundColor = DT::styleEqual(
				names(var_type_colours), var_type_colours
			)
		) %>% 
		DT::formatStyle(
			"AF", "AF",
			backgroundColor = DT::styleInterval(
				brks, clrs
			)
		) %>%
		gtformatter(samples)
	
	return(dt)
}


# ,
# sliderInput("QR", "Reference Quality"),
# sliderInput("QA", "Alternative Quality")
