enableBookmarking(store = "server")
# Global data -----------------------------------------------------------------
# 
gtfnames <- c(
	"sequence", "source", "feature", "start", "end",
	"score", "strand","phase","attributes"
)

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

#' gen_var_eff_DT
#' Generates the datatable of predicted variant effects
#' @param loci a data.frame / tibble of loci
#' @param row_clicked the number of a row in the loci table
#' @return a DT data table of the predicted effects of the variant
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

#' wormbase_view
#' Generates an HTML embed tag for a wormbase 'jbrowse' session at the
#' specified location, optinally with a region highlighted.
#' @param chr the chromosome (no preceeding chr just the numeral)
#' @param vstart view starting position
#' @param vend view ending position
#' @param hstart highlight starting position (optional)
#' @param hend highlight ending position (optional)
#' @return character string of html embed tag for jbrowe at locus
wormbase_view <- function(chr, vstart, vend, hstart = NULL, hend = NULL) {
	# NB chr form: "III"	
	stopifnot(chr %in% c(as.character(utils::as.roman(1:6)),"M","X"))
	url <- paste0("https://wormbase.org/tools/genome/jbrowse-simple/",
		"?data=data%2Fc_elegans_PRJNA13758&",
		"loc=", chr, "%3A", vstart, "..", vend,
		"&tracks=Curated_Genes%2CClassical_alleles" #"%2CYACs_Fosmids_Cosmids"
	)
	
	if(!is.null(hstart)) {
		url <- paste0(url,"&highlight=", chr, "%3A", hstart, "..", hend)
	}
	#url
	paste0("<embed width = '100%' height = '500px' src='", url, "'>")
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
		"DP_filter", "Depth (DP) default max: 3x mean",
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

AF_slider <- function() {
	sliderInput(
		"AF_filter", "Variant Allele Frequence (AF)",
		min = 0, max = 1,
		value = 0.8, step = 0.01
	)
}

quality_sliders <- function(allloci) {
	tagList(
		DP_slider(allloci),
		QUAL_slider(allloci),
		QR_slider(allloci),
		QA_slider(allloci),
		AF_slider()
	)
}

quality_filters <- function(df, input) {
	df %>% dplyr::filter(
		DP >= input$DP_filter[1], DP <= input$DP_filter[2],
		QUAL >= input$QUAL_filter,
		QR >= input$QR_filter,
		QA >= input$QA_filter,
		AF >= input$AF_filter
	)
}

# strand biases?
# mapping quality

## | Multi genotype locus handeling -------------------------------------------

#' split_geno
#' 
#' splits the comma separated fields of multivariate loci, and puts them on
#' separate rows in a new tibble.
#' The EFF field requires special treatment as there is not a 1-1 mapping
#' between number of fields and number of variants, it is split and recombined
#' appropriately.
#' 
#' @param dfr a single row from a data.frame / tibble
#' @return A tibble 
split_geno <- function(dfr) {
	spltcls <- dfr %>% 
		as.list() %>% 
		# only split character vectors - keeps types consistent later
		purrr::map(~{ 
			if(is.character(.x)){
				strsplit(.x, ",")[[1]]
			} else {
				.x
			}
		})
	EFF <- character(length = length(spltcls$ALT))
	# genotypes alternate in the eff string so
	# assigning every Nth value to genotype N ~ __CHECK DIFFERING NUM EFFs__
	for (i in seq_along(spltcls$ALT)) {
		EFF[i] <- paste0(
			spltcls$EFF[seq(i, length(spltcls$EFF), by = length(spltcls$ALT))],
			collapse = ","
		)
	}
	spltcls$EFF <- EFF
	spltcls %>% tibble::as_tibble()
}

#' loci_by_genotype
#' 
#' Spreads genotype information into wider format for genotype filtering
#' splits loci with more than 1 variant into separate rows.
#' 
#' @param df a data.frame / tibble for loci
#' @return A tibble with the genotypes spread and the multivariant loci
#' split into separate rows.
loci_by_genotype <- function(df) {
	nms <- colnames(df)
	gts <- nms[grepl("gt_", nms)]
	gts <- gts[gts != "gt_GT"]
	
	vcftidy_wide <- df %>%
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
			QUAL = as.double(QUAL),
			AC = as.character(AC)
		) %>%
		dplyr::bind_rows(vcftidy_wide %>% dplyr::filter(NUMALT <= 1)) %>%
		dplyr::mutate(
			AF = as.numeric(AF),
			QR = as.integer(QR),
			QA = as.integer(QA),
		) %>%
		dplyr::arrange(desc(AF), desc(QUAL), CHROM, POS)
}

#' get_genotype_class
#' 
#' Classifies genotypes into Homozygous Reference of Alternate alleles or
#' Heterozygous reference (one allele is reference) and Heterozygous alternate 
#' (both alleles are non-reference and are different from one-another)
#' 
#' @param sig genotype signature of the form '0/0'
#' @return a character string with the class of genotype
get_genotype_class <- function(sig) {
	if( (length(sig) != 1) | (!is.character(sig)) ) {
		stop("sig must be a character vector of length 1")
	}
	if(is.na(sig)) {
		return(as.character(NA))
	}
	
	splitg <- strsplit(sig,"/")[[1]]
	if(splitg[1] == splitg[2]) {
		if(splitg[1] == "0") {
			return("HomoRef")
		} else {
			return("HomoAlt")
		}
	} else {
		if(splitg[1] != "0") {
			return("HeteroAlt")
		} else {
			return("HeteroRef")
		}
	}
}

## | Deletion filtering -------------------------------------------------------

#' import_deletions
#' @param filenm path to the .gff deletions file generated by MiModD
#' @return a data.frame of deletion variants
import_deletions <- function(filenm) {
	vroom::vroom(
		filenm, col_names = gtfnames, delim = "\t", col_types = "ccciiiccc"
	) %>%
		as.data.frame() %>% 
		tidyr::extract(
			attributes,
			into = c("sample","p_value"),
			regex = "sample=(\\d+);p_value=(-?[\\d.]+(?:e-?\\d+)?)",
			convert = TRUE
		) %>%
		dplyr::select(-source, -strand, -phase)
}

get_del_gen_class <- function(sig) {
	if(is.na(sig)) {
		return("HomoRef")
	} else {
		return("HomoAlt") #,"HeteroAlt")
	}
}

## | Plots --------------------------------------------------------------------

#' loci_plot
#' 
#' @param df plot data
loci_plot <- function(df) { ## var_type_colours !! global
	var_count_by_chr <- df %>% 
		dplyr::group_by(CHROM) %>%
		dplyr::summarise(n = dplyr::n())
	mylabs <- paste0(
		var_count_by_chr$CHROM, " (",
		format(var_count_by_chr$n, big.mark = ","),")"
	)
	names(mylabs) <- var_count_by_chr$CHROM
	ggplot2::ggplot(df, ggplot2::aes(POS/1e6, AF)) + 
		#ggplot2::geom_point(ggplot2::aes(colour = TYPE, alpha = QUAL)) + 
		ggplot2::geom_point(ggplot2::aes(colour = TYPE)) + 
		#geom_point_interactive(aes(colour = TYPE, alpha = QUAL)) + 
		ggplot2::facet_wrap(
			~CHROM, nrow = 1, scales = "free_x",
			labeller = ggplot2::as_labeller(mylabs)
		) + 
		ggplot2::scale_x_continuous(labels = scales::comma) +
		ggplot2::scale_colour_manual(values = var_type_colours) + 
		#ggplot2::theme_light() +
		ggplot2::theme_bw() +
		ggplot2::lims(y = c(0, 1)) +
		ggplot2::theme(
			axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
		) + 
		ggplot2::labs(
			x = "Position (Mbp)",
			y = "Allele Frequency",
			colour = "", alpha = ""
		)
}

## | SNP_freq_plot ------------------------------------------------------------

SNP_freq_plot <- function(df) {
	ggplot2::ggplot(df, ggplot2::aes(POS)) +
		ggplot2::geom_density() +
		ggplot2::facet_wrap(~CHROM, nrow = 1, scales = "free_x") + 
		ggplot2::scale_x_continuous(labels = scales::comma) +
		ggplot2::theme_light() +
		ggplot2::theme(
			axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
		) + 
		ggplot2::labs(
			x = "Position (bp)",
			y = "Variant Density"#,colour = "", alpha = ""
		)
}

#' layout_ggplotly
#' 
#' Tweaks the layout of the x and y axis labels so they don't overlap
#' the axis ticks
#' from: https://stackoverflow.com/questions/42763280/r-ggplot-and-plotly-axis-margin-wont-change
#' 
#' @param gg and plotly plot object
#' @param x a number by which to adjust the vertical position of the x axis name
#' @param y a number by which to adjust the horizontal position
#'  of the y axis name
#' @return a plotly plot object
layout_ggplotly <- function(gg, x = -0.1, y = -0.04){
	# The 1 and 2 goes into the list that contains the options for
	# the x and y axis labels respectively
	gg[['x']][['layout']][['annotations']][[1]][['y']] <- x
	gg[['x']][['layout']][['annotations']][[2]][['x']] <- y
	gg
}

#' mut_type_freq_plot
#' 
#' A plot of the number of variant of a given type are currently filtered
#' 
#' @param df plot data
#' @return a plotly plot object
mut_type_freq_plot <- function(df) { ## var_type_colours !! global
	df %>%
		dplyr::count(TYPE) %>% 
		dplyr::mutate(
			pc = (n/sum(n)) * 100
		) %>% { 
			ggplot2::ggplot(., ggplot2::aes(TYPE, n)) + 
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
		plotly::ggplotly(tooltip = "text") %>%
		plotly::config(
			displaylogo = FALSE,
			modeBarButtonsToRemove = list(
				"autoScale2d", "resetScale2d",
				"hoverClosestCartesian", "hoverCompareCartesian",
				"select2d", "lasso2d", "zoomIn2d", "zoomOut2d",
				"toggleSpikelines"
			)
		)
}

## | UI generators ------------------------------------------------------------

#' alias_samples
#' 
#' Generates the text input UI elements for aliasing the samples
#' @param samples Vector of sample identifiers
#' @return tagList of textInput elements
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

#' genotype_selector
#' 
#' generates genotype filtering UI elements
#' @param samples Vector of sample identifiers
#' @param input the shiny input object
#' @return UI for set subtraction a tagList of shinyWidgets checkboxGroupButtons
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

#' variant_column_selector
#' 
#' generates the column selector UI element for the main variant table
#' @param meta data.frame column details (vcfR2tidy metadata table)
#' @return Column selector UI element
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

#' chr_selectizer
#' @param chr a charater vector of chromosome names
#' @return a shinyWidgets pickerInput UI element that lets the user select
#' which chromosomes to display
chr_selectizer <- function(chr) {
	shinyWidgets::pickerInput(
		"picked_chr", "Select Chromosomes",
		choices = chr,
		selected = chr,
		multiple = TRUE,
		shinyWidgets::pickerOptions(
			actionsBox = TRUE,
			selectedTextFormat = "count > 3",
		)
	)
}

## | Main Variant Table ----------------------------

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

#' gtformatter
#' 
#' Applies the formatter to each sample column in turn and adds the next
#' formatting rule to the DT object for all samples
#' 
#' @param dt a DTdatatable object
#' @param samples a character vector of sample names
#' @return a function which take a DT datatable object and fills the genotype
#' columns of the DT object with different colours depending on the genotype
gtformatter <- function(dt, samples) {
	Reduce(gtformfun, samples, init = dt)
}

#' barplotformfun
#' @param dt a DT data table 
#' @param id the name of the column to colour by genotype
#' @param df data frame for scale values
barplotformfun <- function(dt, id, df) {
	dt %>%
		DT::formatStyle(
			id, id,
			background = DT::styleColorBar(
				df[[id]],
				'lightblue'
			),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		)
}

#' barplotformatter
#' 
#' Applies the formatter to each sample column in turn and adds the next
#' formatting rule to the DT object for all samples
#' 
#' @param dt a DTdatatable object
#' @param vars a character vector of variables
#' @param df initial object
barplotformatter <- function(dt, vars, df) {
	purrr::reduce(vars, barplotformfun, df, .init = dt)
}

#' gen_var_tab
#' 
#' Formatting the main variant table
#' 
#' @param loci data.frame / tibble of loci
#' @param samples Vector of sample identifiers
#' @param aliases Vector of sample aliases
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
		barplotformatter(vars = c("DP", "QUAL", "QR", "QA"), df = df) %>%
		DT::formatStyle(
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

#' gen_del_tab
#' 
#' Formatting the deletion variant table
#' 
#' @param deltab data.frame / tibble of loci
#' @param samples Vector of sample identifiers
#' @param aliases Vector of sample aliases
gen_del_tab <- function(deltab, samples, aliases) {
	aliases[aliases == ""] <- samples[aliases == ""]
	cnms <- colnames(deltab)
	purrr::walk2(
		samples, aliases, ~{
			cnms <<- gsub(
				paste0("(.*_)", .x ,""),
				paste0("\\1",.y),
				cnms
			)
		}
	)
	deltab %>%
		DT::datatable(
			colnames = cnms,
			rownames = FALSE, selection = "single",
			options = list(
				scrollX = TRUE
			)
		)
}
