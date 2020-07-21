#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Notes -----------------------------------------------------------------------
# Quality Pre-filtering
# 
# - [ ] set membership plot - which sets
# - [ ] mutation types proportion
# - [ ] dropdown button selectize chrs?
# - [ ] pick columns to show (use $meta)
# - [ ] show full gt table in Tab box
# - [ ] keys?
# - [ ] wormbase links?
# - [x] name samples - __display bug on null__
# - [x] colour types in table - match plot
# 
# Libraries -------------------------------------------------------------------
library(shiny)
library(shinydashboard)
library(shinybusy)
library(shinyWidgets)
library(shinyjs)
library(purrr)
library(ggplot2)
library(scales)

library(ggiraph)
library(plotly)

library(dplyr)
library(DT)
library(vcfR)
library(future)
#library(BSgenome.Celegans.UCSC.ce11)
#library(TxDb.Celegans.UCSC.ce11.ensGene)

# Global data -----------------------------------------------------------------
# 
gtfnames <- c(
	"sequence", "source", "feature", "start", "end",
	"score", "strand","phase","attributes"
)
# seq <- ape::as.DNAbin(getSeq(BSgenome.Celegans.UCSC.ce11))
# #vcfin <- vcfR::read.vcfR("~/Downloads/Galaxy171-[SnpSift_Filter_on_data_158].vcf", verbose = FALSE)

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
## |Effect Annotation -------------------------------------------------------

#' effParser
#' parsing snpeff strings EFF field http://snpeff.sourceforge.net/SnpEff_manual.html#eff
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
#' split multiple EFF strings into separate individual EFF strings and convert each to a row in a table
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
				c("0/0", "0/1", "1/0", "1/1"),
				c("#458B00","#FF7F00","#FF7F00", "#CD3333")
			),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		)
}

## | Plots --------------------------------------------------------------------
lociPlot <- function(df) {
	ggplot(df, aes(POS, AF)) + 
		geom_point(aes(colour = TYPE, alpha = QUAL)) + 
		#geom_point_interactive(aes(colour = TYPE, alpha = QUAL)) + 
		facet_wrap(~CHROM, nrow = 1, scales = "free_x") + 
		scale_x_continuous(labels = comma) +
		scale_colour_manual(values = var_type_colours) + 
		theme_light() +
		theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
		labs(
			x = "Position (bp)",
			y = "Allele Frequency"
		)
}
# e.g.
# library(future)
# plan(multiprocess, workers = 6) 
# ## vcftidy$dat
# tmp <- vcftidySel %>% 
# 	mutate(EFF = furrr::future_map(EFF,effExtractor), nEFF = map_int(EFF,nrow))

## | Chr prep ------------------------------------------------------------------
chrprep <- function(inobj, min_DP, max_DP, min_MQ, max_MQ, min_QUAL) {
	chromsMsk <- future.apply::future_lapply(inobj, function(x) {
		mskd <- masker(
			x,
			min_DP = min_DP, max_DP = max_DP,
			min_MQ = min_MQ, max_MQ = max_MQ,
			min_QUAL = min_QUAL
		)
		proc.chromR(mskd)
	})
}


# UI --------------------------------------------------------------------------
## |Sidebar -------------------------------------------------------------------

# Define UI for application that draws a histogram
sidebar <- dashboardSidebar(
	sidebarMenu(
		#menuItem("Inputs", tabName = "inputs", icon = icon("upload")),
		fileInput("vcf", "Select a VCF file", accept = ".vcf"),
		fileInput("gff (Deletions))", "Select a gff file",accept = ".gff"),
		menuItem("Options", tabName = "options", icon = icon("th")),
		menuItem("Overview", tabName = "overview", icon = icon("th")),
		menuItem("Filtering", tabName = "table", icon = icon("table"))#,
		#uiOutput("DPslider"),
		#uiOutput("MQslider"),
		#uiOutput("min_QUAL")
	)
)
## |Body ----------------------------------------------------------------------
body <- dashboardBody(
	shinybusy::add_busy_spinner(spin = "fading-circle"),
	tabItems(
		### ||Options ---------------------------------------------------------
		tabItem(
			tabName = "options",
			fluidRow(
				box(
					title = "Sample Aliases",
					uiOutput("renameUI")
				)
			)
		),
		### ||Overview
		tabItem(
			tabName = "overview", # quality filters?
			fluidRow(
				box(
					title = "Quality Ranges", 
					status = "primary", solidHeader = TRUE, collapsible = TRUE,
					width = 4,
					#plotOutput("qualplot")
				),
				box(
					title = "Chr Plots", 
					status = "primary", solidHeader = TRUE, collapsible = TRUE,
					width = 8,
					#uiOutput("chr"),
					#plotOutput("chrPlot")
				)
			)
		),
		### ||Table - ---------------------------------------------------------
		tabItem(
			tabName = "table",
			fluidRow(
				box(
					title = "Genotype Filters", 
					status = "primary", solidHeader = TRUE, collapsible = TRUE,
					width = 4,
					uiOutput("setSelector")
				),
				tabBox(
					title = "Plots", width = 8,
					tabPanel(
						"Allele Frequency",
						useShinyjs(),
						extendShinyjs(
							text = "shinyjs.resetClick = function() { Shiny.onInputChange('.clientValue-plotly_click-A', 'null'); }"
						),
						#plotOutput("chrplot", click = "chrplotclick")#,
						plotlyOutput("chrplot"),
						#verbatimTextOutput("chrplot_sel")
						actionButton("reset","Reset Point selection")#,
						#verbatimTextOutput("chrplot_click")
						#verbatimTextOutput("tmp")
						#girafeOutput("chrplot")
						#verbatimTextOutput("testpoints")
					),
					tabPanel(
						"Effect",
						DTOutput("effect")
					),
					tabPanel(
						"sets"
					),
					tabPanel(
						"Var types"
					)
				)
			),
			fluidRow(
				box(
					title = "Variants", 
					status = "primary", solidHeader = TRUE, collapsible = TRUE,
					width = 12,
					#DT::DTOutput("tidytab")
					#dropdownButton(
						uiOutput("col_picker"),
					#),
					DT::DTOutput("filtVarsDT")
				)
			)
		)
	)
)

## |UI Wrapper ----------------------------------------------------------------
ui <- dashboardPage(
	dashboardHeader(title = "MutantSets"),
	sidebar,
	body
)

# Server ----------------------------------------------------------------------
server <- function(input, output) {
# Options ---------------------------------------------------------------------
	options(shiny.maxRequestSize = Inf) # Do not limit file size
# Input pre-processing --------------------------------------------------------
	vcf <- reactive({ # Parse VCF
		req(input$vcf)
		vcfin <- vcfR::read.vcfR(input$vcf$datapath, verbose = FALSE)
	})
	
	vcftidy <- reactive({ # VCF to tidy dataframe
		vcfR::vcfR2tidy(vcf(), single_frame = TRUE)
	})
	
	samples <- reactive({
		vcftidy()$dat %>% distinct(Indiv) %>% pull()
	})
	
	chrs <- reactive({
		req(input$vcf)
		unique(vcftidy()$dat$CHROM)
	})
# Genotype filtering ----------------------------------------------------------
	## 
	output$renameUI <- renderUI({
		renamers <- lapply(samples(), function(sampid) {
			textInput(
				inputId = paste0("alias_", sampid),
				label = "Sample Name",
				placeholder = sampid
			)
		})
		tagList(renamers)
	})
	
	output$setSelector <- renderUI({
		selector <- lapply(samples(), function(sampid) {
			label <- paste0("Genotypes to Include for ", sampid)
			if(!is.null(input[[paste0("alias_", sampid)]])) {
				label <- paste0(
					"Genotypes to Include for ", sampid,
					" (", input[[paste0("alias_", sampid)]], ")"
				)
			}
			checkboxGroupButtons(
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
	})
	
	allLoci <- reactive({
		nms <- colnames(vcftidy()$dat)
		gts <- nms[grepl("gt_", nms)]
		gts <- gts[gts != "gt_GT"]
		vcftidy()$dat %>%
			filter(NUMALT <= 1) %>% ###!!!! Need to Handle Multiple Values here!!!
			mutate(AF = as.numeric(AF)) %>%
			tidyr::unite(pos, CHROM, POS, sep = "_", remove = FALSE) %>%
			select(-all_of(gts)) %>%
			tidyr::pivot_wider(
				names_from = Indiv, values_from = gt_GT
			)
	})
	
	locusGenoTypes <- reactive({
		req(input$pick_cols)
		allLoci() %>%
			select(
				pos, CHROM, POS, REF, ALT, QUAL,# DP, QR, ODDS, TYPE, NUMALT, EFF, AF,
				all_of(input$pick_cols),
				EFF,
				all_of(samples())
			)
	})
	
	tolisten <- reactive({
		lapply(samples(), function(x){ paste0("sample_", x) })
	})
	
	filtfun <- eventReactive(tolisten(),{
		function(dff, sampids) {
			for (i in seq_along(sampids)) {
				inval <- input[[paste0("sample_", sampids[i])]]
				if(!is.null(inval)) {
					dff <- dff[dff[[sampids[i]]] %in% inval,]
				}
			}
			dff
		}
	})
	
	loci <- reactive({
		# Assuming 1 variant per position need to have handling of multiple vars in one place!!!!!
		# if(is.null(chrplot_click()) | length(chrplot_click()) == 0) {
		# 	df <- locusGenoTypes() %>%
		# 		filtfun()(samples())
		# } else {
		# 	df <- locusGenoTypes() %>%
		# 		slice(chrplot_click()$pointNumber + 1)
		# }
		
		df <- locusGenoTypes() %>%
		filtfun()(samples())

		# if(!is.null(input$chrplotclick)) {
		# 	df <- nearPoints(df, input$chrplotclick, "POS", "AF")
		# }
	
		#pull(pos)
	})
	
	## Chr plot ----------------------------------------
	#output$testpoints <- renderPrint(print(input$chrplotclick))
	#output$chrplot <- renderPlot({
	output$chrplot <- renderPlotly({
	#output$chrplot <- renderGirafe({
		plot <- loci() %>%
			lociPlot()
		#plot
		ggplotly(plot)
		#girafe(code = print(plot))
	})
	
	# output$chrplot_sel <- renderPrint({
	# 	event_data("plotly_selected")
	# })
	# 
	# chrplot_sel <- reactive({
	# 	event_data("plotly_selected")
	# })
	# 
	# chrplot_click <- reactive({
	# 	event_data("plotly_click")
	# })
	
	observeEvent(input$reset, {
		js$resetClick()
	})
	
	#output$chrplot_click <- renderText(names(chrplot_click()))
	
	filtered_vars <- reactive({
		req(input$filtVarsDT_row_last_clicked)
		df <- loci() %>% 
			dplyr::slice(input$filtVarsDT_row_last_clicked) %>%
			pull(EFF) %>%
			effExtractor() %>%
			select(-Amino_Acid_Length, -ERRORS, -WARNINGS, -Genotype_Number)
		
	})
	## eff --------------------------------
	output$effect <- DT::renderDataTable({
		req(input$filtVarsDT_row_last_clicked)
		df <- loci() %>% 
			dplyr::slice(input$filtVarsDT_row_last_clicked) %>%
			pull(EFF) %>%
			effExtractor() %>%
			select(-Amino_Acid_Length, -ERRORS, -WARNINGS, -Genotype_Number)
		
		DT::datatable(df) %>%
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
	})
	
	gtformatter <- reactive({
		function(dt) {
			Reduce(gtformfun, samples(), init = dt)
		}
	})
	## vars tab -------------------------
	output$col_picker <- renderUI({
		opts <- vcftidy()$meta %>%
			dplyr::filter(!grepl("gt_", ID)) %>%
			dplyr::filter(ID != "EFF")
		optlabeled <- opts$ID
		names(optlabeled) <- paste0("(", opts$ID, ") ", opts$Description)
		pickerInput(
			inputId = "pick_cols", 
			label = "Select/deselect all Columns", 
			choices = optlabeled, 
			selected = c(
				"CHROM", "POS", "REF", "ALT", "QUAL", "DP", "QR", "ODDS",
				"TYPE", "NUMALT", "EFF", "AF"
			),
			options = pickerOptions(
				actionsBox = TRUE,
				selectedTextFormat = "count > 3",
				liveSearch = TRUE,
				size = 10
			),
			# options = list(
			# 	`actions-box` = TRUE, 
			# 	size = 10,
			# 	`selected-text-format` = "count > 3"
			# ), 
			multiple = TRUE
		)
	})
	
	output$filtVarsDT <- DT::renderDataTable({
		#df <- filtVars() 
		df <- loci() %>% select(-EFF, -pos)
		brks <- quantile(df$AF, probs = seq(.05, .95, .05), na.rm = TRUE)
		clrs <- colour_ramp(
			#c("#ffeda0", "#feb24c", "#f03b20")
			c("#efedf5", "#bcbddc", "#756bb1")
		)(c(0, brks))
		clrs <- 
		df %>%
			DT::datatable(
				rownames = FALSE,
				options = list(
					#autoWidth = TRUE,
					#scrollX=TRUE,
					columnDefs = list(
						list(
							render = JS(
								"function(data, type, row, meta) {",
								"return type === 'display' && data.length > 10 ?",
								"'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
								"}"),
							targets = list(2, 3)
						)
					)#"_all"
				)
			) %>%
			## barplots
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
				backgroundColor = styleInterval(
					brks, clrs
				)
			) %>%
			gtformatter()()
	}, selection = 'single', server = TRUE)
	
	
	gff <- reactive({
		req(input$gff)
		vroom::vroom(input$gff$datapath, col_names = gtfnames) %>%
		as.data.frame() %>% 
		tidyr::extract(
			attributes,
			into = c("sample","p_value"),
			regex = "sample=(\\d+);p_value=(-?[\\d.]+(?:e-?\\d+)?)",
			convert = TRUE
		)
	})
	
	
	# output$qualplot <- renderPlot({
	# 	qualityscores() %>%
	# 		ggplot(aes(CHROM, scores)) +
	# 		lvplot::geom_lv(aes(fill = CHROM)) +
	# 		facet_wrap(~metric, scales = "free", ncol = 1)
	# })
}

# Run the application ---------------------------------------------------------
shinyApp(ui = ui, server = server)

#chromoqc(do.call("masker",c(chroms[[1]], list(min_QUAL = 1, max_DP = 200)) ),dp.alpha=20)


#%>% select(-source, -strand, -phase) %>% 
# group_by(sequence, feature, start, end) %>%
# tidyr::pivot_wider(names_from = sample, values_from = c("score", "p_value"))
