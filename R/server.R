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
		vcftidy()$dat %>% dplyr::distinct(Indiv) %>% dplyr::pull()
	})
	
	chrs <- reactive({
		req(input$vcf)
		unique(vcftidy()$dat$CHROM)
	})
	# Genotype filtering ----------------------------------------------------------
	## set sample names
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
	
	# Genotype filters UI
	output$setSelector <- renderUI({
		selector <- lapply(samples(), function(sampid) {
			label <- paste0("Genotypes to Include for ", sampid)
			if(!is.null(input[[paste0("alias_", sampid)]])) { # input is not null prior to user setting values?
				label <- paste0(
					"Genotypes to Include for ", sampid,
					" (", input[[paste0("alias_", sampid)]], ")" # bug parens are appearing without content
				)
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
	})
	
	# multiVarLoci <- reactive({
	# 	vcftidy()$dat %>%
	# 		filter(NUMALT > 1)
	# 		group_by(pos) %>%
	# 		group_split() %>%
	# 		map_dfr(splitGeno)
	# })
	
	allLoci <- reactive({
		nms <- colnames(vcftidy()$dat)
		gts <- nms[grepl("gt_", nms)]
		gts <- gts[gts != "gt_GT"]
		
		vcftidy_wide <- vcftidy()$dat %>%
			###!!!! Need to Handle Multiple Values here!!!
			#filter(NUMALT <= 1) %>% 
			#mutate(AF = as.numeric(AF)) %>%
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
			dplyr::bind_rows(vcftidy_wide %>% dplyr::filter(NUMALT <= 1)) %>%
			dplyr::mutate(AF = as.numeric(AF)) %>%
			dplyr::arrange(AF, QUAL, CHROM, POS)
		
		# vcftidy()$dat %>%
		# 	###!!!! Need to Handle Multiple Values here!!!
		# 	filter(NUMALT <= 1) %>% 
		# 	mutate(AF = as.numeric(AF)) %>%
		# 	tidyr::unite(pos, CHROM, POS, sep = "_", remove = FALSE) %>%
		# 	select(-all_of(gts)) %>%
		# 	tidyr::pivot_wider(
		# 		names_from = Indiv, values_from = gt_GT
		# 	)
	})
	
	locusGenoTypes <- reactive({
		req(input$pick_cols)
		allLoci() %>%
			dplyr::select(
				pos, CHROM, POS, REF, ALT, QUAL,
				# DP, QR, ODDS, TYPE, NUMALT, EFF, AF,
				all_of(input$pick_cols),
				EFF,
				all_of(samples())
			)
	})
	
	tolisten <- reactive({
		#lapply(samples(), function(x){ paste0("sample_", x) })
		lapply(samples(), function(x){ input[[paste0("sample_", x)]] })
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
		
		#dplyr::pull(pos)
	})
	
	# |Download data ---------------------------------------------------------
	output$downloadData <- downloadHandler(
		filename = "selected_mutants.tsv",
		content = function(file) {
			vroom::vroom_write(loci(), file)
		}
	)
	
	## Chr plot ----------------------------------------
	#output$testpoints <- renderPrint(print(input$chrplotclick))
	#output$chrplot <- renderPlot({
	output$chrplot <- plotly::renderPlotly({
		#output$chrplot <- renderGirafe({
		plot <- loci() %>%
			lociPlot()
		#plot
		plotly::ggplotly(plot)
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
		#req(input$filtVarsDT_rows_selected)
		df <- loci() %>% 
			dplyr::slice(input$filtVarsDT_row_last_clicked) %>%
			#dplyr::slice(input$filtVarsDT_rows_selected) %>%
			dplyr::pull(EFF) %>%
			effExtractor() %>%
			dplyr::select(-Amino_Acid_Length, -ERRORS, -WARNINGS, -Genotype_Number)
		
	})
	## eff --------------------------------
	output$effect <- DT::renderDataTable({
		req(input$filtVarsDT_row_last_clicked)
		df <- loci() %>% 
			dplyr::slice(input$filtVarsDT_row_last_clicked) %>%
			dplyr::pull(EFF) %>%
			effExtractor() %>%
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
		shinyWidgets::pickerInput(
			inputId = "pick_cols", 
			label = "Select/deselect all Columns", 
			choices = optlabeled, 
			selected = c(
				"CHROM", "POS", "REF", "ALT", "QUAL", "DP", "QR", "ODDS",
				"TYPE", "NUMALT", "EFF", "AF"
			),
			options = shinyWidgets::pickerOptions(
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
		df <- loci() %>% dplyr::select(-EFF, -pos)
		brks <- quantile(df$AF, probs = seq(.05, .95, .05), na.rm = TRUE)
		clrs <- scales::colour_ramp(
			#c("#ffeda0", "#feb24c", "#f03b20")
			c("#efedf5", "#bcbddc", "#756bb1")
		)(c(0, brks))
		
		clrs <- df %>% 
			DT::datatable(
				rownames = FALSE,
				selection = "single",
				options = list(
					#autoWidth = TRUE,
					#scrollX=TRUE,
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
				backgroundColor = DT::styleInterval(
					brks, clrs
				)
			) %>%
			gtformatter()()
	}, server = TRUE)
	
	# Deletions ---------------------------------------------------------------
	
	gff <- reactive({
		req(input$gff)
		vroom::vroom(input$gff$datapath, col_names = gtfnames, delim = "\t") %>%
			as.data.frame() %>% 
			tidyr::extract(
				attributes,
				into = c("sample","p_value"),
				regex = "sample=(\\d+);p_value=(-?[\\d.]+(?:e-?\\d+)?)",
				convert = TRUE
			) %>%
			dplyr::select(-source, -strand, -phase)
	})
	
	gff_wide <- reactive({
		gff() %>%
			dplyr::group_by(sequence, feature, start, end) %>%
			tidyr::pivot_wider(
				names_from = sample, values_from = c("score", "p_value")
			)
	})
	
	genotypeFilterInputs <- eventReactive(tolisten(), {
		lst <- map(samples(), ~input[[paste0("sample_", .x)]])
		names(lst) <- paste0("sample_", samples())
		lst
	})
	
	delTab <- reactive({
		delFeatFun(gff_wide(), samples(), genotypeFilterInputs())
	})
	
	output$filteredDels <- DT::renderDataTable({
		#gff() %>% DT::datatable(rownames = FALSE)
		delTab() %>% DT::datatable(rownames = FALSE)
	})
	
	# output$qualplot <- renderPlot({
	# 	qualityscores() %>%
	# 		ggplot(aes(CHROM, scores)) +
	# 		lvplot::geom_lv(aes(fill = CHROM)) +
	# 		facet_wrap(~metric, scales = "free", ncol = 1)
	# })
}
