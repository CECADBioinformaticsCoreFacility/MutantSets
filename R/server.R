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
	
	sample_vars_tolisten <- reactive({
		lapply(samples(), function(x){ input[[paste0("sample_", x)]] })
	})
	
	chrs <- reactive({
		req(input$vcf)
		unique(vcftidy()$dat$CHROM)
	})
	# Genotype filtering ----------------------------------------------------------
	## set sample names
	output$renameUI <- renderUI({
		req(sample_vars_tolisten())
		alias_samples(samples())
	})
	
	# Genotype filters UI
	output$setSelector <- renderUI({
		genotype_selector(samples(), input)
	})
	
	# multiVarLoci <- reactive({
	# 	vcftidy()$dat %>%
	# 		filter(NUMALT > 1)
	# 		group_by(pos) %>%
	# 		group_split() %>%
	# 		map_dfr(splitGeno)
	# })
	
	allLoci <- reactive({
		lociByGenotype(vcftidy()$dat)
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
	
	output$total_vb <- renderValueBox({
		req(input$vcf)
		valueBox(
			#width = 8,
			format(nrow(allLoci()), big.mark = ","),
			"Total Variants", icon = icon("dna"),
			color = "purple"
		)
	})
	
	output$nfiltered_vb <- renderValueBox({
		req(input$vcf)
		valueBox(
			#width = 8,
			format(nrow(loci()), big.mark = ","),
			"Filtered Variants", icon = icon("dna")
		)
	})
	

	filtfun <- eventReactive(sample_vars_tolisten(),{
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
	
	## Allele Frequency plot --------------------------------------------------
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
	
	output$mutTypeFreqPlot <- renderPlotly({
		mutTypeFreqPlot(loci())
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
	## Variant effect annotation --------------------------------
	output$effect <- DT::renderDataTable({
		req(input$filtVarsDT_row_last_clicked)
		gen_var_eff_DT(loci(), input$filtVarsDT_row_last_clicked)
	})
	
	## vars tab -------------------------
	output$col_picker <- renderUI({
		variant_column_selector(vcftidy()$meta)
	})
	
	output$filtVarsDT <- DT::renderDataTable({
		gen_var_tab(loci(), samples())
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
	
	genotypeFilterInputs <- eventReactive(sample_vars_tolisten(), {
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