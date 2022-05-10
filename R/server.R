# Server ----------------------------------------------------------------------
server <- function(input, output) {
	# Options -----------------------------------------------------------------
	options(shiny.maxRequestSize = Inf) # Do not limit file size
	# Input pre-processing ----------------------------------------------------
	## | Parse VCF
	vcf <- reactive({ 
		req(input$vcf)
		vcfin <- vcfR::read.vcfR(input$vcf$datapath, verbose = FALSE)
	})
	
	vcftidy <- reactive({ # VCF to tidy dataframe
		vcfR::vcfR2tidy(vcf(), single_frame = TRUE)
	})
	
	## | Get Basic Data properties --------------------------------------------
	
	samples <- reactive({
		vcftidy()$dat %>% dplyr::distinct(Indiv) %>% dplyr::pull()
	})
	
	sample_vars_tolisten <- reactive({
		lapply(samples(), function(x){ input[[paste0("sample_", x)]] })
	})
	
	aliases <- reactive({
		req(sample_vars_tolisten())
		purrr::map_chr(samples(), ~input[[paste0("alias_", .x)]])
	})
	
	chrs <- reactive({
		req(input$vcf)
		unique(vcftidy()$dat$CHROM)
	})
	
	output$chr_selection <- renderUI({ chr_selectizer(chrs()) })
	outputOptions(output, "chr_selection", suspendWhenHidden = FALSE)
	
	# Transform Genotype Data -------------------------------------------------

	allLoci <- reactive({
		loci_by_genotype(vcftidy()$dat)
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
	
	filtfun <- eventReactive(sample_vars_tolisten(), {
		function(dff, sampids) {
			for (i in seq_along(sampids)) {
				inval <- input[[paste0("sample_", sampids[i])]]
				if(!is.null(inval)) {
					dff <- dff[
						unlist(
							sapply(dff[[sampids[i]]], get_genotype_class)
						) %in% 
						unlist(sapply(inval, get_genotype_class)),
					]
				}
			}
			dff
		}
	})
	
	delfiltfun <- eventReactive(sample_vars_tolisten(), {
		function(dff, sampids) {
			for (i in seq_along(sampids)) {
				inval <- input[[paste0("sample_", sampids[i])]]
				if(!is.null(inval)) {
					dff <- dff[
						unlist(
							sapply(
								dff[[paste0("score_",sampids[i])]],
								get_del_gen_class
							)
						) %in% 
						unlist(sapply(inval, get_genotype_class)),
					]
				}
			}
			dff
		}
	})
	
	loci <- reactive({
		req(input$QUAL_filter, input$picked_chr)
		# Assuming 1 variant per position need to have handling of multiple vars in one place!!!!!
		# if(is.null(chrplot_click()) | length(chrplot_click()) == 0) {
		# 	df <- locusGenoTypes() %>%
		# 		filtfun()(samples())
		# } else {
		# 	df <- locusGenoTypes() %>%
		# 		slice(chrplot_click()$pointNumber + 1)
		# }
		
		df <- locusGenoTypes() %>%
			filtfun()(samples()) %>%
			quality_filters(input) %>%
			dplyr::filter(CHROM %in% input$picked_chr)
		
		# if(!is.null(input$chrplotclick)) {
		# 	df <- nearPoints(df, input$chrplotclick, "POS", "AF")
		# }
		
		#dplyr::pull(pos)
	}) %>%
	bindCache(
		input$DP_filter, input$QUAL_filter, input$QR_filter,
		input$QA_filter, input$AF_filter, 
		input$picked_chr
	) %>%
		bindEvent(input$go)
	
	# Genotype filtering ------------------------------------------------------
	## | Set sample names -----------------------------------------------------
	output$renameUI <- renderUI({
		req(sample_vars_tolisten())
		alias_samples(samples())
	})
	
	## | Genotype filters UI --------------------------------------------------
	output$setSelector <- renderUI({
		genotype_selector(samples(), input)
	})
	outputOptions(output, "setSelector", suspendWhenHidden = FALSE)
	
	## | Quality filters UI ---------------------------------------------------
	output$quality_sliders <- renderUI({
		quality_sliders(allLoci())
	})
	outputOptions(output, "quality_sliders", suspendWhenHidden = FALSE)
	
	# | Download data ----------------------------------------------------------
	output$downloadData <- downloadHandler(
		filename = "selected_mutants.tsv",
		content = function(file) {
			vroom::vroom_write(loci(), file)
		}
	)
	
	# Visual Outputs ----------------------------------------------------------
	## | Counts ---------------------------------------------------------------
	output$total_vb <- renderValueBox({
		req(input$vcf)
		valueBox(
			format(nrow(allLoci()), big.mark = ","),
			"Total Variants", icon = icon("dna"),
			color = "purple"
		)
	})
	
	output$nfiltered_vb <- renderValueBox({
		req(input$vcf)
		valueBox(
			format(nrow(loci()), big.mark = ","),
			"Filtered Variants", icon = icon("dna")
		)
	})
	
	
	## | Allele Frequency plot ------------------------------------------------
	#output$testpoints <- renderPrint(print(input$chrplotclick))
	#output$chrplot <- renderPlot({
	output$chrplot <- plotly::renderPlotly({
		#output$chrplot <- renderGirafe({
		plot <- loci() %>%
			loci_plot()
		#plot
		plotly::ggplotly(dynamicTicks = TRUE, plot) %>%
		layout_ggplotly() %>%
		plotly::layout(
			legend = list(
				#y = -0.1
				title = list(text = ""),
				valign = "bottom"#,
				#yanchor = "middle"
			)
		) %>%
		plotly::config(
			displaylogo = FALSE,
			modeBarButtonsToRemove = list("hoverCompareCartesian")
		) 
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
	
	# observeEvent(input$reset, {
	# 	js$resetClick()
	# })
	
	#output$chrplot_click <- renderText(names(chrplot_click()))
	
	## | Mutation type plot ---------------------------------------------------
	
	output$mutTypeFreqPlot <- plotly::renderPlotly({
		mut_type_freq_plot(loci())
	})
	
	## | Variant effect annotation --------------------------------------------
	output$effect <- DT::renderDataTable({
		req(input$filtVarsDT_row_last_clicked)
		gen_var_eff_DT(loci(), input$filtVarsDT_row_last_clicked)
	})
	
	genomeBrowserString <- reactiveValues()
	
	observeEvent(input$filtVarsDT_row_last_clicked, {
		req(input$filtVarsDT_row_last_clicked)
		row <- loci() %>% dplyr::slice(input$filtVarsDT_row_last_clicked)
		genomeBrowserString$string <- 
		wormbase_view(
			gsub("chr", "", row$CHROM),
			row$POS - 5000,
			row$POS + 5000,
			row$POS - 1,
			row$POS
		)
	})
	
	observeEvent(input$filteredDels_row_last_clicked, {
		req(input$filteredDels_row_last_clicked)
		row <- delTab() %>% dplyr::slice(input$filteredDels_row_last_clicked)
		genomeBrowserString$string <- 
			wormbase_view(
				gsub("chr", "", row$sequence),
				row$start - 5000,
				row$end + 5000,
				row$start,
				row$end
			)
	})
	
	output$genomeBrowser <- renderText({
		if(
			all(
				is.null(input$filtVarsDT_row_last_clicked),
				is.null(input$filteredDels_row_last_clicked)
			)
		) {
			return("Please select a variant to view...")
		}
		genomeBrowserString$string
		
	})
	
	## | __Main Variants Table__ ----------------------------------------------
	
	output$col_picker <- renderUI({
		variant_column_selector(vcftidy()$meta)
	})
	
	output$filtVarsDT <- DT::renderDataTable({
		gen_var_tab(loci(), samples(), aliases())
	}, server = TRUE)
	
	# Deletions ---------------------------------------------------------------
	
	gff <- reactive({
		req(input$gff)
		import_deletions(input$gff$datapath)
	})
	
	gff_wide <- reactive({
		gff() %>%
			dplyr::group_by(sequence, feature, start, end) %>%
			tidyr::pivot_wider(
				names_from = sample, values_from = c("score", "p_value")
			) %>%
			ungroup()
	})
	
	delTab <- reactive({
		gff_wide() %>% 
			delfiltfun()(samples()) %>%
			dplyr::filter(sequence %in% input$picked_chr)
	})
	
	output$filteredDels <- DT::renderDataTable({
		delTab() %>% 
			#DT::datatable(rownames = FALSE, selection = "single")
			gen_del_tab(samples(), aliases())
	}, server = TRUE)
	
	# | Download data ----------------------------------------------------------
	output$downloadDelData <- downloadHandler(
		filename = "selected_deletions.tsv",
		content = function(file) {
			vroom::vroom_write(delTab(), file)
		}
	)
}
