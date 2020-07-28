# UI --------------------------------------------------------------------------
## |Sidebar -------------------------------------------------------------------

# Define UI for application that draws a histogram
sidebar <- dashboardSidebar(
	sidebarMenu(
		#menuItem("Inputs", tabName = "inputs", icon = icon("upload")),
		fileInput("vcf", "Select a VCF file", accept = ".vcf"),
		fileInput("gff", "Select a gff file", accept = ".gff"),
		#menuItem("Options", tabName = "options", icon = icon("th")),
		menuItem("Filtering", tabName = "table", icon = icon("table")),
		menuItem("Overview", tabName = "overview", icon = icon("th"))#,
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
		# tabItem(
		# 	tabName = "options",
		# 	fluidRow(
		# 		box(
		# 			title = "Sample Aliases",
		# 			uiOutput("renameUI")
		# 		)
		# 	)
		# ),
		### ||Table - ---------------------------------------------------------
		tabItem(
			tabName = "table",
			fluidRow(
				# box(
				# 	title = "Genotype Filters", 
				# 	status = "primary", solidHeader = TRUE, collapsible = TRUE,
				# 	width = 4,
				# 	uiOutput("setSelector")
				# ),
				tabBox(
					title = "Filters", width = 4,
					tabPanel(
						title = "Sample Aliases",
						uiOutput("renameUI")
					),
					tabPanel(
						title = "Genotype Filters", status = "primary",
						solidHeader = TRUE, collapsible = TRUE,
						uiOutput("chr_selection"),
						uiOutput("setSelector")
					),
					tabPanel(
						title = "Quality Filters",
						uiOutput("quality_sliders")
					)
				),
				tabBox(
					title = "Plots", width = 8,
					tabPanel(
						"Allele Frequency",
						shinyjs::useShinyjs(),
						shinyjs::extendShinyjs(
							text = paste0(
								"shinyjs.resetClick = function() {",
								"Shiny.onInputChange(",
								"'.clientValue-plotly_click-A', 'null'",
								"); ",
								"}"
							)
						),
						#plotOutput("chrplot", click = "chrplotclick")#,
						plotly::plotlyOutput("chrplot"),
						#verbatimTextOutput("chrplot_sel")
						actionButton("reset","Reset Point selection")#,
						#verbatimTextOutput("chrplot_click")
						#verbatimTextOutput("tmp")
						#girafeOutput("chrplot")
						#verbatimTextOutput("testpoints")
					),
					tabPanel(
						"Effect",
						DT::DTOutput("effect")
					),
					tabPanel(
						"Sets"
					),
					tabPanel(
						"Variant Types",
						fluidRow(
							#column(width = 1),
							column(width = 8,
								plotly::plotlyOutput("mutTypeFreqPlot")
							),
							column(width = 3,
								fluidRow(
									valueBoxOutput("total_vb", width = "80%")
								),
								fluidRow(
									valueBoxOutput("nfiltered_vb", width = "80%")
								)
							)
						)
					)
				)
			),
			fluidRow(
				tabBox(
					title = "Variants", width = 12,
					tabPanel(
						"Variants",
						#verbatimTextOutput("test"),
						uiOutput("col_picker"),
						DT::DTOutput("filtVarsDT"),
						downloadButton("downloadData")
					),
					tabPanel(
						"MiModD deletions",
						# verbatimTextOutput("test"),
						DT::DTOutput("filteredDels")
					)
				)
				# box(
				# 	title = "Variants", 
				# 	status = "primary", solidHeader = TRUE, collapsible = TRUE,
				# 	width = 12,
				# 	#DT::DTOutput("tidytab")
				# 	#dropdownButton(
				# 		uiOutput("col_picker"),
				# 	#),
				# 	DT::DTOutput("filtVarsDT"),
				# 	downloadButton("downloadData")
				# )
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
		)
	)
)

## |UI Wrapper ----------------------------------------------------------------
ui <- dashboardPage(
	dashboardHeader(title = "MutantSets"),
	sidebar,
	body
)
