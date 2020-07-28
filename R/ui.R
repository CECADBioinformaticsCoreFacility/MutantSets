# UI --------------------------------------------------------------------------
## |Sidebar -------------------------------------------------------------------

# Define UI for application that draws a histogram
sidebar <- shinydashboard::dashboardSidebar(
	shinydashboard::sidebarMenu(
		#menuItem("Inputs", tabName = "inputs", icon = icon("upload")),
		fileInput("vcf", "Select a VCF file", accept = ".vcf"),
		fileInput("gff", "Select a gff file", accept = ".gff"),
		#menuItem("Options", tabName = "options", icon = icon("th")),
		shinydashboard::menuItem(
			"Filtering", tabName = "table", icon = icon("table")
		),
		shinydashboard::menuItem(
			"Overview", tabName = "overview", icon = icon("th")
		)#,
		#uiOutput("DPslider"),
		#uiOutput("MQslider"),
		#uiOutput("min_QUAL")
	)
)
## |Body ----------------------------------------------------------------------
body <- shinydashboard::dashboardBody(
	shinybusy::add_busy_spinner(spin = "fading-circle"),
	shinydashboard::tabItems(
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
		shinydashboard::tabItem(
			tabName = "table",
			fluidRow(
				# box(
				# 	title = "Genotype Filters", 
				# 	status = "primary", solidHeader = TRUE, collapsible = TRUE,
				# 	width = 4,
				# 	uiOutput("setSelector")
				# ),
				shinydashboard::tabBox(
					title = "Filters", width = 4,
					tabPanel(
						title = "Sample Aliases",
						uiOutput("renameUI")
					),
					tabPanel(
						title = "Genotype Filters", status = "primary",
						solidHeader = TRUE, collapsible = TRUE,
						uiOutput("chr_selection"),
						uiOutput("setSelector"),
						"*Note that at loci with multiple variants 1/1 will also allow homozygous alternate alleles e.g. 2/2 etc."
					),
					tabPanel(
						title = "Quality Filters",
						uiOutput("quality_sliders")
					)
				),
				shinydashboard::tabBox(
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
									shinydashboard::valueBoxOutput(
										"total_vb", width = "80%"
									)
								),
								fluidRow(
									shinydashboard::valueBoxOutput(
										"nfiltered_vb", width = "80%"
									)
								)
							)
						)
					)
				)
			),
			fluidRow(
				shinydashboard::tabBox(
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
		shinydashboard::tabItem(
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
ui <- shinydashboard::dashboardPage(
	shinydashboard::dashboardHeader(title = "MutantSets"),
	sidebar,
	body
)
