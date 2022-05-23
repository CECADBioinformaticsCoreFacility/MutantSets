# UI --------------------------------------------------------------------------
## |Sidebar -------------------------------------------------------------------

# Define UI for application that draws a histogram
sidebar <- bs4Dash::dashboardSidebar(
	bs4Dash::sidebarMenu(
		fileInput("vcf", "Select a VCF file", accept = ".vcf"),
		fileInput("gff", "Select a gff file", accept = ".gff"),
		actionButton("go","Start / Apply Filters"),
		#menuItem("Options", tabName = "options", icon = icon("th")),
		menuItem(
			"Filtering", tabName = "table", icon = icon("table")
		)#,
		# menuItem(
		# 	"Overview", tabName = "overview", icon = icon("th")
		# )#,
	)
)
## |Body ----------------------------------------------------------------------
body <- bs4Dash::dashboardBody(... = 
	shinybusy::add_busy_spinner(spin = "fading-circle"),
	bs4Dash::tabItems(
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
		bs4Dash::tabItem(
			tabName = "table",
			fluidRow(
				# box(
				# 	title = "Genotype Filters", 
				# 	status = "primary", solidHeader = TRUE, collapsible = TRUE,
				# 	width = 4,
				# 	uiOutput("setSelector")
				# ),
				bs4Dash::tabBox(
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
						#verbatimTextOutput("test"),
						"*Note that at loci with multiple variants 1/1 will also allow homozygous alternate alleles e.g. 2/2 etc."
					),
					tabPanel(
						title = "Quality Filters",
						uiOutput("quality_sliders")
					)
				),
				## || Plots tabs ----
				bs4Dash::tabBox(
					title = "Plots", width = 8,
					sidebar = bs4Dash::boxSidebar(
						id = "allele_plot_sidebar",
						width = 25,
						numericInput("width", "width", 9, 0, 96000),
						numericInput("height", "height", 4, 0, 96000),
						selectInput(
							"type", "file type",
							selected = "png", choices = c(
								"png", "svg", "pdf", "jpeg", "eps", "ps",
								"tiff", "bmp"
							)
						),
						selectInput(
							"units", "units", selected = "in",
							choices = c("in", "cm", "mm", "px")
						),
						numericInput(
							"dpi", "dpi",
							value =  300, min = 1, max = 96000, step = 1
						),
						downloadButton("save_image", "Save Image")
						#textInput()
					),
					tabPanel(
						"Allele Frequency",
						shinyjs::useShinyjs(),
						# shinyjs::extendShinyjs(
						# 	text = paste0(
						# 		"shinyjs.resetClick = function() {",
						# 		"Shiny.onInputChange(",
						# 		"'.clientValue-plotly_click-A', 'null'",
						# 		"); ",
						# 		"}"
						# 	)
						# ),
						#plotOutput("chrplot", click = "chrplotclick")#,
						plotly::plotlyOutput("chrplot")#,
						#verbatimTextOutput("chrplot_sel")
						# actionButton("reset","Reset Point selection")#,
						#verbatimTextOutput("chrplot_click")
						#verbatimTextOutput("tmp")
						#girafeOutput("chrplot")
						#verbatimTextOutput("testpoints")
					),
					tabPanel(
						"Variant Density",
						plotly::plotlyOutput("vdplot")
					),
					tabPanel(
						"Effect",
						DT::DTOutput("effect")
					),
					tabPanel(
						"Genome Browser",
						htmlOutput("genomeBrowser")
					),
					# tabPanel(
					# 	"Sets"
					# ),
					tabPanel(
						"Variant Types",
						fluidRow(
							#column(width = 1),
							column(width = 8,
								plotly::plotlyOutput("mutTypeFreqPlot")
							),
							column(width = 3,
								fluidRow(
									valueBoxOutput(
										"total_vb", width = "80%"
									)
								),
								fluidRow(
									valueBoxOutput(
										"nfiltered_vb", width = "80%"
									)
								)
							)
						)
					)
				)
			),
			fluidRow(
				bs4Dash::tabBox(
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
						DT::DTOutput("filteredDels"),
						downloadButton("downloadDelData")
					)
				)
			)
		)#,
		### ||Overview
		# tabItem(
		# 	tabName = "overview", # quality filters?
		# 	fluidRow(
		# 		box(
		# 			title = "Quality Ranges", 
		# 			status = "primary", solidHeader = TRUE, collapsible = TRUE,
		# 			width = 4,
		# 			#plotOutput("qualplot")
		# 		),
		# 		box(
		# 			title = "Chr Plots", 
		# 			status = "primary", solidHeader = TRUE, collapsible = TRUE,
		# 			width = 8,
		# 			#uiOutput("chr"),
		# 			#plotOutput("chrPlot")
		# 		)
		# 	)
		# )
	)
)

## |UI Wrapper ----------------------------------------------------------------
ui <- bs4Dash::dashboardPage(
	bs4Dash::dashboardHeader(title = "MutantSets"),
	sidebar,
	body
)
