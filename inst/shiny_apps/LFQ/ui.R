shinyUI(
	dashboardPage(
		dashboardHeader(title = "proteomeR - LFQ"),
		dashboardSidebar(
		  sidebarMenu(
		    fileInput('file1', 'ProteinGroups.txt',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
  		  menuItem("Experimental Design",
  		    fileInput('file2', 'ExperimentalDesign.txt',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
  		    radioButtons("anno", "Sample annotation", choices = list("Parse from columns" = "columns", "Use Experimental Design" = "expdesign"), selected = "columns")
  		  ),
		    menuItemOutput("columns"),
  			menuItem("Imputation & Stats options",
  			  radioButtons("imputation", "Imputation type", c("QRILC method", "Minimal probability", "Perseus", "k-nearest neighbors"), selected = "QRILC method"),
  			  radioButtons("stat", "Statistical method", c("Linear model", "ANOVA"), selected = "Linear model")
  			),
  			actionButton("analyze", "Analyze"),
  			tags$hr(),
        uiOutput("downloadTable"),
  			uiOutput("downloadButton")
  		  )
			),
		dashboardBody(
		  helpText("Please cite: Smits et al (2017)"),
		  fluidRow(
				box(numericInput("p", "P value cut off", min = 0.0001, max = 0.1, value = 0.05), width = 2),
				box(numericInput("lfc", "Fold change cut off (log2)", min = 0, max = 10, value = 1), width = 2),
				infoBoxOutput("signBox"),
				box(radioButtons("pres", "Data presentation", c("contrast", "centered"), selected = "contrast"), width = 2),
				box(radioButtons("contrasts", "Contrasts", c("control", "all"), selected = "control"), width = 2)
			),
			fluidRow(
			  box(title = "Significant proteins", box(uiOutput("select"), width = 6), box(uiOutput("exclude"), width = 6), DT::dataTableOutput("table"), width = 7),
				tabBox(title = "Plots", width = 5,
				  tabPanel(title = "Selected Protein", plotOutput("selected_plot"), downloadButton('downloadPlot', 'Download plot')),
				  tabPanel(title = "Heatmap",
				    fluidRow(
				      box(numericInput("k", "Kmeans clusters", min = 0, max = 15, value = 7), width = 4),
				      box(numericInput("limit", "Color limit (log2)", min = 0, max = 16, value = 6), width = 4),
				      box(numericInput("size", "Heatmap size (4-30)", min = 4, max = 30, value = 10), width = 4)
				    ),
				    fluidRow(
				      uiOutput("plot"),
				      downloadButton('downloadHeatmap', 'Download heatmap'))
				  ),
				  tabPanel(title = "Volcano plot",
				    fluidRow(
				      box(uiOutput("volcano_cntrst"), width = 6),
				      box(numericInput("fontsize", "Font size", min = 0, max = 8, value = 4), width = 3),
				      box(checkboxInput("check_names", "Display names", value = TRUE), width = 3)
				    ),
				    fluidRow(
				      plotOutput("volcano", height = 600),
				      downloadButton('downloadVolcano', 'Download volcano')
				    )
				  ),
				  tabPanel(title = "Normalization",
               plotOutput("norm", height = 600),
	             downloadButton('downloadNorm', 'Download')
	        ),
				  tabPanel(title = "Missing values",
				           plotOutput("missval", height = 600),
				           downloadButton('downloadMissval', 'Download')
				  )
				)
			)
		)
	)
)
