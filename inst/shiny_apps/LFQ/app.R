library(tidyverse)
library(magrittr)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(Rlibstree)
library(DT)
library(Biobase)
library(MSnbase)
library(limma)
library(vsn)
library(fdrtool)
library(proteomeR)
library(shiny)
library(shinydashboard)

ui <- shinyUI(
	dashboardPage(
		dashboardHeader(title = "proteomeR - LFQ"),
		dashboardSidebar(
		  sidebarMenu(
		    menuItem("Files", selected = TRUE,
		      fileInput('file1', 'ProteinGroups.txt',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
  		    fileInput('file2', 'ExperimentalDesign.txt',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
  		    radioButtons("anno", "Sample annotation", choices = list("Parse from columns" = "columns", "Use Experimental Design" = "expdesign"), selected = "columns")
  		  ),
		    menuItemOutput("columns"),
  			menuItem("Imputation & Stats options",
  			  radioButtons("imputation", "Imputation type", choices = list("QRILC" = "QRILC", "Minimal probability" = "MinProb", "Manual" = "man", "k-nearest neighbors" = "knn"), selected = "MinProb"),
  			  radioButtons("stat", "Statistical method", c("Linear model", "ANOVA"), selected = "Linear model")
  			),
  			actionButton("analyze", "Analyze"),
  			tags$hr(),
        uiOutput("downloadTable"),
  			uiOutput("downloadButton")
  		  )
			),
		dashboardBody(
		  helpText("Please cite: "),
		  fluidRow(
				box(numericInput("p", "P value cut off", min = 0.0001, max = 0.1, value = 0.05), width = 2),
				box(numericInput("lfc", "Fold change cut off (log2)", min = 0, max = 10, value = 1), width = 2),
				infoBoxOutput("signBox"),
				box(radioButtons("pres", "Data presentation", c("contrast", "centered"), selected = "contrast"), width = 2),
				box(radioButtons("contrasts", "Contrasts", c("control", "all"), selected = "control"), width = 2)
			),
			fluidRow(
			  column(width = 7,
			    box(title = "Top Table", box(uiOutput("select"), width = 6), box(uiOutput("exclude"), width = 6), DT::dataTableOutput("table"), width = 12)
			  ),
			  column(width = 5,
  				tabBox(title = "Result Plots", width = 12,
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
  				  )
				  ),
  				tabBox(title = "QC Plots", width = 12,
  				  tabPanel(title = "Protein Numbers",
  				    plotOutput("numbers", height = 600),
  				    downloadButton('downloadNumbers', 'Download')
  				  ),
  				  tabPanel(title = "Sample coverage",
  				    plotOutput("coverage", height = 600),
  				    downloadButton('downloadCoverage', 'Download')
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
)

server <- shinyServer(function(input, output) {
  options(shiny.maxRequestSize=60*1024^2)

  ### UI functions ### ----------------------------------------------------------------------------------------------------------------
  output$columns <- renderMenu({
    menuItem("Columns",
             selectizeInput("name", "Name column", choices=colnames(data()), selected = "Gene.names"),
             selectizeInput("id", "ID column", choices=colnames(data()), selected = "Protein.IDs"),
             selectizeInput("filt", "Filter on columns" , colnames(data()), multiple = TRUE, selected = c("Reverse","Potential.contaminant")),
             if (input$anno == "columns" & !is.null(data())) {
               cols <- grep("^LFQ", colnames(data()))
               prefix <- getCommonPrefix(data()[,cols] %>% colnames())
               selectizeInput("control", "Control", choices=make.names(colnames(data())[cols] %>% gsub(prefix,"",.) %>% substr(., 1, nchar(.)-1)))
             },
             if (input$anno == "expdesign" & !is.null(expdesign())) { selectizeInput("control", "Control", choices=make.names(expdesign()$condition)) }
    )
  })

  ### Reactive functions ### ----------------------------------------------------------------------------------------------------------
  data <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header=T, sep="\t", stringsAsFactors = F) %>% mutate(id = row_number())
  })
  expdesign <- reactive({
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header=T, sep="\t", stringsAsFactors = F) %>% mutate(id = row_number())
  })

  filt <- reactive({
    data <- data()
    cols <- grep("^LFQ", colnames(data))
    cols_filt <- grep(paste("^", input$filt, "$", sep = "", collapse = "|"), colnames(data))

    if (!is.null(cols_filt)) {
      if (length(cols_filt) == 1) {
        data %<>% filter(.[,cols_filt] != "+")
      } else {
        data %<>% filter(!apply(.[,cols_filt] == "+", 1, any))
      }
    }
    data %<>% unique_names(., input$name, input$id)

    if (input$anno == "columns") {
      data %<>% make_se_parse(., cols)
    }
    if (input$anno == "expdesign") {
      data %<>% make_se(., cols, expdesign())
    }
    data %>% filter_missval(., 0)
  })

  norm <- reactive({
    data <- filt()
    norm_vsn(data)
  })

  imp <- reactive({
    norm <- norm()
    imputation(norm, input$imputation)
  })

  df <- reactive({
    imp <- imp()
    if(input$stat == "ANOVA") {
      withProgress(message = 'Analysis', value = 0.1, {
        df <- anova_tukey(imp, input$control, input$contrasts)
      })
    }
    if(input$stat == "Linear model") {
      df <- linear_model(imp, input$control, input$contrasts)
    }
    df
  })

  sign <- reactive({
    df <- df()
    cutoffs(df, input$p, input$lfc)
  })

  ### All object and functions upon 'Analyze' input  ### ---------------------------------------------------------------------------

  observeEvent(input$analyze, {

    ### Interactive UI functions ### ----------------------------------------------------------------------------------------------
    output$downloadTable <- renderUI({
      selectizeInput("dataset", "Choose a dataset to download" , c("results","significant_proteins","displayed_subset","full_dataset"))
    })

    output$downloadButton <- renderUI({
      downloadButton('downloadData', 'Download')
    })

    output$signBox <- renderInfoBox({
      infoBox("Significant proteins", paste(sign() %>% .[rowData(.)$sign == "+"] %>% nrow(), " out of", sign() %>% nrow(), sep = " "), icon = icon("thumbs-up", lib = "glyphicon"), color = "green", width = 4)
    })

    output$select <- renderUI({
      row_data <- rowData(sign())
      cols <- grep("_sign",colnames(row_data))
      names <- colnames(row_data)[cols]
      names %<>% gsub("_sign","",.)
      selectizeInput("select", "Select direct comparisons", choices=names, multiple = TRUE)
    })

    output$exclude <- renderUI({
      row_data <- rowData(sign())
      cols <- grep("_sign",colnames(row_data))
      names <- colnames(row_data)[cols]
      names %<>% gsub("_sign","",.)
      selectizeInput("exclude", "Exclude direct comparisons", choices=names, multiple = TRUE)
    })

    output$volcano_cntrst <- renderUI({
      if (!is.null(selected())) {
        df <- rowData(selected())
        cols <- grep("_sign$",colnames(df))
        selectizeInput("volcano_cntrst", "Contrast", choices = gsub("_sign", "", colnames(df)[cols]))
      }
    })

    ### Reactive functions ### ----------------------------------------------------------------------------------------------
    excluded <- reactive({
      if(is.null(input$exclude)) {
        excluded <- sign()
      } else {
        if(length(input$exclude) == 1) {
          df <- rowData(sign())
          col <- grep(paste(input$exclude, "_sign", sep = ""), colnames(df))
          excluded <- sign()[df[,col] != "+",]
        } else {
          df <- rowData(sign())
          cols <- grep(paste(input$exclude, "_sign", sep = "", collapse = "|"), colnames(df))
          excluded <- sign()[apply(df[,cols] != "+", 1, all)]
        }
      }
      excluded
    })

    selected <- reactive({
      if(is.null(input$select)) {
        selected <- excluded()
      } else {
        if(length(input$select) == 1) {
          df <- rowData(excluded())
          col <- grep(paste(input$select, "_sign", sep = ""), colnames(df))
          selected <- sign()[df[,col] == "+",]
        } else {
          df <- rowData(excluded())
          cols <- grep(paste(input$select, "_sign", sep = "", collapse = "|"), colnames(df))
          selected <- sign()[apply(df[,cols] == "+", 1, all)]
        }
      }
      selected
    })

    res <- reactive({
      results(selected())
    })

    table <- reactive({
      res <- res() %>% filter(sign == "+") %>% select(-sign)
      if(input$pres == "centered") {
        cols <- grep("_ratio", colnames(res))
        table <- res[,-cols]
        colnames(table)[1:2] <- c("Protein Name", "Protein ID")
        colnames(table)[grep("sign", colnames(table))] %<>% gsub("[.]", " - ", .)
        colnames(table) %<>% gsub("_centered", "", .) %>% gsub("[_]", " ", .)
      }
      if(input$pres == "contrast") {
        cols <- grep("_centered", colnames(res))
        table <- res[,-cols]
        colnames(table)[1:2] <- c("Protein Name", "Protein ID")
        colnames(table)[grep("sign", colnames(table))] %<>% gsub("[.]", " - ", .)
        colnames(table) %<>% gsub("_ratio", "", .) %>% gsub("[_]", " ", .)
      }
      table
    })

    selected_plot_input <- reactive ({
      if(!is.null(input$table_rows_selected)) {
        selected_id <- table()[input$table_rows_selected,1]
        plot_single(selected(), selected_id, input$pres)
      }
    })

    heatmap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        plot_heatmap(selected(), input$pres, input$k, input$limit)
      })
    })

    volcano_input <- reactive({
      if(!is.null(input$volcano_cntrst)) {
        plot_volcano(selected(), input$volcano_cntrst, input$fontsize, input$check_names)
      }
    })

    norm_input <- reactive({
      plot_norm(filt(), norm())
    })

    missval_input <- reactive({
      plot_missval(norm())
    })

    numbers_input <- reactive({
      plot_numbers(norm())
    })

    coverage_input <- reactive({
      plot_coverage(norm())
    })

    ### Output functions ### ----------------------------------------------------------------------------------------------
    output$table <- DT::renderDataTable({
      table()
    }, options = list(pageLength = 25, scrollX = T), selection = list(mode = 'single', selected = c(1)))

    output$selected_plot <- renderPlot({
      selected_plot_input()
    })

    output$heatmap <- renderPlot({
      heatmap_input()
    })

    output$volcano <- renderPlot({
      volcano_input()
    })

    output$norm <- renderPlot({
      norm_input()
    })

    output$missval <- renderPlot({
      missval_input()
    })

    output$numbers <- renderPlot({
      numbers_input()
    })

    output$coverage <- renderPlot({
      coverage_input()
    })

    observe({
      output$plot <- renderUI({
        plotOutput("heatmap", height = (100 * as.numeric(input$size)))
      })
    })

    ### Download objects and functions ### ---------------------------------------------------------------------------------
    datasetInput <- reactive({
      switch(input$dataset,
             "results" = results(sign()),
             "significant_proteins" = results(sign()) %>% filter(sign == "+") %>% select(-sign),
             "displayed_subset" = res() %>% filter(sign == "+") %>% select(-sign),
             "full_dataset" = left_join(rownames_to_column(exprs(sign()) %>% data.frame()), fData(sign()), by = c("rowname" = "name")))
    })

    output$downloadData <- downloadHandler(
      filename = function() { paste(input$dataset, ".txt", sep = "") },
      content = function(file) { write.table(datasetInput(), file, col.names = T, row.names = F, sep ="\t") }
    )

    output$downloadPlot <- downloadHandler(
      filename = function() { paste("Barplot_", table()[input$table_rows_selected,1], ".pdf", sep = "") },
      content = function(file) {
        pdf(file)
        print(selected_plot_input())
        dev.off()
      }
    )

    output$downloadHeatmap <- downloadHandler(
      filename = 'Heatmap.pdf',
      content = function(file) {
        pdf(file, height = 10, paper = "a4")
        print(heatmap_input())
        dev.off()
      }
    )

    output$downloadVolcano <- downloadHandler(
      filename = function() { paste("Volcano_", input$contrast, ".pdf", sep = "") },
      content = function(file) {
        pdf(file)
        print(volcano_input())
        dev.off()
      }
    )

    output$downloadNorm <- downloadHandler(
      filename = "normalization.pdf",
      content = function(file) {
        pdf(file)
        print(norm_input())
        dev.off()
      }
    )

    output$downloadMissval <- downloadHandler(
      filename = "missing_values_heatmap.pdf",
      content = function(file) {
        pdf(file)
        print(missval_input())
        dev.off()
      }
    )

    output$downloadNumbers <- downloadHandler(
      filename = "numbers.pdf",
      content = function(file) {
        pdf(file)
        print(numbers_input())
        dev.off()
      }
    )

    output$downloadCoverage <- downloadHandler(
      filename = "coverage.pdf",
      content = function(file) {
        pdf(file)
        print(coverage_input())
        dev.off()
      }
    )
  })
})


shinyApp(ui = ui, server = server)
