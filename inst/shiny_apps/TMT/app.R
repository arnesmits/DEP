library(dplyr)
library(SummarizedExperiment)
library(DEP)
library(shiny)
library(shinydashboard)

ui <- shinyUI(
  dashboardPage(
    dashboardHeader(title = "DEP - TMT"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Files", selected = TRUE,
                 fileInput('file1', 'ProteinGroups.txt',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                 fileInput('file2', 'ExperimentalDesign.txt',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
        ),
        menuItemOutput("columns"),
        menuItem("Imputation options",
                 radioButtons("imputation", "Imputation type", choices = c("man", MSnbase::imputeMethods())[1:9], selected = "MinProb"),
                 p(a("Detailed information link", href = "https://www.rdocumentation.org/packages/MSnbase/versions/1.20.7/topics/impute-methods", target="_blank"))
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
        box(numericInput("p", "adj. P value", min = 0.0001, max = 0.1, value = 0.05), width = 2),
        box(numericInput("lfc", "Log2 fold change", min = 0, max = 10, value = 1), width = 2),
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
                      tabPanel(title = "Selected Protein", plotOutput("selected_plot"), downloadButton('downloadPlot', 'Save plot')),
                      tabPanel(title = "Heatmap",
                               fluidRow(
                                 box(numericInput("k", "Kmeans clusters", min = 0, max = 15, value = 7), width = 4),
                                 box(numericInput("limit", "Color limit (log2)", min = 0, max = 16, value = 6), width = 4),
                                 box(numericInput("size", "Heatmap size (4-30)", min = 4, max = 30, value = 10), width = 4)
                               ),
                               fluidRow(
                                 uiOutput("plot"),
                                 downloadButton('downloadHeatmap', 'Save heatmap'))
                      ),
                      tabPanel(title = "Volcano plot",
                               fluidRow(
                                 box(uiOutput("volcano_cntrst"), width = 6),
                                 box(numericInput("fontsize", "Font size", min = 0, max = 8, value = 4), width = 3),
                                 box(checkboxInput("check_names", "Display names", value = TRUE), width = 3)
                               ),
                               fluidRow(
                                 plotOutput("volcano", height = 600),
                                 downloadButton('downloadVolcano', 'Save volcano')
                               )
                      )
               ),
               tabBox(title = "QC Plots", width = 12,
                      tabPanel(title = "Protein Numbers",
                               plotOutput("numbers", height = 600),
                               downloadButton('downloadNumbers', 'Save')
                      ),
                      tabPanel(title = "Sample coverage",
                               plotOutput("coverage", height = 600),
                               downloadButton('downloadCoverage', 'Save')
                      ),
                      tabPanel(title = "Normalization",
                               plotOutput("norm", height = 600),
                               downloadButton('downloadNorm', 'Save')
                      ),
                      tabPanel(title = "Missing values - Quant",
                               plotOutput("detect", height = 600),
                               downloadButton('downloadDetect', 'Save')
                      ),
                      tabPanel(title = "Missing values - Heatmap",
                               plotOutput("missval", height = 600),
                               downloadButton('downloadMissval', 'Save')
                      ),
                      tabPanel(title = "Imputation",
                               plotOutput("imputation", height = 600),
                               downloadButton('downloadImputation', 'Save')
                      )
               )
        )
      )
    )
  )
)

server <- shinyServer(function(input, output) {
  options(shiny.maxRequestSize=200*1024^2)

  ### UI functions ### ----------------------------------------------------------------------------------------------------------------
  output$columns <- renderMenu({
    menuItem("Columns",
             selectizeInput("name", "Name column", choices=colnames(data()), selected = "gene_name"),
             selectizeInput("id", "ID column", choices=colnames(data()), selected = "protein_id"),
             if (!is.null(expdesign())) { selectizeInput("control", "Control", choices=make.names(expdesign()$condition)) }
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
    rows <- grep("###", data[,grep(input$name, colnames(data))])
    cols <- grep("^signal_sum", colnames(data))
    data[-rows,] %>% make_unique(., input$name, input$id, delim = "[|]") %>% make_se(., cols, expdesign()) %>% filter_missval(., 0)
  })

  norm <- reactive({
    data <- filt()
    normalize_vsn(data)
  })

  imp <- reactive({
    norm <- norm()
    DEP::impute(norm, input$imputation)
  })

  df <- reactive({
    imp <- imp()
    test_diff(imp, input$control, input$contrasts)
  })

  sign <- reactive({
    df <- df()
    add_rejections(df, input$p, input$lfc)
  })

  ### All object and functions upon 'Analyze' input  ### ---------------------------------------------------------------------------

  observeEvent(input$analyze, {

    ### Interactive UI functions ### ----------------------------------------------------------------------------------------------
    output$downloadTable <- renderUI({
      selectizeInput("dataset", "Choose a dataset to save" , c("results","significant_proteins","displayed_subset","full_dataset"))
    })

    output$downloadButton <- renderUI({
      downloadButton('downloadData', 'Save')
    })

    output$signBox <- renderInfoBox({
      infoBox("Significant proteins", paste(sign() %>% .[rowData(.)$significant, ] %>% nrow(), " out of", sign() %>% nrow(), sep = " "), icon = icon("thumbs-up", lib = "glyphicon"), color = "green", width = 4)
    })

    output$select <- renderUI({
      row_data <- rowData(sign())
      cols <- grep("_significant",colnames(row_data))
      names <- colnames(row_data)[cols]
      names <- gsub("_significant", "", names)
      selectizeInput("select", "Select direct comparisons", choices=names, multiple = TRUE)
    })

    output$exclude <- renderUI({
      row_data <- rowData(sign())
      cols <- grep("_significant",colnames(row_data))
      names <- colnames(row_data)[cols]
      names <- gsub("_significant", "", names)
      selectizeInput("exclude", "Exclude direct comparisons", choices=names, multiple = TRUE)
    })

    output$volcano_cntrst <- renderUI({
      if (!is.null(selected())) {
        df <- rowData(selected())
        cols <- grep("_significant$",colnames(df))
        selectizeInput("volcano_cntrst", "Contrast", choices = gsub("_significant", "", colnames(df)[cols]))
      }
    })

    ### Reactive functions ### ----------------------------------------------------------------------------------------------
    excluded <- reactive({
      if(is.null(input$exclude)) {
        excluded <- sign()
      } else {
        if(length(input$exclude) == 1) {
          df <- rowData(sign())
          col <- grep(paste(input$exclude, "_significant", sep = ""), colnames(df))
          excluded <- sign()[!df[,col],]
        } else {
          df <- rowData(sign())
          cols <- grep(paste(input$exclude, "_significant", sep = "", collapse = "|"), colnames(df))
          excluded <- sign()[apply(!df[,cols], 1, all)]
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
          col <- grep(paste(input$select, "_significant", sep = ""), colnames(df))
          selected <- excluded()[df[,col],]
        } else {
          df <- rowData(excluded())
          cols <- grep(paste(input$select, "_significant", sep = "", collapse = "|"), colnames(df))
          selected <- excluded()[apply(df[,cols], 1, all)]
        }
      }
      selected
    })

    res <- reactive({
      get_results(selected())
    })

    table <- reactive({
      res <- res() %>% filter(significant) %>% select(-significant)
      if(input$pres == "centered") {
        cols <- grep("_ratio", colnames(res))
        table <- res[,-cols]
        colnames(table)[1:2] <- c("Protein Name", "Protein ID")
        colnames(table)[grep("significant", colnames(table))] <-
          gsub("[.]", " - ", colnames(table)[grep("significant", colnames(table))])
        colnames(table) <- gsub("_centered", "", colnames(table)) %>% gsub("[_]", " ", .)
      }
      if(input$pres == "contrast") {
        cols <- grep("_centered", colnames(res))
        table <- res[,-cols]
        colnames(table)[1:2] <- c("Protein Name", "Protein ID")
        colnames(table)[grep("significant", colnames(table))] <-
          gsub("[.]", " - ", colnames(table)[grep("significant", colnames(table))])
        colnames(table) <- gsub("_ratio", "", colnames(table)) %>% gsub("[_]", " ", .)
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
      plot_normalization(filt(), norm())
    })

    missval_input <- reactive({
      plot_missval(norm())
    })

    detect_input <- reactive({
      plot_detect(norm())
    })

    imputation_input <- reactive({
      plot_imputation(norm(), df())
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

    output$detect <- renderPlot({
      detect_input()
    })

    output$imputation <- renderPlot({
      imputation_input()
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
             "results" = get_results(sign()),
             "significant_proteins" = get_results(sign()) %>% filter(significant) %>% select(-significant),
             "displayed_subset" = res() %>% filter(significant) %>% select(-significant),
             "full_dataset" = left_join(rownames_to_column(assay(sign()) %>% data.frame()), data.frame(rowData(sign())), by = c("rowname" = "name")))
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

    output$downloadDetect <- downloadHandler(
      filename = "missing_values_quant.pdf",
      content = function(file) {
        pdf(file)
        gridExtra::grid.arrange(detect_input())
        dev.off()
      }
    )

    output$downloadImputation <- downloadHandler(
      filename = "imputation.pdf",
      content = function(file) {
        pdf(file)
        print(imputation_input())
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
