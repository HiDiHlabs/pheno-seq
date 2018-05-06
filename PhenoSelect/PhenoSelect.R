library(ggplot2)
library(shiny)
library(shinyFiles)
library(DT)
library(plyr)
library(reshape)


Sys.setenv(LANG = "en")
rm(list = ls())



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

ui <- fluidPage(sidebarLayout(
  sidebarPanel(
    title = 'Samples to display',
    conditionalPanel(
      'input.dataset === "Select wells"',
      uiOutput("sample_name_and_chip_number"),
      radioButtons(
        "cells_or_controls",
        label = "Spheroids or controls?",
        c('Spheroids' = "cells", 'Controls' = "controls")
      ),
      conditionalPanel(
        'input.cells_or_controls === "cells"',
        fluidRow(column(2,
                        uiOutput("choose_samples")),
                 column(
                   2,
                   radioButtons("single_or_all", label =
                                  "Wells displayed", c('All' = "ds.2", 'Single' = "ds.2s"))
                 )),
        fluidRow(
          column(
            10,
            plotOutput(
              "hist_bars",
              width = 400,
              dblclick = "hist_bars_dblclick",
              brush = brushOpts(id = "hist_bars_brush",
                                resetOnNew = FALSE)
            )
          ),
          column(
            5,
            radioButtons(
              "parameters_to_plot",
              label = "Parameters to plot",
              c('Size' = "Area", 'Circularity' = "Circularity")
            ),
            uiOutput("choose_Area"),
            uiOutput("choose_circularity"),
            sliderInput(
              "numberofbins",
              "Number of bins:",
              min = 1,
              max = 200,
              value = 100
            )
          )
        )
      )
    ),
    conditionalPanel(
      'input.dataset === "Selected for dispense"',
      fileInput(
        'uploaded_filter_file',
        "Choose filter file to upload",
        accept = c('text/csv', 'text/comma-seperated-values, text/plain', '.csv')
      ),
      actionButton("unselect_button_all", "Reset selection"),
      plotOutput("dispense_map", width = 400),
      h4(textOutput("number_selected")),
      br(),
      downloadButton("save_filterfile", "Save filter file"),
      br(),
      br()
    ),
    conditionalPanel(
      'input.dataset === "t-sne plot"',
      uiOutput("choose_RData"),
      uiOutput("load_RData"),
      br(),
      div(uiOutput("tsne_features"), style = "font-size: 75%; width: 75%")
    )
  ),
  mainPanel(
    tabsetPanel(
      id = 'dataset',
      tabPanel(
        'Select experiment',
        fluidRow(
          h1("PhenoSelect"),
          h5("Please select input data first."),
          br(),
          shinyDirButton("dir", "Choose directory of your experiment", "Set new location"),
          br(),
          br(),
          h4("Choosen home path"),
          textOutput("path"),
          br(),
          h4("Please choose the .csv file containing image features"),
          uiOutput("choose_csv"),
          br(),
          h4("Please enter the chip number")
        ),
        fluidRow(
          column(
            6,
            numericInput("num", label = "Chip number", value =
                           NA),
            textInput("A1", label = "Sample in well A1", value = NA),
            textInput("B1", label = "Sample in well B1", value = NA),
            textInput("C1", label = "Sample in well C1", value = NA),
            textInput("D1", label = "Sample in well D1", value = NA),
            uiOutput("load_dataset")
          ),
          column(
            6,
            selectInput(
              "cell_type",
              label = "Select cell type",
              c("HD1495", "MCF10CA"),
              multiple = FALSE
            ),
            textInput("A2", label = "Sample in well A2", value = NA),
            textInput("B2", label = "Sample in well B2", value = NA),
            textInput("C2", label = "Sample in well C1", value = NA),
            textInput("D2", label = "Sample in well D2", value = NA)
          )
        )
      ),
      tabPanel(
        'Select wells',
        fluidRow(column(6,
               h4(
                 textOutput('well_name')
               ))),
        br(),
        fluidRow(column(4,
                        div(DT::dataTableOutput('wells'),style="font-size: 75%; width: 75%"))),
        fluidRow(column(
          6,
          imageOutput("image_well_selection")
        ),
        column(6,
               actionButton("select_button_single", "Select for dispense"),
               br(),
               br(),
               radioButtons(
                 "channel_selection",
                 label = "Displayed channel",
                 c(
                   'Hoechst' = "Hoechst",
                   'CellTracker' = "CellTracker",
                   'Overlay' = "Merged"
                 )
               )
        ))
      ),
      tabPanel(
        'Selected for dispense',
        fluidRow(
          column(
            3,
            textInput("comment", h4("Add comments here"), value =
                        "Add a comment."),
            actionButton("add_comment_button", "Add")
          ),
          column(
            3,
            radioButtons(
              "cells_or_controls_selected",
              label = h4("Select a dataset"),
              c('Spheroids' = "cells", 'Controls' = "controls")
            )
          ),
          column(
            3,
            h4("Save selection"),
            downloadButton("save_selection", "Save selection as .csv")
          )
        ),
        fluidRow(
          br(),
          h5("Please click on the row number to view images and statistics."),
          br(),
          column(12,
                 div(DT::dataTableOutput('selected_table')),style="font-size: 75%; width: 75%")
          ),
        fluidRow(
          br(),
          column(
            5,
            h4(textOutput("well_name_selected")),
            imageOutput("image_well_selected")
          ),
          column(
              2,
              br(),
              actionButton("unselect_button_single", "Unselect well"),
              br(),
              br(),
              radioButtons(
                "channel_selected",
                label = "Displayed channel",
                c(
                  'Hoechst' = "Hoechst",
                  'CellTracker' = "CellTracker",
                  'Overlay' = "Merged"
                )
              ),
              radioButtons(
                "parameters_to_plot_selection",
                label = "Parameters to plot",
                c(
                  'Size' = "Area",
                  'Cell number' = "estimated_cell_number",
                  'Circularity' = "Circularity"
                )
              ),
              radioButtons(
                "dataset_to_plot_selection",
                label = "Dataset to plot",
                c(
                  'Selected wells' = "selected_wells",
                  'Single spheroid wells' = "single_cell_wells",
                  'All spheroids' = "all_cells"
                )
              )
          ),
          column(5,
              plotOutput("dotplot_highlighted", width =
                           200)
          )
        )
      ),
      tabPanel(
        't-sne plot',
        fluidRow(column(
          5,
          h3("T-SNE plot as imported from Pagoda."),
          h4("Mark region and double click to zoom into plot. Click on point to view corresponding image."),
          uiOutput("choose_parameter")
        ),
          column(
            1,
            br(),
            radioButtons(
              "tsne_channel_selected",
              label = "Displayed channel",
              c(
                'Hoechst' = "Hoechst",
                'CellTracker' = "CellTracker",
                'Overlay' = "Merged"
              ))
          ),
        column(
          1,
          br(),
          downloadButton("downloadPlot","Save current plot")
        )),
        fluidRow(
          column(
            7,
            plotOutput("tsne_plot", width =
                         500,
                       dblclick = "tsne_plot_dblclick",
                       brush = brushOpts(id = "tsne_plot_brush",
                                         resetOnNew = FALSE),
                       click = "tsne_plot_click")
          ),
          column(5,
            imageOutput("tsne_image")
                 )
        )
      )
    ) #close tabsetPanel
  ) #close mainPanel
) #close sidebarLayout
) #close FluidPage


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

server <- function(input, output) {
  # dir
  shinyDirChoose(input,
                 'dir',
                 roots = c(home = "~"))
  dir <- reactive(input$dir)
  output$dir <- renderPrint(dir())
  
  # path
  path <- reactive({
    home <- normalizePath("~")
    file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  })
  output$path <- renderPrint(path())
  
  # display .csv files in chosen directory
  output$choose_csv <- renderUI({
    tagList(
      selectInput(
        'csv_file',
        label = NULL,
        list.files(path(), pattern = ".csv"),
        multiple = FALSE,
        selectize = FALSE
      )
    )
  })
  
  
  # load dataset button, allows loading of data for the remaining tabs
  output$load_dataset <- renderUI({
    validate(
      need(input$dir != '', 'Please choose a file location!'),
      need(input$csv_file != '', 'Please choose a .csv file!'),
      need(input$cell_type != '', 'Please choose a cell type!'),
      need(input$num != '', 'Please enter a chip ID!'),
      need(
        input$A1 != '' &&
          input$A2 != '' &&
          input$B1 != '' &&
          input$B2 != '' &&
          input$C1 != '' && input$C2 != '' && input$D1 != '' &&
          input$D2 != '',
        'Please indicate all samples!'
      )
      
    )
    tagList(actionButton("go", "Load Dataset"))
  })
  
  
  
  ### compute other tabs only after pressing the "Load dataset" button ###
  
  observeEvent(input$go, {
    output$exp_name <- renderPrint(path())
    csv_file <- file.path(path(), input$csv_file)
    name_csv_file <- gsub(".csv", "", csv_file)
    cell_type <- input$cell_type
    if (cell_type == "MCF10CA") {
      slope <- 42.97
      intercept <- 2409.59
    } else if (cell_type == "HD1495") {
      slope <- 40.85
      intercept <- 1706.35
    } else{
      stop("Please select a valid cell type for cell number approximation.")
    }
    
    
    
    ds <-
      read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)
    ds$row <-
      as.integer(gsub("Image_(\\d*)_.+", "\\1", ds$image_name))
    ds$col <-
      as.integer(gsub("Image_.*_(\\d*)", "\\1", ds$image_name))
    ds <- ds[with(ds, order(row, col)),]
    
    
    ########################################################################
    
    
    ################ Make up grid for sample annotation ####################
    #1 or 2 after the sample letter (ex. 1 for sample A1)
    colranges_1 <- c(0:11, 24:35, 48:59)
    colranges_2 <- c(12:23, 36:47, 60:71)
    subset_1 <-
      ds[apply(ds [c("col")], 1, function(x)
        any (x %in% colranges_1)),]
    subset_2 <-
      ds[apply(ds [c("col")], 1, function(x)
        any (x %in% colranges_2)),]
    
    #the letter before the number (ex.A for sample A1)
    rowranges_A <- c(0:8, 36:44)
    rowranges_B <- c(9:17, 45:53)
    rowranges_C <- c(18:26, 54:62)
    rowranges_D <- c(27:35, 63:71)
    
    #individual subsetting of letters based on rowranges and colranges
    #Yes, this can surely be done in a more elegant way ...
    subset_1$sample <- NA
    subset_2$sample <- NA
    #A1
    subset_1[subset_1$row %in% rowranges_A,]$sample <- input$A1
    #A2
    subset_2[subset_2$row %in% rowranges_A,]$sample <- input$A2
    #B1
    subset_1[subset_1$row %in% rowranges_B,]$sample <- input$B1
    #B2
    subset_2[subset_2$row %in% rowranges_B,]$sample <- input$B2
    #C1
    subset_1[subset_1$row %in% rowranges_C,]$sample <- input$C1
    #C2
    subset_2[subset_2$row %in% rowranges_C,]$sample <- input$C2
    #D1
    subset_1[subset_1$row %in% rowranges_D,]$sample <- input$D1
    #D2
    subset_2[subset_2$row %in% rowranges_D,]$sample <- input$D2
    
    #put datasets together and order
    ds.2 <- rbind(subset_1, subset_2)
    ds.2 <- ds.2[with(ds.2, order(row, col)),]
    
    
    #Scale Major Axis to micrometer
    ds.2$Major.Axis <- as.integer(ds.2$Major.Axis / 0.3956)
    ds.2$Minor.Axis <- as.integer(ds.2$Minor.Axis / 0.3956)
    ds.2$Area <- ds.2$Major.Axis * ds.2$Minor.Axis
    
    
    # #Set Circularity to numeric, this is needed for correct plotting/selection
    # ds.2$Circularity <- as.numeric(ds.2$Circularity)
    
    #chip number
    chip_number <- input$num
    #########################################################################
    
    # #exclude wells that contain cells close to the border (Dimension.0 or Dimension.1 < 20)
    # wells_to_exclude <- ds.2[which(ds.2$Dimension.0 < 20 | ds.2$Dimension.1 < 20),"image_name"]
    # ds.2 <- subset(ds.2, !(ds.2$image_name %in% wells_to_exclude))
    
    
    #negative and positive controls
    controls_rows <- c(4, 13, 22, 31, 40, 49, 58, 67)
    negativecontroles_cols <- c(2, 14, 26, 38, 50, 62)
    positivecontroles_cols <- c(9, 21, 33, 45, 57, 69)
    negative <- merge(controls_rows, negativecontroles_cols)
    positive <- merge(controls_rows, positivecontroles_cols)
    names_negative <- NULL
    names_positive <- NULL
    
    for (i in seq(length(negative[, 1]))) {
      names_negative_i <-
        paste0("Image_", negative[i, "x"], "_", negative[i, "y"])
      names_positive_i <-
        paste0("Image_", positive[i, "x"], "_", positive[i, "y"])
      names_negative <- rbind(names_negative, names_negative_i)
      names_positive <- rbind(names_positive, names_positive_i)
    }
    
    rownames(names_negative) <- NULL
    colnames(names_negative) <- "image_name"
    names_negative <-
      data.frame(names_negative, stringsAsFactors = FALSE)
    names_negative$row <-
      as.numeric(gsub("Image_(.*)_.+", "\\1", names_negative$image_name))
    names_negative$col <-
      as.numeric(gsub("Image_.*_(.*)", "\\1", names_negative$image_name))
    names_negative <- merge(names_negative, ds.2, all.x = TRUE)
    names_negative$sample <- "Negative control"
    names_negative[is.na(names_negative)] <- 0
    
    rownames(names_positive) <- NULL
    colnames(names_positive) <- "image_name"
    names_positive <-
      data.frame(names_positive, stringsAsFactors = FALSE)
    names_positive$row <-
      as.numeric(gsub("Image_(.*)_.+", "\\1", names_positive$image_name))
    names_positive$col <-
      as.numeric(sub("Image_.*_(.*)", "\\1", names_positive$image_name))
    names_positive <- merge(names_positive, ds.2, all.x = TRUE)
    names_positive$sample <- "Positive control"
    names_positive[is.na(names_positive)] <- 0
    
    
    #take controls out
    ds.2 <- ds.2[!ds.2$image_name %in% names_positive$image_name,]
    ds.2 <- ds.2[!ds.2$image_name %in% names_negative$image_name,]
    
    #merge controls
    ds.2 <- rbind(ds.2, names_negative, names_positive)
    
    #load dummy welllist file for barcodes
    dummy_welllist <-
      read.csv(file.path(path(), "dummy_WellList.TXT"), sep = "\t")
    dummy_welllist$image_name <-
      paste0("Image_", dummy_welllist$Row, "_", dummy_welllist$Col)
    ds.2 <-
      merge(dummy_welllist[, c("image_name", "Barcode")], ds.2, by = "image_name")
    
    
    
    
    #add approximated cellnumber column
    ds.2$estimated_cell_number <-
      round((ds.2$Area - intercept) / slope, 0)
    ds.2[ds.2$estimated_cell_number <= 0,]$estimated_cell_number <-
      1
    
    
    #sample names without controls
    sample_names <-
      unique(ds.2[!ds.2$sample %in% c("Negative control", "Positive control"),]$sample)
    
    #save complete ds.2 as ds.2_all
    ds.2_all <- ds.2
    
    #only 3 parameters for ds.2 for convenient display in app
    toMatch <-
      c(
        "image_name",
        "row",
        "col",
        "Area",
        "Circularity",
        "sample",
        "estimated_cell_number"
      )
    toMatch <- gsub("(.+)", "^\\1$", toMatch, perl = TRUE)
    ds.2 <-
      subset(ds.2, select = grep(paste(toMatch, collapse = "|"), colnames(ds.2), value =
                                   TRUE))
    #change column order
    ds.2 <-
      ds.2[c(
        "image_name",
        "row",
        "col",
        "Area",
        "Circularity",
        "sample",
        "estimated_cell_number"
      )]
    
    #add comment column
    comment <- "Empty"
    ds.2 <- data.frame(ds.2, comment, stringsAsFactors = FALSE)
    
    #only wells with single spheroids
    duplicates <-
      ds.2$image_name[which(duplicated(ds.2$image_name))]
    ds.2s <- ds.2[!ds.2$image_name %in% duplicates, ]
    
    
    ########################################## RenderUI depending on input in first Tab ###########################
    output$choose_samples <- renderUI({
      tagList(
        checkboxGroupInput(
          'included_samples',
          'Included samples:',
          sample_names,
          selected = sample_names
        )
      )
    })
    
    
    output$choose_Area <- renderUI({
      tagList(sliderInput(
        "selectedbins_Area",
        "Selected size:",
        min = 0,
        max = max(round(ds.2$Area)),
        value = c(0, max(round(ds.2$Area)))
      ))
    })
    
    output$choose_circularity <- renderUI({
      tagList(
        sliderInput(
          "selectedbins_Circularity",
          "Selected circularity:",
          min = 0,
          max = 1,
          value = c(0, 1)
        )
      )
    })
    
    output$sample_name_and_chip_number <- renderUI({
      tagList(
        h4(renderText("Name of the Experiment:")),
        h4(renderText(basename(path(
          
        )))),
        h4(renderText("Chip Number:")),
        h4(renderText(chip_number)),
        br()
        
      )
    })
    
    
    ###############################################################################################################
    
    
    
    
    ########################################## functions ##########################################################
    dispense_map <- function(datafr) {
      ggplot(data = datafr, aes(x = col, y = row)) +
        geom_point(
          data = expand.grid(seq(0, 71), seq(0, 71)),
          aes(x = Var1, y = Var2),
          color = "grey90",
          shape = 15,
          size = 1
        ) +
        geom_point(aes(color = sample), size = 2, shape = 15) +
        scale_shape_manual(values = c(18, 8)) +
        scale_y_reverse(breaks = seq(0, 71, 5)) +
        scale_x_continuous(breaks = seq(0, 71, 5)) +
        labs(paste0(
          length(selected_wells$df$image_name),
          "_wells_selected"
        )) +
        theme_classic() +
        coord_fixed()
      
    }
    
    highlightbars <-
      function(datafr,
               parameter,
               numbbins,
               value1,
               value2,
               col.value,
               col = NA,
               ...) {
        combined <- paste0(deparse(substitute(datafr)), "$", parameter)
        hst <- ggplot(datafr, aes(eval(parse(text = parameter)))) +
          geom_histogram(bins = numbbins, na.rm = TRUE) +
          xlim(0, max(eval(parse(text = combined))))
        hstdat <- ggplot_build(hst)
        breaks <- hstdat$data[[1]]$x
        idxmin <- findInterval(value1, breaks)
        idxmax <- findInterval(value2, breaks)
        cols <- rep(col, length(breaks))
        cols[idxmin:idxmax] <- col.value
        ggplot(datafr, aes(eval(parse(text = parameter)))) +
          geom_histogram(bins = numbbins,
                         na.rm = TRUE,
                         col = cols) +
          geom_histogram(bins = numbbins,
                         na.rm = TRUE,
                         fill = cols) +
          xlim(0, max(eval(parse(text = combined)))) +
          ylim(0, max(hstdat$data[[1]]$count)) +
          xlab(parameter) +
          theme_bw() +
          coord_cartesian(xlim = ranges$x, ylim = ranges$y)
        
      }
    
    highlightdot <- function(datafr,parameter, tag){
      ggplot(data=datafr, aes_string(x="sample", y=parameter))+
        geom_jitter(data=datafr, color="skyblue1",alpha=0.3)+
        geom_jitter(data=tag, color="tomato", size=3)+
        geom_boxplot(color = "black", outlier.color =NA, alpha=0.1)+
        theme_bw()
    }
    
    highlightdot_new <-
      function(datafr,
               parameter,
               ...) {
        ggplot(data = datafr, aes_string(x = "sample", y = parameter)) +
          geom_boxplot(color = "darkgreen",
                       outlier.color = NA,
                       alpha = 0.1)+
          geom_jitter(aes(colour=highlight),alpha=0.7, size=3) +
          scale_colour_manual(values=c("skyblue1","tomato"))+
          theme_bw()+
          theme(legend.position="none")
      }
    
    
    interactive_tsne <- function(featbartsne, parameter) {
      mid <-
        mean(as.numeric(eval(parse(
          text = paste0("featbartsne$", parameter)
        ))))
      ggplot(featbartsne, aes(xtsne, ytsne, color = eval(parse(text = parameter)))) +
        geom_point(size = 7) +
        theme_classic() +
        labs(x = "x-tsne", y = "y-tsne") +
        scale_color_gradient2(low = "blue",
                              high = "red",
                              midpoint = mid) +
        theme(legend.text = element_text(size = 20),
              legend.title = element_text(size = 0)) +
        scale_y_reverse()+
        coord_cartesian(xlim = ranges_tsne$x, ylim = ranges_tsne$y)+
        ggtitle(paste("Parameter:",parameter, "||", chip_number,"||",length(featbartsne[,1])," spheroids"))
    }
    ################################################################################################################
    
    
    ##################################### reactive values ##########################################################
    #### reactive data frame for selected wells ######
    selected_wells <- reactiveValues()
    selected_wells$df <-
      data.frame(
        image_name = numeric(0),
        row = numeric(0),
        col = numeric(0),
        Area = numeric(0),
        Circularity = numeric(0),
        sample = numeric(0),
        comment = numeric(0)
      )
    #### reactive data frame for filter file ######
    filter <- reactiveValues(df = NULL)
    #### reactive data frame for tsne file ######
    tsne_file <- reactiveValues()
    tsne_file$df <-
      data.frame(
        xtsne = numeric(0),
        ytsne = numeric(0),
        name = numeric(0),
        sample = numeric(0)
      )
    #### reactive ranges for zoom in plot #####
    ranges <- reactiveValues(x = NULL, y = NULL)
    ranges_tsne <- reactiveValues(x = NULL, y = NULL)
    ################################################################################################################
    
    
    ################################################ TABLES ########################################################
    
    ### table with wells to select (cotanining spheroids) ####
    output$wells <- DT::renderDataTable({
      if (input$cells_or_controls == "cells") {
        dss <- eval(parse(text = input$single_or_all))
        selected_ds <- dss[dss$sample %in% input$included_samples, ]
        value1_size <- input$selectedbins_Area[1]
        value2_size <- input$selectedbins_Area[2]
        value1_circularity <- input$selectedbins_Circularity[1]
        value2_circularity <- input$selectedbins_Circularity[2]
        dss_selected <-
          selected_ds[which(
            selected_ds$Area > value1_size &
              selected_ds$Area < value2_size &
              selected_ds$Circularity > value1_circularity &
              selected_ds$Circularity < value2_circularity
          ),]
        dss_selected$Circularity <-
          round(dss_selected$Circularity, 2)
        dss_selected$Area <- round(dss_selected$Area, 2)
        DT::datatable(dss_selected, select = "single")
      } else{
        dss_selected <-
          ds.2[ds.2$sample %in% c("Positive control", "Negative control"),]
        dss_selected$Circularity <-
          round(dss_selected$Circularity, 2)
        dss_selected$Area <- round(dss_selected$Area, 2)
        DT::datatable(dss_selected, select = "single")
      }
    })
    
    ### table with selected wells###
    output$selected_table <-  DT::renderDataTable({
      if (input$cells_or_controls_selected == "cells") {
        dss_selected <-
          selected_wells$df[selected_wells$df$sample %in% sample_names,]
        dss_selected$Circularity <-
          round(dss_selected$Circularity, 2)
        dss_selected$Area <- round(dss_selected$Area, 2)
      } else{
        dss_selected <-
          selected_wells$df[selected_wells$df$sample %in% c("Positive control", "Negative control"),]
        dss_selected$Circularity <-
          round(dss_selected$Circularity, 2)
        dss_selected$Area <- round(dss_selected$Area, 2)
      }
      DT::datatable(dss_selected, select = "single")%>% formatStyle(0, cursor = 'pointer')
    })
    
    ################################################################################################################
    
    
    ################################################ BUTTONS #######################################################
    
    ### Button to select wells for RT dispense ###
    #selection based on image_name, so that all spheroids contained in one well will appear
    observeEvent(input$select_button_single, {
      if (!is.null(input$wells_rows_selected)) {
        if (input$cells_or_controls == "cells") {
          dss <- eval(parse(text = input$single_or_all))
          selected_ds <-
            dss[dss$sample %in% input$included_samples, ]
          value1_size <- input$selectedbins_Area[1]
          value2_size <- input$selectedbins_Area[2]
          value1_circularity <- input$selectedbins_Circularity[1]
          value2_circularity <- input$selectedbins_Circularity[2]
          dss_selected <-
            selected_ds[which(
              selected_ds$Area > value1_size &
                selected_ds$Area < value2_size &
                selected_ds$Circularity > value1_circularity &
                selected_ds$Circularity < value2_circularity
            ),]
          selected_row_name_in_selected_ds <-
            isolate(dss_selected[input$wells_rows_selected, 'image_name'])
          selected_row_in_selected_ds <-
            isolate(dss_selected[dss_selected$image_name %in% selected_row_name_in_selected_ds,])
          duprow <-
            selected_row_name_in_selected_ds %in% selected_wells$df$image_name
        } else{
          dss_selected <-
            ds.2[ds.2$sample %in% c("Positive control", "Negative control"),]
          selected_row_name_in_selected_ds <-
            isolate(dss_selected[input$wells_rows_selected, 'image_name'])
          selected_row_in_selected_ds <-
            isolate(dss_selected[dss_selected$image_name %in% selected_row_name_in_selected_ds,])
          duprow <-
            selected_row_name_in_selected_ds %in% selected_wells$df$image_name
        }
        if (duprow == FALSE) {
          selected_wells$df <-
            rbind(selected_wells$df, selected_row_in_selected_ds)
        } else{
          showModal(
            modalDialog(
              title = "Aldeady selected",
              "This well is already selected.",
              easyClose = TRUE
            )
          )
        }
      } else{
        showModal(modalDialog(
          title = "Select a well",
          "Please select a well!",
          easyClose = TRUE
        ))
      }
    })
    
    
    ### Button to unselect wells for RT dispense ###
    observeEvent(input$unselect_button_single, {
      if (input$cells_or_controls_selected == "cells") {
        unselected_row_in_selected_ds <-
          isolate(input$selected_table_rows_selected)
        unselected_row_name_in_selected_ds <-
          selected_wells$df[unselected_row_in_selected_ds, 'image_name']
      } else{
        dss_selected <-
          selected_wells$df[selected_wells$df$sample %in% c("Positive control", "Negative control"),]
        unselected_row_in_selected_ds <-
          isolate(input$selected_table_rows_selected)
        unselected_row_name_in_selected_ds <-
          dss_selected[unselected_row_in_selected_ds, 'image_name']
      }
      if (!is.null(input$selected_table_rows_selected)) {
        selected_wells$df <-
          selected_wells$df[-c(grep(
            unselected_row_name_in_selected_ds,
            selected_wells$df$image_name
          )),]
      } else {
        showModal(modalDialog(
          title = "Select a well",
          "Please select a well!",
          easyClose = TRUE
        ))
      }
    })
    
    ### Button to reset selection for RT dispense ###
    observeEvent(input$unselect_button_all, {
      selected_wells$df <-
        data.frame(Column1 = numeric(0), Column2 = numeric(0))
    })
    
    
    ### Button to add a comment ###
    observeEvent(input$add_comment_button, {
      if (input$comment != "Add a comment.") {
        selected_row_in_selected_ds <-
          isolate(input$selected_table_rows_selected)
        selected_row_name_in_selected_ds <-
          selected_wells$df[selected_row_in_selected_ds, 'image_name']
        selected_wells$df$comment <-
          as.character(selected_wells$df$comment)
        selected_wells$df[selected_wells$df$image_name %in% selected_row_name_in_selected_ds, "comment"] <-
          input$comment
      } else{
        showModal(
          modalDialog(
            title = "Enter comment",
            "Please enter a comment!",
            easyClose = TRUE
          )
        )
      }
    })
    
    
    ### Button to save selected wells ###
    output$save_selection <- downloadHandler(
      filename = function () {
        paste0(
          chip_number,
          "_",
          length(selected_wells$df$image_name),
          "_wells_selected",
          ".csv"
        )
      },
      content = function (file) {
        write.csv(ds.2_all[ds.2_all$image_name %in% selected_wells$df$image_name,], file, row.names =
                    FALSE)
      }
    )
    
    
    
    
    ### Button to save filter file ###
    observe({
      row <- selected_wells$df$row + 1
      col <- selected_wells$df$col + 1
      filter_matrix <- matrix(0, nrow = 72, ncol = 72)
      for (i in seq(length(col))) {
        filter_matrix[row[i], col[i]] <- 1
      }
      filter$df <- as.data.frame(filter_matrix)
    })
    output$save_filterfile <- downloadHandler(
      filename = function () {
        paste0(
          chip_number,
          "_",
          length(selected_wells$df$image_name),
          "_wells_selected_filter_file",
          ".csv"
        )
      },
      content = function (file) {
        write.table(
          filter$df,
          file,
          row.names = FALSE,
          col.names = FALSE,
          sep = ","
        )
      }
    )
    
    
    ### Button to Upload already saved filter file with selected wells ###
    observeEvent(input$uploaded_filter_file, {
      inFile <- input$uploaded_filter_file
      if (is.null(inFile))
        selected_wells$df <- NULL
      filterfile <- read.csv(inFile$datapath, header = FALSE)
      match <- as.data.frame(which(filterfile == 1, arr.ind = TRUE))
      well_list <- NULL
      for (i in seq(length(match$row))) {
        row <- match$row[i] - 1
        col <- match$col[i] - 1
        temp <- paste0("Image_", row, "_", col)
        well_list <- c(well_list, temp)
      }
      selected_rows <-
        isolate(ds.2[ds.2$image_name %in% well_list,])
      selected_wells$df <- selected_rows
    })
    
    ### Button to upload tsne file ###
    observeEvent(input$uploaded_tsne_file, {
      inFile <- input$uploaded_tsne_file
      if (is.null(inFile))
        tsne_file$df <- NULL
      tsne_file$df <- read.csv(inFile$datapath, header = TRUE)
    })
    
    ################################################################################################################
    
    ################################################ PLOTS #########################################################
    
    ### Dispense Map ###
    output$dispense_map <- renderPlot({
      if (nrow(selected_wells$df) != 0) {
        dispense_map(selected_wells$df)
      }
    })
    
    ### Histogram (using highlightbars function) ###
    output$hist_bars <- renderPlot({
      dss <- eval(parse(text = input$single_or_all))
      parameter <- input$parameters_to_plot
      input_parameter <-
        paste0("input$selectedbins", "_", parameter)
      bins <- input$numberofbins
      datafr <- dss[dss$sample %in% input$included_samples, ]
      value1 <- eval(parse(text = input_parameter))[1]
      value2 <- eval(parse(text = input_parameter))[2]
      highlightbars(datafr, parameter, bins, value1, value2, "darkred")
    })
    
    ### Zoom function in histogram plot ###
    observeEvent(input$hist_bars_dblclick, {
      brush <- input$hist_bars_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    
    ### Dotplot (using highlightdot function) ###
    output$dotplot_highlighted <- renderPlot({
      info <- input$selected_table_cell_clicked
      if (nrow(selected_wells$df) == 0 || is.null(selected_wells$df) ||is.null(info$value) || info$col != 0) return()
          if (input$cells_or_controls_selected == "cells") {
            set.seed(1)
            parameter <- input$parameters_to_plot_selection
            selected_row_in_selected_ds <-
              selected_wells$df[info$value,]
            image_name <- selected_row_in_selected_ds[, "image_name"]
            if (input$dataset_to_plot_selection == "all_cells") {
              datafr <- ds.2[ds.2$sample %in% sample_names,]
            } else if (input$dataset_to_plot_selection == "single_cell_wells") {
              datafr <- ds.2s[ds.2s$sample %in% sample_names,]
            } else if (input$dataset_to_plot_selection == "selected_wells") {
              datafr <-selected_wells$df[selected_wells$df$sample %in% sample_names,]
            }
            if (parameter == "Circularity"){
              datafr$Circularity <- datafr$Circularity*10
            }
            tag <- datafr[datafr[, "image_name"] %in% image_name,]
            highlightdot(datafr, parameter, tag)+
              ggtitle(image_name)
          }
    })
    
    
    ################################################################################################################
    
    ################################################ IMAGES ########################################################
    
    ### Images for selection page ###
    output$image_well_selection <- renderImage({
      selected_channel <- input$channel_selection
      if (input$cells_or_controls == "cells") {
        dss <- eval(parse(text = input$single_or_all))
        selected_ds <- dss[dss$sample %in% input$included_samples, ]
        value1_size <- input$selectedbins_Area[1]
        value2_size <- input$selectedbins_Area[2]
        value1_circularity <- input$selectedbins_Circularity[1]
        value2_circularity <- input$selectedbins_Circularity[2]
        dss_selected <-
          selected_ds[which(
            selected_ds$Area > value1_size &
              selected_ds$Area < value2_size &
              selected_ds$Circularity > value1_circularity &
              selected_ds$Circularity < value2_circularity
          ),]
        selected_row_in_selected_ds <-
          dss_selected[input$wells_rows_selected,]
        image_name <- selected_row_in_selected_ds[, "image_name"]
        filename_well <-
          paste0(selected_channel, "_", image_name, ".png")
        list(src = file.path(path(), "files", filename_well),
             width = 350)
      } else{
        dss_selected <-
          ds.2[ds.2$sample %in% c("Positive control", "Negative control"),]
        image_name <-
          dss_selected[input$wells_rows_selected, 'image_name']
        filename_well <-
          paste0(selected_channel, "_", image_name, ".png")
        list(src = file.path(path(), "files", filename_well),
             width = 350)
      }
    }, deleteFile = FALSE)
    
    ### Images for selected wells ###
    output$image_well_selected <- renderImage({
      info <- input$selected_table_cell_clicked
      if (nrow(selected_wells$df) == 0 || is.null(selected_wells$df) || is.null(info$value) || is.na(info$value) || info$col != 0){
        stop("Please select a well.")
      }else{
      selected_channel <- input$channel_selected
        if (input$cells_or_controls_selected == "cells") {
          selected_row_in_selected_ds <-
            selected_wells$df[info$value,]
          image_name <- selected_row_in_selected_ds[, "image_name"]
          filename_well <-
            paste0(selected_channel, "_", image_name, ".png")
          list(src = file.path(path(), "files", filename_well),
               width = 350)
        } else{
          dss_selected <-
            selected_wells$df[selected_wells$df$sample %in% c("Positive control", "Negative control"),]
          selected_row_in_selected_ds <-
            dss_selected[input$selected_table_rows_selected,]
          image_name <- selected_row_in_selected_ds[, "image_name"]
          filename_well <-
            paste0(selected_channel, "_", image_name, ".png")
          list(src = file.path(path(), "files", filename_well),
               width = 350)
        }
      }
    }, deleteFile = FALSE)
    
    
    
    ################################################################################################################
    
    ################################################ TEXTS ########################################################
    
    ### Well name on selection page ###
    output$well_name <- renderText({
      if (input$cells_or_controls == "cells") {
        dss <- eval(parse(text = input$single_or_all))
        selected_ds <- dss[dss$sample %in% input$included_samples, ]
        value1_size <- input$selectedbins_Area[1]
        value2_size <- input$selectedbins_Area[2]
        value1_circularity <- input$selectedbins_Circularity[1]
        value2_circularity <- input$selectedbins_Circularity[2]
        dss_selected <-
          selected_ds[which(
            selected_ds$Area > value1_size &
              selected_ds$Area < value2_size &
              selected_ds$Circularity > value1_circularity &
              selected_ds$Circularity < value2_circularity
          ),]
        selected_row_in_selected_ds <-
          dss_selected[input$wells_rows_selected,]
        paste0("Well: ", selected_row_in_selected_ds[, "image_name"])
      } else{
        dss_selected <-
          ds.2[ds.2$sample %in% c("Positive control", "Negative control"),]
        selected_row_in_selected_ds <-
          dss_selected[input$wells_rows_selected,]
        paste0("Well: ", selected_row_in_selected_ds[, "image_name"])
        
      }
    })
    
    
    ### Well name on selected page ###
    output$well_name_selected <- renderText({
      info <- input$selected_table_cell_clicked
      if (nrow(selected_wells$df) == 0 || is.null(selected_wells$df) ||is.null(info$value) || info$col > 0) return()
        selected_channel <- input$channel_selected
        if (input$cells_or_controls_selected == "cells") {
          selected_row_in_selected_ds <-
            selected_wells$df[info$value,]
          image_name <- selected_row_in_selected_ds[, "image_name"]
          paste0("Well: ", image_name)
        } else{
          dss_selected <-
            selected_wells$df[selected_wells$df$sample %in% c("Positive control", "Negative control"),]
          selected_row_in_selected_ds <-
            dss_selected[input$selected_table_rows_selected,]
          image_name <- selected_row_in_selected_ds[, "image_name"]
          image_name
        }
    })
    
    
    ### Number and type of selected wells for dispense map plot ###
    output$number_selected <- renderText({
      if (nrow(selected_wells$df) != 0) {
        only_unique <- subset(selected_wells$df, !duplicated(image_name))
        samples_count <- count(only_unique$sample)
        samples_names <- as.character(samples_count[, "x"])
        samples_numbers <- as.integer(samples_count[, "freq"])
        names_and_numbers <-
          paste(samples_numbers, samples_names, "wells", sep = " ")
        do.call(paste, c(as.list(names_and_numbers), sep = " and "))
      }
    })
    
    
    
    
    ################################################ Rendering of tSNE Data (independent of the rest) #########################################################
    
    output$choose_RData <- renderUI({
      tagList(
        selectInput(
          'RData_file',
          label = NULL,
          list.files(path(), pattern = ".RData"),
          multiple = FALSE,
          selectize = FALSE
        )
      )
    })
    output$choose_parameter <- renderUI({
      tagList(
        selectInput(
          'parameter',
          label = "Select a feature to be displayed in tSNE.",
          colnames(ds.2_all[,!names(ds.2_all) %in% c("Barcode","image_name","sample")]),
          multiple = FALSE,
          selectize = TRUE
        )
      )
    })

    # load dataset button, allows loading of data for the remaining tabs
    output$load_RData <- renderUI({
      validate(need(
        input$RData_file != '',
        'Please choose a RData file as exported from PAGODA'
      ))
      tagList(actionButton("go2", "Load RData File"))
    })
    
    
    
    observeEvent(input$go2, {
      RData_file <- file.path(path(), input$RData_file)
      tsne <- load(RData_file)
      tsne <- data.frame(tSNE.pagoda$Y)
      tsne$image_name <-
        paste(
          "Image",
          gsub(".+\\-R(\\d+)C.+", "\\1", rownames(tsne)),
          gsub(".+\\-R\\d+C(\\d+)", "\\1", rownames(tsne)),
          sep = "_"
        )
      featbartsne <- merge(tsne, ds.2_all, by = "image_name")
      colnames(featbartsne)[colnames(featbartsne) %in% c("X1", "X2")] <-
        c("xtsne", "ytsne")
      featbartsne$xtsne <- as.numeric(featbartsne$xtsne)
      featbartsne$ytsne <- as.numeric(featbartsne$ytsne)
      

      observeEvent(input$tsne_plot_dblclick, {
        brush <- input$tsne_plot_brush
        if (!is.null(brush)) {
          ranges_tsne$x <- c(brush$xmin, brush$xmax)
          ranges_tsne$y <- c(brush$ymin, brush$ymax)
        } else {
          ranges_tsne$x <- NULL
          ranges_tsne$y <- NULL
        }
      })
      
      observeEvent(input$parameter, {
        output$tsne_plot <- renderPlot({
          parameter <- input$parameter
          interactive_tsne(featbartsne, parameter)
        })
        output$downloadPlot <- downloadHandler(
          parameter <- input$parameter,
          filename = function () {
            paste0(
              chip_number,
              "_",
              nrow(featbartsne),
              "_spheroids_",
              parameter,
              ".png"
            )
          },
          content = function(file) {
            ggsave(file, plot = interactive_tsne(featbartsne, parameter), device = "png")
          })
      })
        
      
      
      observeEvent(input$tsne_plot_click,{
        if (input$tsne_plot_click != ''){
        output$tsne_features <- renderUI({
          tagList(
            DT::renderDataTable({
              if (!is.null(input$tsne_plot_click)) {
                features <- nearPoints(
                  featbartsne,
                  input$tsne_plot_click,
                  addDist = FALSE,
                  maxpoints = 1,
                  threshold = 5
                )[1,]
                features <- melt(features,id="image_name")
                DT::datatable(features[,c("image_name","variable","value")], select = "single")
              } else{
                stop("Select a spheroid")
              }
            })
          )
        })
        }
      })
      
      
      ## Images for tsne plot ###
      observeEvent(input$tsne_plot_click,{
        if (input$tsne_plot_click != ''){
      output$tsne_image <- renderImage({
          selected_channel <- input$tsne_channel_selected
          name <-nearPoints(
            featbartsne,
            input$tsne_plot_click,
            addDist = FALSE,
            maxpoints = 1,
            threshold = 5
          )[, "image_name"]
          filename_well <- paste0(selected_channel, "_", name, ".png")
          list(src = file.path(path(), "files", filename_well),
               width = 350)
      }, deleteFile = FALSE)
        }
      })
      
      
      
    })#closes go2
    
  })#closes go
  
  
}
shinyApp(ui, server)
