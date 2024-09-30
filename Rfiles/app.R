# scRNA Analyser

# Imports functions
source('global.R')

# Define User Interface for the app
ui <- dashboardPage(skin = 'purple',
                    dashboardHeader(title = "scRNAseq Analysis"), 
                    
                    # Sidebar Layout for Navigation
                    dashboardSidebar(
                      # Hides the sidebar toggle button
                      tags$head(
                        tags$style(HTML(".skin-blue .main-header .sidebar-toggle {display: none;}"))
                      ),
                      
                      # Sidebar menu with navigation links
                      sidebarMenu(id='tab',
                                  useShinyjs(),
                                  # Home page tab
                                  menuItem("Home Page", tabName = "home", icon = icon("list")),  
                                  # Analyzer tab
                                  menuItem("scRNAseq Analyzer", tabName = "input", icon = icon("edit")),  
                                  
                                  # Conditional panel: Shown only if the scRNAseq Analyser tab is active
                                  conditionalPanel(condition = "input.tab == 'input'",
                                                   div(
                                                     # File upload for .rds files
                                                     fileInput("file", "Upload File", multiple=TRUE, accept=c('.rds')),  
                                                     # Reset button to clear inputs
                                                     actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                                                     # Run button to start analysis
                                                     actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                                                   )
                                  )
                      )
                    ), 
                    
                    # Main body of the app containing different tabs
                    dashboardBody(
                      tabItems(
                        # scRNAseq Analyser tab content
                        tabItem(tabName = "input", 
                                tabsetPanel(id = 'main_tabs',
                                            # Loads markdown file with instructions
                                            tabPanel("Instructions",
                                                     includeMarkdown("./markdown/instructions.md")
                                            )
                                )
                        ),
                        
                        # Home Page content
                        # Link to preprocessing tutorial on GitHub
                        tabItem(tabName = "home",
                                tags$h1(HTML("<b>Welcome to the scRNAseq Suerat analysis RShiny app</b>")),
                                tags$a(href="https://github.com/LouisCHague/scRNAseq-Analyser-Project/blob/main/preprocessing_tutorial.R", 
                                       "Preprocessing Steps") 
                        )
                      )
                    )         
)

# Server logic for the app
server <- function(input, output, session) {
  # Allow file uploads up to 300MB
  options(shiny.maxRequestSize = 300 * 1024^2)
  
  values <- reactiveValues()  # Store reactive values
  
  # Disable the "Run" button by default until a file is uploaded
  shinyjs::disable("run")
  
  # Observer: Enable the "Run" button when a file is uploaded, disable otherwise
  observe({
    if (is.null(input$file) != TRUE) {
      shinyjs::enable("run")  # Enable button if a file is selected
    } else {
      shinyjs::disable("run")  # Disable button if no file is selected
    }
  })
  
  # Handles what happens when the "Run" button is clicked
  observeEvent(input$run, {
    shinyjs::disable("run")  # Disable the run button to prevent re-running until complete
    
    # Remove previously generated tabs (UMAP, Gene Expression) if "Run" is clicked again
    removeTab("main_tabs", "UMAP")
    removeTab("main_tabs", "Gene Expression")
    
    # Show a modal spinner while processing the data
    show_modal_spinner(text = "Preparing plots...")
    
    # Load the Seurat object from the uploaded file
    obj <- load_seurat_obj(input$file$datapath)
    
    # Check for errors in the uploaded file
    if (is.vector(obj)) {  # If there's an error, show a modal with the error details
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
             paste(unlist(obj), collapse = "<br><br>"))
      ))
      shinyjs::enable("run")  # Re-enable the run button if there was an error
      
    } else {
      # If the file is correct, generate plots and add new tabs to display them
      
      # Render the UMAP plot (2D visualization of cells)
      output$umap <- renderPlot({
        if (!is.null(input$metadata_col)) {
          create_metadata_UMAP(obj, input$metadata_col)  
        }
      })
      
      # Render the feature plot Gene expression plot
      output$featurePlot <- renderPlot({
        if (!is.null(input$gene)) {
          create_feature_plot(obj, input$gene) 
        }
      })
      
      # File download handler for the feature plot
      output$downloadFeaturePlot <- downloadHandler(
        filename = function(){
          paste0(input$gene, '_feature_plot', '.png')  
        },
        content = function(file){
          plot <- create_feature_plot(obj, input$gene) 
          ggsave(filename = file, width = 10, height = 5, type = "cairo")
        }
      )
      
      # File download handler for the UMAP plot
      output$download_umap <- downloadHandler(
        filename = function(){
          paste0(input$metadata_col, '_UMAP', '.png') 
        },
        content = function(file){
          plot <- create_metadata_UMAP(obj, input$metadata_col)
          ggsave(filename = file, width = 10, height = 5, type = "cairo")
        }
      )
      
      # Add new "UMAP" tab to display the UMAP plot and input options
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "UMAP",
          fluidRow(
            column(
              width = 8,
              plotOutput(outputId = 'umap'),  # Display UMAP plot
              downloadButton("download_umap", "Download UMAP")  # Download button for UMAP
            ),
            column(
              width = 4,
              selectizeInput("metadata_col",  # Dropdown to select metadata column for UMAP
                             "Metadata Column", 
                             colnames(obj@meta.data)  # Get metadata columns from the Seurat object
              )
            )
          ),
          style = "height: 90%; width: 95%; padding-top: 5%;"
        ),
        select = TRUE  
      )
      
      # Add new "Gene Expression" tab to display the feature plot and input options
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "Gene Expression",
          fluidRow(
            column(
              width = 8,
              plotOutput(outputId = 'featurePlot'),  # Display feature plot
              downloadButton("downloadFeaturePlot", "Download Feature Plot")  # Download button for feature plot
            ),
            column(
              width = 4,
              selectizeInput("gene",  # Dropdown to select gene for feature plot
                             "Genes", 
                             rownames(obj)  # Get gene names from Seurat object
              )
            )
          ),
          style = "height: 90%; width: 95%; padding-top: 5%;"
        )
      )
      
      # Remove spinner after all plots are created
      remove_modal_spinner()
      
      # Re-enable the run button for subsequent analysis
      shinyjs::enable("run")
    }
  })
  
  # Observer: Handles the reset button, clears all inputs and generated plots
  observeEvent(input$reset, {
    shinyjs::reset("file")  # Reset the file input
    removeTab("main_tabs", "UMAP")  # Remove UMAP tab if present
    removeTab("main_tabs", "Gene Expression")  # Remove Gene Expression tab if present
    shinyjs::disable("run")  # Disable the run button after reset
  })
  
}

# Launch the app
shinyApp(ui, server)
