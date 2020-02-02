#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

jsResetCode <- "shinyjs.reset = function() {history.go(0)}"
shinyUI(
   fluidPage(
  #set backgroud color
  setBackgroundColor(
    color = c("#F7FBFF", "#2171B5"),
    gradient = "radial",
    direction = c("top", "left")
  ),
  
  #reset session by reset button
  useShinyjs(),                                           # Include shinyjs in the UI
  extendShinyjs(text = jsResetCode, functions="shinyjs.reset"), 
  
  #panels
  tabsetPanel(
    ##tabPanel-Input
    tabPanel("Input", fluid = TRUE,
             
             # tab title ----
             titlePanel("Upload data"),
             
             # sidebar layout with input and output tables ----
             sidebarLayout(
               
               # sidebar panel for inputs ----
               sidebarPanel(
                 #show ct demo
                 actionBttn("runexample", "Import demo data", style="simple", size="sm", color = "primary"),
                 
                 # input1: Select a file ----
                 fileInput("file1", "Count matrix File (.xlsx)",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 #input2: select a file ----
                 fileInput("file2", "Manifest File (.xlsx)",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 #select column name
                 selectInput("design", "Column name for analysis", " "),
                 
                 #select ref group
                 uiOutput("level0"),
                 
                 #select study group
                 uiOutput("level1"),
                 
                 #select column name
                 selectInput("species", "Species", c("Human"="Human", "Mouse"="Mouse")),
                 
                 #action run
                 actionBttn("runbutton", "GO", style="simple", size="sm", color = "primary"),
                 
                 actionBttn("reset", "RESET", style="simple", size="sm", color = "warning"),
                 
                 #comment message
                 helpText("For demonstration, click `Import demo data` and `GO`"),
                 helpText("Click `RESET` to upload your data and `GO` to launch analyses"),
                 
                 #README link
                 h6("Need support in generating count matrix?"),
                 a(actionBttn(inputId = "email1", 
                              label = "Contact developer", 
                              icon = icon("envelope", lib = "font-awesome"),
                              size="xs",
                              color="success"),
                   href="mailto:doubleomics@gmail.com"),
                 
                 a(actionBttn(inputId = "twitter_share",
                              label = "Follow us",
                              icon = icon("twitter", lib = "font-awesome"),
                              size="xs",
                              color="success"),
                   href = "https://twitter.com/DoubleOmics")
                 
               ),
               # Main panel for displaying outputs ----
               mainPanel(
                 # Output: Data file ----
                 span(textOutput("ngene"),style="color:blue"),
                 span(textOutput("nsample"),style="color:blue"),
                 tableOutput("matrix"),
                 tableOutput("pdat")
               )
             )
    ),
    
    #tabPanel-Results
    tabPanel("DGE results", fluid = TRUE,
             # App title ----
             titlePanel("Download results"),
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               # Sidebar panel for inputs ----
               sidebarPanel(
                 # Input: Choose dataset ----
                 selectInput("results", "Choose a dataset:",
                             choices = c("Results", "Normalized matrix")),
                 # Button
                 downloadButton("downloadData", "Download")
                 
               ),
               # Main panel for displaying outputs ----
               mainPanel(
                 tableOutput("table")
               )
             )
    ),
    
    #tabPanel-Plots
    tabPanel("Volcano plot", fluid = TRUE,
             fluidRow(
               column(width = 8,
                      plotOutput("plot1", height = 800,
                                 # Equivalent to: click = clickOpts(id = "plot_click")
                                 click = "plot1_click",
                                 brush = brushOpts(
                                   id = "plot1_brush"
                                 )
                      )
               ),
               column(width = 4,
                      h4("Brushed points"),
                      verbatimTextOutput("brush_info")
               )
             )
    ),
    
    #tabPanel-GO Results
    tabPanel("Pathways results", fluid = TRUE,
             # App title ----
             titlePanel("Download results"),
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               # Sidebar panel for inputs ----
               sidebarPanel(
                 # Input: Choose dataset ----
                 selectInput("gopathway", "Choose a dataset:",
                             choices = c("Biological process greater",
                                         "Biological process less",
                                         "Molecular function greater",
                                         "Molecular function less",
                                         "Cellular component greater",
                                         "Cellular component less",
                                         "KEGG greater",
                                         "KEGG less")),
                 # Button
                 downloadButton("downloadGo", "Download")
               ),
               # Main panel for displaying outputs ----
               mainPanel(
                 tableOutput("gores")
               )
             )             
    )
  )
)
