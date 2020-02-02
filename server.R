#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

shinyServer(function(input, output, session) {
  #tabPanel-Input
  ###display demo count matrix
  ctobj <- isolate(reactiveVal())
  pobj <- isolate(reactiveVal())
  
  observeEvent(input$runexample, {
    ctobj <- (NULL)
    pobj <- (NULL)
    
    set.seed(123456)
    n=2000; m=9
    ctobj(simulateRnaSeqData(output="matrix", n=n, m=m)%>%data.frame%>%
            mutate(Gene=sample(unique(hgnc.table$Approved.Symbol), n))%>%
            dplyr::select(Gene, everything(.)))
    
    pobj(data.frame(ID=paste0("sample", 1:m), 
                    Treatment=rep(c("Dose10", "Control", "Dose20"), each=3), 
                    Gender=sample(c("F", "M"), m, T)))
    
    #ngenes
    output$ngene <- renderText({paste("Number of genes: ", dim(ctobj())[1], " [First 10 rows displayed]")})
    #nsamples
    output$nsample <- renderText({paste("Number of samples: ", (dim(ctobj())[2])-1, " [First 10 rows displayed]")})
    #display 10rows count matrix
    output$matrix <- renderTable({
      head(ctobj(), 10)
    })
    #display10rows manifest
    output$pdat <- renderTable({
      head(pobj(), 10)
    })
    #model variables
    ##comparison variable
    observe({
      updateSelectInput(session, "design", choices="Treatment")
    })
    ##ref0
    output$level0 <- renderUI({
      selectInput("ref0", "Reference group", "Control")
    })
    ##ref1
    output$level1 <- renderUI({
      selectInput("ref1", "Study group", "Dose20")
    })
    ##species
    observe({
      updateSelectInput(session, "species", choices="Human")
    })
  })
  
  observeEvent(input$file1, {
    ctobj <- (NULL)
    ctobj(read_excel(input$file1$datapath))
    
    ##SHOW SUMMARY
    output$ngene <- renderText({paste("Number of genes: ", dim(ctobj())[1], ". [First 10 rows displayed]")})
    
    output$nsample <- renderText({paste("Number of samples: ", (dim(ctobj())[2])-1, ". [First 10 rows displayed]")})
    
    ##DISPLAY 10 ROWS
    output$matrix <- renderTable({
      head(ctobj(), 10)
    })
  })

  observeEvent(input$file2, {
    pobj <- (NULL)
    pobj(read_excel(input$file2$datapath))
    
    output$pdat <- renderTable({
      head(pobj(),  10)
    })
    ##MODEL VARIABLES
    ###COMPARISON VARIALBE
    observe({
      updateSelectInput(session, "design", choices=names(pobj()))
    })
    ###CONTROL
    output$level0 <- renderUI({
      selectInput("ref0", "Reference group", pobj()[[input$design]])
    })
    ###TARGET
    output$level1 <- renderUI({
      selectInput("ref1", "Study group", pobj()[[input$design]])
    })
  })
  
  ##ANALYSIS
  resobj <- reactiveVal()
  volplot <- reactiveVal()
  gores <- reactiveVal()
  
  observeEvent(input$runbutton, {
    resobj <- (NULL)
    volplot <- (NULL)
    gores <- (NULL)
    withProgress(message = 'Running ...', value=0, style = "old",{
      ###DGE
      Sys.sleep(1)
      resobj(getdgeres(ctobj(), pobj(),
                       comparison=input$design, 
                       level1=input$ref1,
                       level0=input$ref0))
      incProgress(0.4, detail="Differeital expression ... ")
      ###PATHWAY
      Sys.sleep(2)
      gores(getgores(resobj()[["results"]], species=input$species))
      incProgress(0.4, detail="Pahtway analysis ... ")
      
      ###PLOTTING
      Sys.sleep(3)
      volplot(getvolcano(resobj()[["results"]]))
      incProgress(0.2, detail="making plot ... ")
    })
  })
  
  ##tabPanel-RESULTS
  todowndat <- reactive({
    switch(input$results,
           "Results" = resobj()[["results"]],
           "Normalized matrix" = resobj()[["normal"]]
    )
  })
  
  output$table <- renderTable({
    todowndat()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$results, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(todowndat(), file, row.names = FALSE)
    }
  ) 
  
  ##tabPanel-PLOT
  output$plot1 <- renderPlot({
    volplot()[["plot"]]
  })
  
  output$brush_info <- renderPrint({
    showdf<-volplot()[["brush"]]%>%dplyr::select(Gene, log2FoldChange, pvalue, padj, log10padj) 
    brushedPoints(showdf, input$plot1_brush)
  })
  
  ##tabPanel-PATHWAY
  godowndat <- reactive({
    switch(input$gopathway,
           "Biological process greater"=gores()[["bpupper"]],
           "Biological process less"=gores()[["bpless"]],
           "Molecular function greater"=gores()[["mfupper"]],
           "Molecular function less"=gores()[["mfless"]],
           "Cellular component greater"=gores()[["ccupper"]],
           "Cellular component less"=gores()[["ccless"]],
           "KEGG greater"=gores()[["kgupper"]],
           "KEGG less"=gores()[["kgless"]])
  })
  
  output$gores <- renderTable({
    godowndat()
  })
  
  output$downloadGo <- downloadHandler(
    filename = function() {
      paste(input$gopathway, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(todowndat(), file, row.names = FALSE)
    }
  ) 
  #RESET for new analysis
  observeEvent(input$reset, {js$reset()}) 
})
