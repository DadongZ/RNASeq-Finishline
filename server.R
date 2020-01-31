#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
shinyServer(function(input, output, session) {
  
  #tabPanel-Input
  
  ###demo data
  ####count
  set.seed(123456)
  n=10000; m=9
  ctdemo<- simulateRnaSeqData(output="matrix", n=n, m=m)%>%data.frame%>%
      mutate(Gene=sample(unique(hgnc.table$Approved.Symbol), n))%>%
      dplyr::select(Gene, everything(.))
  ####manifest
  pdemo<-data.frame(ID=paste0("sample", 1:m), 
                    Treatment=rep(c("Dose10", "Control", "Dose20"), each=3), 
                    Gender=sample(c("F", "M"), m, T))
  ###display demo count matrix
  observeEvent(input$runexample, {
    output$matrix <- renderTable({
      head(ctdemo, 10)
    })
    output$pdat <- renderTable({
      head(pdemo, 10)
    })
    
    observe({
      updateSelectInput(session, "design", choices="Treatment")
    })
    
    output$level0 <- renderUI({
      selectInput("ref0", "Reference group", "Control")
    })
    
    output$level1 <- renderUI({
      selectInput("ref1", "Study group", "Dose20")
    })
    
    observe({
      updateSelectInput(session, "species", choices="Human")
    })
    
    resdemo <- reactiveVal()
    voldemo <- reactiveVal()
    godemo <- reactiveVal()
    observeEvent(input$runbutton, {
      resdemo <- (NULL)
      voldemo <- (NULL)
      godemo <- (NULL)
      withProgress(message = 'Running...', style="old", value=0, {
        Sys.sleep(1)
        resdemo(getdgeres(ctdemo, pdemo,
                       comparison=input$design, 
                       level1=input$ref1,
                       level0=input$ref0))
        incProgress(0.4, detail="Differeital expression ... ")
        
        Sys.sleep(2)
        godemo(getgores(resdemo()[["results"]], species=input$species))
        incProgress(0.4, detail="Pahtway analysis ... ")
        
        Sys.sleep(3)
        voldemo(getvolcano(resdemo()[["results"]]))
        incProgress(0.2, detail="making plot ... ")
      })
    })
    todowndat <- reactive({
      switch(input$results,
             "Results" = resdemo()[["results"]],
             "Normalized matrix" = resdemo()[["normal"]]
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
    #tabPanel-Results
    output$plot1 <- renderPlot({
      voldemo()[["plot"]]
    })
    
    output$brush_info <- renderPrint({
      showdf<-voldemo()[["brush"]]%>%dplyr::select(Gene, log2FoldChange, pvalue, padj, log10padj) 
      brushedPoints(showdf, input$plot1_brush)
    })
    
    godowndat <- reactive({
      switch(input$gopathway,
             "Biological process greater"=godemo()[["bpupper"]],
             "Biological process less"=godemo()[["bpless"]],
             "Molecular function greater"=godemo()[["mfupper"]],
             "Molecular function less"=godemo()[["mfless"]],
             "Cellular component greater"=godemo()[["ccupper"]],
             "Cellular component less"=godemo()[["ccless"]],
             "KEGG greater"=godemo()[["kgupper"]],
             "KEGG component less"=godemo()[["kgless"]])
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
    
  })
  
  observeEvent(input$reset, {
    
    ctobj <- reactive({
      req(input$file1)
      count <- read_excel(input$file1$datapath)
      return(count)
    })
    
    pobj <- reactive({
      req(input$file2)
      pheno <- read_excel(input$file2$datapath)
      return(pheno=pheno)
    })
    
    
    ### matrix file
    output$matrix <- renderTable({
      head(ctobj(), 10)
    })
    
    output$pdat <- renderTable({
      head(pobj(),  10)
    })
    
    ### pheno file 
    observe({
      updateSelectInput(session, "design", choices=names(pobj()))
    })
    
    output$level0 <- renderUI({
      selectInput("ref0", "Reference group", pobj()[[input$design]])
    })
    
    output$level1 <- renderUI({
      selectInput("ref1", "Study group", pobj()[[input$design]])
    })
    

      resobj <- reactiveVal()
      volplot <- reactiveVal()
      gores <- reactiveVal()
      observeEvent(input$runbutton, {
        resobj <- (NULL)
        volplot <- (NULL)
        gores <- (NULL)
        withProgress(message = 'Running ...', style="old", value=0, {
          Sys.sleep(1)
          resobj(getdgeres(ctobj(), pobj(),
                        comparison=input$design, 
                        level1=input$ref1,
                        level0=input$ref0))
          incProgress(0.4, detail="Differeital expression ... ")
          
          Sys.sleep(2)
          gores(getgores(resobj()[["results"]], species=input$species))
          incProgress(0.4, detail="Pahtway analysis ... ")
          
          Sys.sleep(3)
          volplot(getvolcano(resobj()[["results"]]))
          incProgress(0.2, detail="making plot ... ")
        })
      })
      
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
      
      #tabPanel-Results
      output$plot1 <- renderPlot({
        volplot()[["plot"]]
      })
      
      output$brush_info <- renderPrint({
        showdf<-volplot()[["brush"]]%>%dplyr::select(Gene, log2FoldChange, pvalue, padj, log10padj) 
        brushedPoints(showdf, input$plot1_brush)
      })
      
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
  })
  
  readme <- a("README", href="https://github.com/DadongZ/RNASeqDGE/blob/master/README.md")
  output$README <- renderUI({
    tagList("Need help? ", readme)
  })
  
  issue <- a("issues", href="https://github.com/DadongZ/RNASeqDGE/issues")
  output$issue <- renderUI({
    tagList("Please report issues at: ", issue)
  })
  
})

