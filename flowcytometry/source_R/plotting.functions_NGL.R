asinh.viz <- function(input.dat) {
  
  dat <- input.dat
  
  shinyApp(
    ui = fluidPage(
      
      titlePanel("asinh transform"),
      
      mainPanel(
        plotOutput("asinh.plot"),
        inputPanel(
          selectInput(
            inputId = 'x',
            label = "X",
            choices = names(dat),
            selected = NULL
          ),
          sliderInput(
            inputId = 'a',
            label = "a",
            min = -50,
            max = 50,
            value = 1,
            step = 1
          ),
          sliderInput(
            inputId = 'b',
            label = "b",
            min = -5,
            max = 5,
            value = 1,
            step = .02
          ),
          sliderInput(
            inputId = 'b.frac',
            label = "b fraction",
            min = 0,
            max = 2000,
            value = 1,
            step = 1
          ),
          sliderInput(
            inputId = 'c',
            label = "c",
            min = -10,
            max = 10,
            value = 1,
            step = 0.5
          )
        )
      )
    ),
    
    server <- function(input, output) {
      
      datasetInput <- reactive({
        tmp <- subset(dat, select = input$x)
      })
      
      output$asinh.plot <- renderPlot({
        tmp <- subset(dat, select = input$x)
        ggplot(asinh(input$a + (input$b/input$b.frac) * datasetInput()) + input$c, aes_string(x = names(datasetInput()))) + geom_density()
      })
    }
  )
}


fcs.explorer <- function(input.dat, name.X1, name.Y1, points) {
  
  dat <- input.dat
  colnames(dat)[grep("cluster", colnames(dat), ignore.case = T)] <- "MCluster"
  colnames(dat)[grep("node|nodes", colnames(dat), ignore.case = T)] <- "nodes"
  
  shinyApp(
    ui = fluidPage(
      
      titlePanel("FCS Explorer"),
      
      sidebarLayout(
        sidebarPanel(
          
          sliderInput('sampleSize', 'Sample Size', min=0, max=1000000,
                      value=points),
          
          selectInput('x', 'X.1', names(dat), name.X1),
          selectInput('y', 'Y.1', names(dat), name.Y1),
          
          selectInput('x2', 'X.2', names(dat), names(dat)[[1]]),
          selectInput('y2', 'Y.2', names(dat), names(dat)[[2]]),
          
          selectInput('color', 'Color', c('None', names(dat))),
          
          selectInput('facet_row', 'Facet Row', c(None='.', names(dat))),
          #selectInput('facet_col', 'Facet Column', c(None='.', names(dataset)))
          
          selectInput('sample', 'Sample', c('All', unique(dat$dat))),
          
          if(length(grep("node", colnames(dat), ignore.case = T)) == 1) {
            selectInput('node', 'Node', c('None', sort(unique(dat[, grep("node", colnames(dat), ignore.case = T)]))))
          } else {
            selectInput('node', 'Node', c('None'))
          },
          
          if(length(grep("cluster", colnames(dat), ignore.case = T)) == 1) {
            selectInput('mcluster', 'Meta-cluster', c('None', sort(unique(dat[, grep("cluster", colnames(dat), ignore.case = T)]))))
          } else {
            selectInput('mcluster', 'Meta-cluster', c('None'))
          },
          
          selectInput('cid', 'Cluster.ID', c('None', sort(unique(dat$cid))))
          
        ),
        
        mainPanel(
          plotOutput('plot'),
          plotOutput('plot2')
        )
      )),
    server = function(input, output) {
      
      dataset <- reactive({
        dat[sample(nrow(dat), input$sampleSize), ]
      })
      
      output$plot <- renderPlot({
        
        p <- ggplot(dataset(), aes_string(x=input$x, y=input$y)) + 
          geom_point(size = .5) +
          theme(axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18))
        
        if (input$color != 'None')
          p <- p + aes_string(color=input$color) + scale_color_viridis(option = "plasma")
        
        #facets <- paste(input$facet_row, '~', input$facet_col)
        facets <- input$facet_row
        #if (facets != '. ~ .')
        if (facets != '.')
          p <- p + facet_wrap(facets)
        
        if (input$node != 'None')
          p <- p + geom_point(data = subset(dataset(), nodes==input$node), shape = 21, size = 3, stroke = 1, color = "blue")
        
        if (input$mcluster != 'None')
          p <- p + geom_point(data = subset(dataset(), MCluster == input$mcluster), shape = 21, size = 3, stroke = 1, color = "red")
        
        if (input$cid != 'None')
          p <- p + geom_point(data = subset(dataset(), cid==input$cid), shape = 21, size = 3, stroke = 1, color = "red")
        
        if (input$sample != 'All')
          p <- p + geom_point(data = subset(dataset(), dat==input$sample), shape = 21, size = 3, stroke = 1, color = "green")
        
        print(p)
        
      })
      
      output$plot2 <- renderPlot({
        
        p <- ggplot(dataset(), aes_string(x=input$x2, y=input$y2)) + 
          geom_point(size = .5) +
          theme(axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18))
        
        if (input$color != 'None')
          p <- p + aes_string(color=input$color) + scale_color_viridis(option = "plasma")
        
        #facets <- paste(input$facet_row, '~', input$facet_col)
        facets <- input$facet_row
        #if (facets != '. ~ .')
        if (facets != '.')
          p <- p + facet_wrap(facets)
        
        if (input$node != 'None')
          p <- p + geom_point(data = subset(dataset(), nodes==input$node), shape = 21, size = 3, stroke = 1, color = "blue")
        
        if (input$mcluster != 'None')
          p <- p + geom_point(data = subset(dataset(), MCluster==input$mcluster), shape = 21, size = 3, stroke = 1, color = "red")
        
        if (input$cid != 'None')
          p <- p + geom_point(data = subset(dataset(), cid==input$cid), shape = 21, size = 3, stroke = 1, color = "red")
        
        print(p)
      })
      
    },
    options = list()
  )
}


meta.explorer <- function(input.dat, name.X1, name.Y1) {
  
  dat <- input.dat
  
  shinyApp(
    ui = fluidPage(
      
      titlePanel("Populations/Meta Explorer"),
      
      sidebarLayout(
        sidebarPanel(
          
          selectInput('x', 'X', names(dat), name.X1),
          selectInput('y', 'Y', names(dat), name.Y1),
          selectInput('color', 'Color', c('None', names(dat)), names(dat)[[2]]),
          
          checkboxInput('jitter', 'Jitter'),
          checkboxInput('smooth', 'Smooth'),
          checkboxInput('boxplot', 'Boxplot'),
          checkboxInput('violin', 'Violin'),
          
          selectInput('facet_row', 'Facet Row', c(None='.', names(dat))),
          selectInput('facet_col', 'Facet Column', c(None='.', names(dat)))
        ),
        
        mainPanel(
          plotOutput('plot')
        )
      )),
    
    server = function(input, output) {
      
      dataset <- reactive({
        dat
      })
      
      output$plot <- renderPlot({
        
        p <- ggplot(dataset(), aes_string(x = input$x, y = input$y)) + 
          geom_jitter(size = 2, width = 0.075, show.legend = TRUE) +
          #labs(title = "% of CD4+CD69+ (activated)") +
          theme(axis.title.x = element_text(size = 18),
                axis.text.x = element_text(size = 18),
                axis.title.y = element_text(size = 20, margin = margin(r = 20)),
                axis.text.y = element_text(size = 18),
                title = element_text(size = 20))
        
        if (input$color != 'None')
          p <- p + aes_string(color=input$color) +
            theme(#axis.title.x = element_blank(),
              #axis.text.x = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 18))
        
        facets <- paste(input$facet_row, '~', input$facet_col)
        if (facets != '. ~ .')
          p <- p + facet_grid(facets) + 
          theme(strip.text.x = element_text(size = 18),
                strip.text.y = element_text(size = 16))
        
        if (input$jitter)
          p <- p + geom_jitter()
        if (input$smooth)
          p <- p + geom_smooth()
        if(input$boxplot)
          p <- p + geom_boxplot()
        if(input$violin)
          p <- p + geom_violin()
        
        print(p)
        
      })
    })
}


fsom.clusters.heatmap <- function(fSOM.object){
  
  meta_MFI <- matrix(NA, nrow = length(levels(fSOM.object$metaclustering)),
                     ncol = length(fSOM.object$map$colsUsed),
                     dimnames = list(paste0("Meta.Cluster_", levels(fSOM.object$metaclustering)), 
                                     fSOM.object$markers[fSOM.object$map$colsUsed]))
  
  for(i in levels(fSOM.object$metaclustering)){
    #mfis <- apply(subset(input, input$fSOM.MetaCluster==i)[, c(1:30)], 2, FUN = median)
    if (length(grep(TRUE,fSOM.object$metaclustering==i)) > 1){
      mfis <- apply(fSOM.object$map$medianValues[fSOM.object$metaclustering==i,][,c(fSOM.object$map$colsUsed)], 2, FUN = median)
    }
    else {
      mfis <-  fSOM.object$map$medianValues[fSOM.object$metaclustering==i,][c(fSOM.object$map$colsUsed)]
    }
    #mfis <- apply(fSOM$data[fSOM$metaclustering[fSOM$map$mapping[,1]]==i,][,c(fSOM$map$colsUsed)],2,median)
    names(mfis) <- fSOM.object$markers[fSOM.object$map$colsUsed]
    meta_MFI[paste0("Meta.Cluster_",i), ] <- mfis
    #fSOM$meta_MFI <- meta_MFI
  } 
  return(meta_MFI)
}

fsom.dat <- function(fSOM.object) {
  dat <- as.data.frame(cbind(fSOM.object$data,
                             nodes = fSOM.object$map$mapping[, 1], 
                             MCluster = fSOM.object$metaclustering[fSOM.object$map$mapping[, 1]]))
  names(dat)[1:length(fSOM.object$markers)] <- fSOM.object$markers
  dat
}
