#
library(shiny)
library(ggplot2)
library(viridis)
library(plotly)
library(dplyr)

####

# FCS input data with mapped fSOM nodes/meta-clusters and tSNE X,Y embeddings
# heatmap annotation for meta-cluster, click to display population/cohort metrics

####

# is the app faster if you first load all datasets?

TPHE.CD4pos <- readRDS("./results_R/fSOMs/TPHE/TPHE.CD4pos.unified.RDS")

# edit input based on some QC; Mcluster 19 is 'junk'
TPHE.CD4pos$input.data <- TPHE.CD4pos$input.data[-grep("19", TPHE.CD4pos$input.data$MCluster), ]
TPHE.CD4pos$fsom.mclust.agg <- TPHE.CD4pos$fsom.mclust.agg[-grep("19", TPHE.CD4pos$fsom.mclust.agg$MCluster), ]

# UI
ui <- navbarPage("Displays:",
                 tabPanel("tSNE + FlowSOM",
                          fluidPage(
                            
                            fluidRow(
                              column(width = 12,
                                     align = "center",
                                     style = "background-color: lightgrey;",
                                     titlePanel(
                                       h3("tSNE/fSOM Explorer")),
                                     h5("Neonatal PBMC T-cell phenotyping;"),
                                     h5("Phenotype of individual cells measured using flow cytometry;"),
                                     h5("Dimensionality reduction using tSNE and objective clustering using FlowSOM allows detection of distinct T-cell populations")
                              )
                            ),
                            
                            fluidRow(
                              column(width = 4,
                                     inputPanel(
                                       selectInput(inputId = 'tsne_fsom.set',
                                                   label = "Choose a tSNE/fSOM dataset:", 
                                                   choices = c("CD4+" = "cd4"),
                                                   selected = "CD4+"
                                       )
                                     )
                              )
                            ),
                            
                            fluidRow(
                              column(width = 6,
                                     verticalLayout(
                                       plotOutput("tsnexy", 
                                                  click = "plot_click", 
                                                  hover = hoverOpts(id = "plot_hover", delay = 100)
                                       ),
                                       plotlyOutput("mcheatmap"),
                                       radioButtons(inputId = 'scale', 
                                                    label = "Scale heatmap:",
                                                    choices = c("Yes" = "yes", "No" = "no"),
                                                    inline = TRUE
                                       )
                                     )
                              ),
                              column(width = 6,
                                     tabsetPanel(
                                       tabPanel("Population/Sample Data",
                                                plotOutput("mcplot_jitter"),
                                                plotOutput("mcplot_line")
                                       ),
                                       tabPanel("FCS Data",
                                                plotOutput("fcsplot1"),
                                                fluidRow(
                                                  column(width = 6,
                                                         selectInput('x1', "X", "update")
                                                  ),
                                                  column(width = 6,
                                                         selectInput('y1', "Y", "update")
                                                  )
                                                ),
                                                plotOutput("fcsplot2"),
                                                fluidRow(
                                                  column(width = 6,
                                                         selectInput('x2', "X", "update")
                                                  ),
                                                  column(width = 6,
                                                         selectInput('y2', "Y", "update")
                                                  )
                                                )
                                       )
                                     )
                              )
                            )
                          )
                          ),
                 tabPanel("tSNE with marker overlays",
                          fluidPage(
                            fluidRow(
                              column(width = 6,
                                     verticalLayout(
                                       plotOutput("tsne_color.by.marker"),
                                       selectInput('color', "Color", "update")
                                     )
                              )
                            )
                          )
                 )
)

# server logic
server <- function(input, output, session) {
  
  observe({
    fcs.colnames <- colnames(datasetInput()$input.data)
    updateSelectInput(session, inputId = 'x1', label = "X", choices = fcs.colnames, selected = "CD31")
    updateSelectInput(session, inputId = 'y1', label = "Y", choices = fcs.colnames, selected = "CD197")
    updateSelectInput(session, inputId = 'x2', label = "X", choices = fcs.colnames, selected = "CD127")
    updateSelectInput(session, inputId = 'y2', label = "Y", choices = fcs.colnames, selected = "FOXP3")
    updateSelectInput(session, inputId = 'color', label = "Color", choices = fcs.colnames, selected = "CD31")
  })
  
  datasetInput <- reactive({
    switch(input$tsne_fsom.set,
           "cd4" = TPHE.CD4pos
    )
  })
  
  axis.yInput <- reactive({
    switch(input$tsne_fsom.set,
           "cd4" = "CD4+ T-cells"
    )
  })
  
  scaleInput <- reactive({
    switch(input$scale,
           "yes" = 1,
           "no"  = 0)
  })
  
  # generate tSNE map; Nav Panel 1
  output$tsnexy <- renderPlot({
    p <- ggplot(sample_n(datasetInput()$input.data, 50000), aes(x = tsne.X, y = tsne.Y, color = factor(MCluster))) +
      labs(title = "Click on colored meta-clusters") +
      geom_point(size = .75, show.legend = FALSE) +
      theme_void() +
      theme(panel.background = element_rect(fill = "grey95"),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            title = element_text(size = 18)) +
      geom_text(data = datasetInput()$fsom.mclust.agg, inherit.aes = FALSE, 
                aes(x = tsne.X, y = tsne.Y, label = rownames(datasetInput()$fsom.mclust.agg)), size = 6, nudge_x = -1)
    p
  })
  
  # generate tSNE map; Nav Panel 2
  output$tsne_color.by.marker <- renderPlot({
    p <- ggplot(sample_n(datasetInput()$input.data, 50000), aes(x = tsne.X, y = tsne.Y)) +
      geom_point(size = .75, show.legend = FALSE) +
      theme_void() +
      aes_string(color = input$color) + scale_color_viridis()
    p
  })
  
  observeEvent(input$plot_hover, {
    
    hover.val <- unique(nearPoints(datasetInput()$input.data, input$plot_hover, threshold = 5, maxpoints = 1)$MCluster)
    
    output$mcheatmap <- renderPlotly({
      plot_ly(x = datasetInput()$dims.used, 
              y = paste0("Meta.Cluster_", rownames(datasetInput()$fsom.mclust.agg)),
              z = if(scaleInput() == 0) {
                as.matrix(datasetInput()$fsom.mclust.agg)[, datasetInput()$dims.used]
              } else {
                scale(as.matrix(datasetInput()$fsom.mclust.agg)[, datasetInput()$dims.used])
              },
              type = "heatmap") %>%
        #layout(yaxis=list(autorange=F, range=c((hover.val-1.5),(hover.val-0.5))))
        layout(shapes =
                 list(type = "rect",
                      line = list(color = "purple", width = 5),
                      x0 = -0.5, x1 = length(datasetInput()$dims.used) - 0.5,
                      y0 = grep(paste0("^", hover.val, "$"), rownames(datasetInput()$fsom.mclust.agg)) - 1.5, 
                      y1 = grep(paste0("^", hover.val, "$"), rownames(datasetInput()$fsom.mclust.agg)) - 0.5)
        )
    })
  })
  
  # click events
  # click tSNE 'structures' (colored by fSOM meta-cluster) to determine metacluster value; pass to "MCplot" for plotting sample points
  observeEvent(input$plot_click, {
    #use nearPoints to set nearest meta-cluster value
    mc.val <- unique(nearPoints(datasetInput()$input.data, input$plot_click, threshold = 5, maxpoints = 1)$MCluster)
    tmp <- sample_n(datasetInput()$input.data, 50000)
    #pass meta-cluster value to second plot
    output$mcplot_jitter <- renderPlot({
      ggplot(datasetInput()$counts.meta, aes_string(x = "TermALIAS", y = paste0("Meta.Cluster_", mc.val))) + 
        geom_jitter(aes(color = TermALIAS), width = 0.2, show.legend = FALSE) + 
        facet_wrap(~VisitALIAS) +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
              strip.text = element_text(size = 20),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 15, b = 0, l = 0)),
              title = element_text(size = 18)) +
        labs(#title = "Grouped by blood draw",
          x = "Pre-term VS Full-term",
          y = paste0("% of ", paste0("Meta.Cluster ", mc.val, " "), "in ", axis.yInput()))
    })
    output$mcplot_line <- renderPlot({
      ggplot(datasetInput()$counts.meta, aes_string(x = "GAatBirth", y = paste0("Meta.Cluster_", mc.val))) + 
        geom_point(aes(color = factor(TermALIAS)), show.legend = FALSE) + 
        geom_smooth(method = loess) + 
        facet_wrap(~VisitALIAS) +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
              strip.text = element_text(size = 20),
              axis.text.x = element_text(size = 16),
              axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 10, l = 0)),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 15, b = 0, l = 0)),
              title = element_text(size = 18)) +
        labs(#title = "Grouped by blood draw",
          x = "Gestational Age at Birth (in weeks)",
          y = paste0("% of ", paste0("Meta.Cluster ", mc.val, " "), "in ", axis.yInput()))
    })
    output$fcsplot1 <- renderPlot({
      p <- ggplot(tmp, aes_string(x=input$x1, y=input$y1)) +
        geom_point(size = .5) +
        geom_point(data = subset(tmp, MCluster == mc.val), shape = 21, size = 3, stroke = 1, alpha = 0.4, color = "red")
        # theme(axis.title.x = element_text(size = 18),
        #       axis.title.y = element_text(size = 18))
      p

    })
    output$fcsplot2 <- renderPlot({
      p <- ggplot(tmp, aes_string(x=input$x2, y=input$y2)) +
        geom_point(size = .5) +
        geom_point(data = subset(tmp, MCluster == mc.val), shape = 21, size = 3, stroke = 1, alpha = 0.4, color = "red")
      # theme(axis.title.x = element_text(size = 18),
      #       axis.title.y = element_text(size = 18))
      p
      
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)