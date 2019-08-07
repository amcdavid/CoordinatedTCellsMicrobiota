ui <- navbarPage("Is this the real life?",
                 tabPanel("Is this just fantasy?"),
                 tabPanel("Caught in a landslide"),
                 tabPanel("FCS Subsampling",
                          fluidPage(
                            inputPanel(
                              shinyDirButton(id = 'folder',
                                             label = "Folder select", 
                                             title = "Please select a folder", 
                                             FALSE), #dirGetter 'hidden' logical value
                              selectInput(inputId = 'unique.dirs',
                                          label = "Unique folders",
                                          choices = "update",
                                          multiple = TRUE)
                            )
                          ))
)

server <- function(input, output, session){
  
  roots = c(wd = '.')
  
  shinyDirChoose(input, 'folder', roots = roots, filetypes = c('', 'fcs'))
  
  observe({
    unique.paths <- basename(unique(dirname(list.files(path = parseDirPath(roots, input$folder), full.names = TRUE, recursive = TRUE, pattern = ".fcs"))))
    updateSelectInput(session, inputId = 'unique.dirs', label = "Unique folders", choices = unique.paths)
  })
  
  fcs.input <- reactive({list.files(path = input$unique.dirs, full.names = TRUE, pattern = ".fcs")})
  
  # FCS <- reactiveVal()
  # 
  # observeEvent(ignoreNULL = TRUE, eventExpr = {input$fcsDirButton},
  #              handlerExpr = {
  #                if (input$fcsDirButton > 0){
  #                  fp <- input$fcsDir
  #                }
  #                FCS(dir(fp, recursive=TRUE, full.names=TRUE, pattern=".fcs$"))
  #                #output$updates <- renderText({FCS()})
  #              })
  # 
  # observeEvent(ignoreNULL = TRUE, eventExpr = {input$subsampleButton},
  #              handlerExpr = {
  #                if (input$subsampleButton == 1){
  # 
  #                  FCS.paths <- FCS()
  #                  samp.size <- input$subsampleNum
  #                  newfolder.name <- input$subDirectory
  #                  post.name <- "subsampled.fcs"
  # 
  #                  FCS.subsampled(FCS.paths, samp.size, newfolder.name, post.name)
  # 
  #                }
  #              })
}#Close Server

shinyApp(ui = ui, server = server)



