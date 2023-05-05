#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load required packages
library(shiny)
require(tidyverse)
require(visNetwork)

# Source functions
source("CTP_helpers.R")

# Load compound paths targets full data
# CompoundTargetsPaths data
SJ_CTP_list <- readRDS("SJ_CTP_list.rds")
# PathsTargets data 
TP_list <- readRDS("PathsTargets_list.rds")
# Node types
node_types <- readRDS("node_types.rds")

# Define UI
ui <- fluidPage(
  titlePanel("Compound Targets Pathways App"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload file"),
      br(),
      textAreaInput("paste", "Or paste SJNs (compound IDs) here"),
      br(),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      h4("Pathway IDs hit by genes from selected compounds:"),
      verbatimTextOutput("output"),
      visNetworkOutput("mynetworkid")
    )
  )
)

# Define server
server <- function(input, output) {
  
  # Extract pathway IDs from file or paste input
  observeEvent(input$submit, {
    
    # Check if file is uploaded
    if (!is.null(input$file)) {
      data <- readLines(input$file$datapath)
    } else {
      data <- strsplit(input$paste, "\n")[[1]]
    }
    
    # Extract pathway IDs
    #pathway_ids <- data[grep("^PW\\d+", data, value = TRUE)]
    SJN_ids <- data
    
    curr_CTP_list <- SJ_CTP_list[SJN_ids]
    
    # Use data and get edgelist and nodes
    curr_CTP_list_edgelist<-data.frame()
    for(i in curr_CTP_list){
      curr_CTP_list_edgelist<-rbind(curr_CTP_list_edgelist,get_CompoundTargetsPaths_edgelist(i))
    }
    
    curr_CTP_list_nodes <- nodes_from_edge_list(curr_CTP_list_edgelist)
    curr_CTP_list_nodes <- dplyr::left_join(data.frame(node=curr_CTP_list_nodes),node_types)
    curr_CTP_list_nodes$color <- case_when(curr_CTP_list_nodes$type == 'gene' ~ '#8DD3C7', curr_CTP_list_nodes$type == 'pathway' ~ '#FFFFB3')
    
    curr_CTP_list_paths <- curr_CTP_list_nodes %>% 
      filter(curr_CTP_list_nodes$type == 'pathway') %>%
      pull('node')
    
    # Make network from data
    output$output <- renderPrint(curr_CTP_list_paths)
    
    output$mynetworkid <- renderVisNetwork({
      # minimal example
      nodes <- data.frame(id = curr_CTP_list_nodes$node,
                          label = curr_CTP_list_nodes$node,
                          group = curr_CTP_list_nodes$type,
                          title = curr_CTP_list_nodes$node,
                          color = curr_CTP_list_nodes$color)
      edges <- data.frame(from=curr_CTP_list_edgelist$edge1,to=curr_CTP_list_edgelist$edge2)
      
      visNetwork(nodes, edges) %>% 
        visGroups(groupname = "gene", color = "#8DD3C7") %>%
        visGroups(groupname = "pathway", color = "#FFFFB3") %>%
        visLegend(width = 0.1, position = "right", main = "group")
    })
    
  })
  
}

# Run the app
shinyApp(ui, server)

