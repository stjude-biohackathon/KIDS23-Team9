# Load required packages
library(shiny)
library(shinyjs)
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
node_types_full <- readRDS("node_types_full.rds")
# Gene compound strengths
gene_compound_strengths <- readRDS("gene_compound_strengths.rds")

ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$img(src="logo.png", width="150px")
  ),
  
  # Tabset
  tabsetPanel(
    # First tab for compound IDs
    tabPanel("Compound ID Query",
             sidebarLayout(
               sidebarPanel(
                 fileInput("compound_id_file", "To view a network of all compounds, upload .txt file containing compound IDs (each line as a compound id):"),
                 br(),
                 textAreaInput("compound_id_text", "To view a network of all compounds or individual compounds, paste compound IDs here:", height = "200px"),
                 br(),
                 actionButton("compound_id_submit_btn", "Submit"),
                 actionButton("reset_compound", "Reset Compounds"),
                 #br(),
                 #uiOutput("compound_vector_list_ui")
                 verbatimTextOutput("output_tab1"),
                 downloadButton("download_pathway_list", "Download Pathway List"),
               ),
               
               mainPanel(
                 h4("Pathway IDs hit by genes from selected compounds:"),
                 #verbatimTextOutput("output_tab1"),
                 visNetworkOutput("mynetworkid_tab1",
                                  height = "1000px")
               )
             )),
    
    # Second tab for pathway IDs
    tabPanel("Pathway ID Query",
             sidebarLayout(
               sidebarPanel(
                 fileInput("pathway_id_file", "Upload .txt file containing pathway IDs (each line as a Pathway id):"),
                 br(),
                 textAreaInput("pathway_id_text", "Or paste pathway IDs here:", height = "200px"),
                 br(),
                 actionButton("pathway_id_submit_btn", "Submit"),
                 actionButton("reset_pathway", "Reset Pathways"),
                 #br(),
                 #uiOutput("pathway_vector_list_ui")
                 verbatimTextOutput("output_tab2"),
                 downloadButton("download_gene_list", "Download Gene Target List"),
               ),
               
               mainPanel(
                 h4("Compound IDs (SJNs) that target genes from the input pathways:"),
                 #verbatimTextOutput("output_tab2"),
                 visNetworkOutput("mynetworkid_tab2",
                                  height = "1000px")
               )
             ))
  )
)

server <- function(input, output) {
  #compound_id_vector <- reactiveVal(NULL)
  #pathway_id_vector <- reactiveVal(NULL)
  # adding reset button
  observeEvent(input$reset_compound, {
    # reset everything and set everything to NULL
    reset("compound_id_file")
    reset("compound_id_text")
    # set everything to NULL
    id_file <- NULL
    
    # add the pasted text to the input
    paste_text <- NULL
    id_vector_val <- NULL
    
    SJN_ids <- NULL
    #print(SJN_ids)
    curr_CTP_list <- NULL
    # Use data and get edgelist and nodes
    curr_CTP_list_edgelist<-NULL
    
    # Output the list of KEGG/HALLMARK pathways in a text-friendly format
    output$output_tab1 <- renderText({
      pathway_list <- NULL
      print("Current Compound list Reset")
    })
    
    # add in the download KEGG/HALLMARK pathways in a text-friendly format button logic
    output$download_pathway_list <- downloadHandler(
      filename = function () {
        paste("pathway_hits","txt",sep = "\n")
      },
      content = function(file) {
        pathway_list <- c(curr_CTP_list_paths)
        writeLines(paste(pathway_list,collapse = "\n"), file)
      }
    )
    
    #output$
    # Make network from data
    output$mynetworkid_tab1 <- renderVisNetwork({
      # minimal example
      print("Current Compound list Reset")
    })
    
  })
  
  observeEvent(input$reset_pathway, {
    # reset everything and set everything to NULL
    reset("pathway_id_file")
    reset("pathway_id_text")
    # set everything to NULL
    id_file <- NULL
    
    # add the pasted text to the input
    paste_text <- NULL
    id_vector_val <- NULL
    
    #SJN_ids <- NULL
    KEGG_paths <- NULL
    #print(SJN_ids)
    curr_CTP_list <- NULL
    # Use data and get edgelist and nodes
    curr_CTP_list_edgelist<-NULL
    
    # Output the list of KEGG/HALLMARK pathways in a text-friendly format
    output$output_tab2 <- renderText({
      gene_list <- NULL
      print("Current Pathway List Reset")
    })
    
    # add in the download KEGG/HALLMARK pathways in a text-friendly format button logic
    output$download_gene_list <- downloadHandler(
      filename = function () {
        paste("gene_hits","txt",sep = "\n")
      },
      content = function(file) {
        gene_list <- c(curr_TP_list_genes)
        writeLines(paste(gene_list,collapse = "\n"), file)
      }
    )
    
    #output$
    # Make network from data
    output$mynetworkid_tab2 <- renderVisNetwork({
      # minimal example
      print("Current Compound list Reset")
    })
    
  })
  
  observeEvent(input$compound_id_submit_btn, {
    # read in the uploaded file
    id_file <- input$compound_id_file
    if (!is.null(id_file)) {
      id_text <- readLines(id_file$datapath)
    } else {
      id_text <- ""
    }
    
    # add the pasted text to the input
    paste_text <- input$compound_id_text
    # split the pasted text by carriage return or newline
    paste_text <- unlist(strsplit(paste_text, "\r?\n"))
    # append the new pasted text to the ID text
    id_text <- append(id_text, paste_text)
    # remove any leading or trailing whitespace
    id_text <- trimws(id_text)
    # split the input string into a vector of IDs
    id_vector_val <- unlist(strsplit(id_text, "\r?\n"))
    
    # filter out any IDs that do not match the correct format
    id_vector_val <- id_vector_val[grep("^SJ\\d+", id_vector_val)]
    
    SJN_ids <- id_vector_val
    #print(SJN_ids)
    curr_CTP_list <- SJ_CTP_list[SJN_ids]
    
    # Use data and get edgelist and nodes
    curr_CTP_list_edgelist<-data.frame()
    for(i in curr_CTP_list){
      if(!is.null(i)) {
        curr_CTP_list_edgelist<-rbind(curr_CTP_list_edgelist,get_CompoundTargetsPaths_edgelist(i))
      }
    }
    
    curr_CTP_list_nodes <- nodes_from_edge_list(curr_CTP_list_edgelist)
    curr_CTP_list_nodes <- dplyr::left_join(data.frame(node=curr_CTP_list_nodes),node_types_full)
    curr_CTP_list_nodes$color <- case_when(curr_CTP_list_nodes$type == 'gene' ~ '#8DD3C7', curr_CTP_list_nodes$type == 'pathway' ~ '#FFFFB3')
    curr_CTP_list_nodes <- dplyr::left_join(curr_CTP_list_nodes,gene_compound_strengths,by=c("node"="gene"))
    
    curr_CTP_list_paths <- curr_CTP_list_nodes %>% 
      filter(curr_CTP_list_nodes$type == 'pathway') %>%
      pull('node')
    #curr_CTP_list_paths <- curr_CTP_list_edgelist
    
    # Output the list of KEGG/HALLMARK pathways in a text-friendly format
    output$output_tab1 <- renderText({
      pathway_list <- c(curr_CTP_list_paths)
      paste(pathway_list,collapse="\n")
    })
    
    # add in the download KEGG/HALLMARK pathways in a text-friendly format button logic
    output$download_pathway_list <- downloadHandler(
      filename = function () {
        paste("pathway_hits","txt",sep = "\n")
      },
      content = function(file) {
        pathway_list <- c(curr_CTP_list_paths)
        writeLines(paste(pathway_list,collapse = "\n"), file)
      }
    )
    
    #output$
    # Make network from data
    output$mynetworkid_tab1 <- renderVisNetwork({
      # minimal example
      if(!is.null(curr_CTP_list_nodes)) {
        nodes <- data.frame(id = curr_CTP_list_nodes$node,
                            label = curr_CTP_list_nodes$node,
                            group = curr_CTP_list_nodes$type,
                            title = paste0("<p><b>", curr_CTP_list_nodes$node,"</b><br>Compounds and Active Strengths<br>", curr_CTP_list_nodes$compound_strengths,"</p>"),
                            color = curr_CTP_list_nodes$color)
        edges <- data.frame(from=curr_CTP_list_edgelist$edge1,to=curr_CTP_list_edgelist$edge2)
        
        visNetwork(nodes, edges) %>% 
          visGroups(groupname = "gene", color = "#8DD3C7") %>%
          visGroups(groupname = "pathway", color = "#FFFFB3") %>%
          visLegend(width = 0.1, position = "right", main = "group")
      }
    })
    
  })
  
  observeEvent(input$pathway_id_submit_btn, {
    # read in the uploaded file
    id_file <- input$pathway_id_file
    if (!is.null(id_file)) {
      id_text <- readLines(id_file$datapath)
    } else {
      id_text <- ""
    }
    
    # add the pasted text to the input
    paste_text <- input$pathway_id_text
    # remove carriage return or newline and make new vector
    paste_text <- unlist(strsplit(paste_text, "\r?\n"))
    # append the id_text with the new paste text
    id_text <- append(id_text, paste_text)
    
    # remove any leading or trailing whitespace
    id_text <- trimws(id_text)
    #print(id_text)
    # split the input string into a vector of IDs
    id_vector_val <- unlist(strsplit(id_text, "\r?\n"))
    
    # set the KEGG paths to our id vector values
    KEGG_paths <- id_vector_val
    curr_TP_list <- TP_list[KEGG_paths]
    
    # Use data and get edgelist and nodes
    curr_TP_list_edgelist<-data.frame()
    for(i in curr_TP_list){
      if(!is.null(i)) {
        curr_TP_list_edgelist<-rbind(curr_TP_list_edgelist,get_PathsTargets_edgelist(i))
      }
    }
    
    curr_TP_list_nodes <- nodes_from_edge_list(curr_TP_list_edgelist)
    curr_TP_list_nodes <- dplyr::left_join(data.frame(node=curr_TP_list_nodes),node_types_full)
    curr_TP_list_nodes$color <- case_when(curr_TP_list_nodes$type == 'gene' ~ '#8DD3C7', curr_TP_list_nodes$type == 'pathway' ~ '#FFFFB3', curr_TP_list_nodes$type == 'gene_no_compound' ~ '#A1CAF1')
    curr_TP_list_nodes <- dplyr::left_join(curr_TP_list_nodes,gene_compound_strengths,by=c("node"="gene"))
    
    curr_TP_list_genes <- curr_TP_list_nodes %>% 
      filter(curr_TP_list_nodes$type == 'gene') %>%
      pull('node')
    
    # Output the list of gene hits/symbols in a text-friendly format
    output$output_tab2 <- renderText({
      gene_list <- c(curr_TP_list_genes)
      paste(gene_list,collapse="\n")
    })
    
    # add in the download KEGG/HALLMARK pathways in a text-friendly format button logic
    output$download_gene_list <- downloadHandler(
      filename = function () {
        paste("gene_hits","txt",sep = "\n")
      },
      content = function(file) {
        gene_list <- c(curr_TP_list_genes)
        writeLines(paste(gene_list,collapse = "\n"), file)
      }
    )
    
    # Make network from data
    #output$output_tab2 <- renderPrint(curr_TP_list_genes)
    
    # Make network from data
    
    output$mynetworkid_tab2 <- renderVisNetwork({
      # minimal example
      if(!is.null(curr_TP_list_nodes)) {
        nodes <- data.frame(id = curr_TP_list_nodes$node,
                            label = curr_TP_list_nodes$node,
                            group = curr_TP_list_nodes$type,
                            title = paste0("<p><b>", curr_TP_list_nodes$node,"</b><br>Compounds and Active Strengths<br>", curr_TP_list_nodes$compound_strengths,"</p>"),
                            color = curr_TP_list_nodes$color)
        edges <- data.frame(from=curr_TP_list_edgelist$edge1,to=curr_TP_list_edgelist$edge2)
        
        visNetwork(nodes, edges) %>% 
          visGroups(groupname = "gene", color = "#8DD3C7") %>%
          visGroups(groupname = "pathway", color = "#FFFFB3") %>%
          visGroups(groupname = 'gene_no_compound', color = '#A1CAF1') %>%
          visLegend(width = 0.1, position = "right", main = "group")
      }
    })
    
  })
}

shinyApp(ui, server)