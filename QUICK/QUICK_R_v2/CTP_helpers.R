### S4 CLASSES ###
setClass("CompoundTargetsPaths", slots = c(
  compound_description = "list",
  genes = "list",
  pathways="list",
  chemical_id = "list",
  MOA_info = "list"
))

# Make S4 class for PathsTargets
setClass("PathsTargets",representation(Pathway="character",
                                       Description="character",
                                       targets="list"
))


### FUNCTIONS ###
# Get edgelist from CompoundTargetsPaths object
get_CompoundTargetsPaths_edgelist<-function(CTP_object){
  CTP_SJN<-CTP_object@compound_description$SJN
  gene_paths<-CTP_object@pathways$pathways[[1]]
  edge_list<-data.frame()
  counter<-0
  for(i in gene_paths){
    counter<-counter+1
    if(length(i)>0){
      gene_name<-names(gene_paths)[counter]
      paths<-unlist(i)
      curr_edges<-data.frame(edge1=rep(gene_name,length(paths)),edge2=paths)
      edge_list<-rbind(edge_list,curr_edges)
    }
  }
  return(edge_list)
}

# Get PathsTargets edgelist
get_PathsTargets_edgelist<-function(TP_object){
  pathway<-TP_object@Pathway
  genes<-TP_object@targets[[1]]
  n_genes<-length(genes)
  edge_list<-data.frame(edge1=rep(pathway,n_genes),edge2=genes)
  return(edge_list)  
}

# Get nodes from edgelist
nodes_from_edge_list<-function(edgelist){
  nodes<-unique(c(edgelist$edge1,edgelist$edge2))
  return(nodes)
}

# Get gene-compound edges, from either CTP or TP object
get_gene_compound_edges<-function(S4_object){
  CTP_SJN<-CTP_object@compound_description$SJN
  genes<-CTP_object@genes$gene_symbols
  edge_list<-data.frame(edge1=genes,edge2=CTP_SJN)
  return(edge_list)
}