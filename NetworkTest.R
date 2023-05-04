library(igraph) # Load the igraph package
g1 <- graph( edges=c(1,2, 2,3, 3, 1), n=3, directed=F ) 

plot(g1) # A simple plot of the network - we'll talk more about plots later
g2 <- graph( edges=c(1,2, 2,3, 3, 1), n=10 )

plot(g2) 
g3 <- graph( c("John", "Jim", "Jim", "Jill", "Jill", "John")) # named vertices

# When the edge list has vertex names, the number of nodes is not needed

plot(g3)
g4 <- graph( c("John", "Jim", "Jim", "Jack", "Jim", "Jack", "John", "John"), 

             isolates=c("Jesse", "Janis", "Jennifer", "Justin") )  

# In named graphs we can specify isolates by providing a list of their names.



plot(g4, edge.arrow.size=.5, vertex.color="gold", vertex.size=15, 

     vertex.frame.color="gray", vertex.label.color="black", 

     vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.2) 
     
eg <- make_empty_graph(40)

plot(eg, vertex.size=10, vertex.label=NA)

nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)

links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)

library(igraph)



net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 

class(net)

E(net)       # The edges of the "net" object

V(net)       # The vertices of the "net" object

E(net)$type  # Edge attribute "type"

V(net)$media # Vertex attribute "media"

