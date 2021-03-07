
##======================================================##
##                                                      ##
##    POLNET 2015 Network Visualization Workshop        ##
##    Porland OR, June 18, 2015                         ##
##                                                      ##
##    Katya Ognyanova, katya@ognyanova.net              ##
##    www.kateto.net/polnet2015                         ##
##                                                      ##
##======================================================##


# Download handouts and example data from www.kateto.net/polnet2015


# Key packages to install: 

install.packages("igraph") 
install.packages("network") 
install.packages("ndtv")
library(igraph)
library(network)
library(ndtv)

# Set the working directory to the folder containing the example data
setwd("C:/Data") 


# ================ Read the example data ================


# DATASET 1: edgelist 

nodes <- read.csv("c:\\Users\\user\\Desktop\\新增資料夾\\Data\\Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("c:\\Users\\user\\Desktop\\新增資料夾\\Data\\Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)

# Examine the data:
head(nodes)
head(links)
nrow(nodes); length(unique(nodes$id))
nrow(links); nrow(unique(links[,c("from", "to")]))

# Collapse multiple links of the same type between the same two nodes
# by summing their weights, using aggregate() by "from", "to", & "type":

links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL


# DATASET 2: matrix 

#nodes2 <- read.csv("c:\\Users\\user\\Desktop\\新增資料夾\\Data\\Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
#links2 <- read.csv("c:\\Users\\user\\Desktop\\新增資料夾\\Data\\Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)

# Examine the data:
#head(nodes2)
#head(links2)

# links2 is an adjacency matrix for a two-mode network:
#links2 <- as.matrix(links2)
#dim(links2)
#dim(nodes2)

# ================ Plotting networks with igraph ================
 
library(igraph)

# Converting the data to an igraph object:
net <- graph.data.frame(links, nodes, directed=T) 

# Examine the resulting object:
class(net)
net 

# It's easy to access nodes, edges, and their attributes:
E(net)
V(net)
E(net)$type
V(net)$media

# You can also manipulate the network matrix:
net[1,]
net[5,7]

# First attempt to plot the graph:
plot(net) # not pretty!

# Removing loops from the graph:
net <- simplify(net, remove.multiple = F, remove.loops = T) #刪除自己跟自己

# Let's and reduce the arrow size and remove the labels:
plot(net, edge.arrow.size=.4,vertex.label=NA)



# ================ Back to network plots with igraph ================


#  ------->> Plot parameters in igraph --------

# Plotting with igraph: node options (starting with 'vertex.) and edge options
# (starting with 'edge.'). A list of options is included in your handout.
?igraph.plotting

# We can set the node & edge options in two ways - one is to specify
# them in the plot() function, as we are doing below.

# Plot with curved edges (edge.curved=.1) and reduce arrow size:
plot(net, edge.arrow.size=.4, edge.curved=.1)

# Set node color to orange and the border color to hex #555555
# Replace the vertex label with the node names stored in "media"
plot(net, edge.arrow.size=.4, edge.curved=0,
     vertex.color="orange", vertex.frame.color="#555555",
     vertex.label=V(net)$media, vertex.label.color="black",
     vertex.label.cex=.7) 


# The second way to set attributes is to add them to the igraph object.

# Generate colors base on media type:
colrs <- c("gray50", "tomato", "gold")
V(net)$color <- colrs[V(net)$media.type]

# Compute node degree (#links) and use it to set node size:

deg <- degree(net, mode="all")
V(net)$size <- deg*3
V(net)$size <- V(net)$audience.size*0.6

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label.color <- "black"
V(net)$label <- NA

# Set edge width based on weight:
E(net)$width <- E(net)$weight/6

#change arrow size and edge color:
E(net)$arrow.size <- .2
E(net)$edge.color <- "gray80"

plot(net) 


# We can also add a legend explaining the meaning of the colors we used:
plot(net) 
legend(x=-1.1, y=-1.1, c("Newspaper","Television", "Online News"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2.5, bty="n", ncol=1)


#






# ------->> Highlighting specific nodes or links --------


# Sometimes we want to focus the visualization on a particular node
# or a group of nodes. Let's represent distance from the NYT:
shortest.paths(net, algorithm="unweighted")
dist.from.NYT <- shortest.paths(net, algorithm="unweighted")[1,]

oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.NYT)+1)
col <- col[dist.from.NYT+1]

plot(net, vertex.color=col, vertex.label=dist.from.NYT, edge.arrow.size=.6, 
     vertex.label.color="white")


# Or we can show all the immediate neighbors of the WSJ:
col <- rep("grey40", vcount(net))
col[V(net)$media=="Wall Street Journal"] <- "#ff5100"

# The neighbors function finds all nodes one step out from the focal actor:
# (the corresponding function that finds all edges for a node is "incident")
neigh.nodes <- neighbors(net, V(net)[media=="Wall Street Journal"], mode="out")

col[neigh.nodes] <- "#ff9d00"
plot(net, vertex.color=col)


