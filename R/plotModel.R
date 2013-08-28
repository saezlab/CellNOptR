#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$
plotModel <- function(model, CNOlist=NULL, bString=NULL, indexIntegr=NULL,
    signals=NULL, stimuli=NULL, inhibitors=NULL, NCNO=NULL, compressed=NULL,
    output="STDOUT", filename=NULL,graphvizParams=list(), show=TRUE, remove_dot=TRUE){
# Quick example:
# ---------------
#   filename = "ToyPKNMMB.sif"
#   g = plotModel(model, cnolist=cnolist)
#   # g$graph contains the model transformed into a graph object

  # user parameters to refine the layout, color, ...
  if (is.null(graphvizParams$arrowsize)==TRUE) {
    graphvizParams$arrowsize=2
    }
  if (is.null(graphvizParams$size)==TRUE) {
    graphvizParams$size="15,15"
  }
  if (is.null(graphvizParams$fontsize)==TRUE) {
    graphvizParams$fontsize=22
  }
  if (is.null(graphvizParams$edgecolor)==TRUE) {
    graphvizParams$edgecolor="black"
  }
  if (is.null(graphvizParams$nodeLabels)==TRUE) {
    graphvizParams$nodeLabels=NULL
  }
  if (is.null(graphvizParams$nodeWidth)==TRUE) {
    graphvizParams$nodeWidth=2
  }
  if (is.null(graphvizParams$nodeHeight)==TRUE) {
    graphvizParams$nodeHeight=1
  }
  if (is.null(graphvizParams$viewEmptyEdges)==TRUE) {
    graphvizParams$viewEmptyEdges = TRUE
  }
  # mode must be before andWidth and andHeight
  if (is.null(graphvizParams$mode)==TRUE) {
    graphvizParams$mode = "sbgn"
  }

  if (graphvizParams$mode %in% c("sbgn", "classic") == FALSE){
      stop("mode must be in 'classic' or 'sbgn'")
    }

  if (is.null(graphvizParams$andWidth)==TRUE) {
    if (graphvizParams$mode == "classic"){
        graphvizParams$andWidth = 0.2
    }
    if (graphvizParams$mode == "sbgn"){
        graphvizParams$andWidth = 0.5
    }

  }
  if (is.null(graphvizParams$andHeight)==TRUE) {
    if (graphvizParams$mode == "classic"){
        graphvizParams$andHeight = 0.2
    }
    if (graphvizParams$mode == "sbgn"){
        graphvizParams$andHeight = 0.5
    }
  }

    # Some required library to build the graph and plot the results using
    # graphviz.
    library(Rgraphviz)
    library(RBGL)

    # Set the output filename
    if (is.null(filename)){
        filename="model"
        output_dot = "model.dot"
    }
    else{
        output_dot = paste(filename, ".dot", sep="")
    }

    # Set the image format if any
    if (output %in% c("STDOUT", "PDF", "PNG", "SVG") != TRUE){
        stop("wrong output format.Must be in PDF, PNG, SVG")
    }
    else{
       if (output=="PDF"){
            pdf(paste(filename, ".pdf", sep=""))
        }
       else if (output=="PNG"){
            png(paste(filename, ".png", sep=""))
        }
       else if (output=="SVG"){
            svg(paste(filename, ".svg", sep=""), pointsize=22, width=10, height=10)
        }
    }

    # If the cnolist is a NULL, nothing to do
    if (is.null(CNOlist)==TRUE){
        cnolist = NULL
    } else { # if not a CNOlist, try to convert it 
        if ((class(CNOlist)=="CNOlist")==FALSE){
            cnolist = CellNOptR::CNOlist(CNOlist)
        } else{
            cnolist = CNOlist
        }
    }
   
    # TODO: The following piece of code should be made modular
 

    # Input data. If the model is a character, we guess that the user provided
    # the MODEL filename (sif format) that we can read directly.
    if (typeof(model) == "character"){
        raw = read.table(model)  # read the PKN data
        # build the unique vertices from the column 1 and 3 of the SIF file
        vertices = unique(c(as.character(raw$V1), as.character(raw$V3)))
        v1 = raw$V1 # alias to vertices in column 1
        v2 = raw$V3 # alias to vertices in column 2
        edges = raw$V2 # alias to the edges
        if (is.null(bString)){# default is only 1 so all edges are accepted
            BStimes<-rep(1,length(edges))
        }else{
            BStimes<-bString
        }
        Integr<-rep(0,length(edges))  #
    }
    # otherwise, the user probably provided the model already read by readSIF
    # in which case, and gates must be extracted from strings such as
    # "node1+node2=node3"
    # This block is the tricky part of the function. Change with care.
    else if (typeof(model)=="list" && any("namesSpecies" == names(model))){
        if (length(bString)!=length(model$reacID) & is.null(bString)==FALSE){
            stop(paste("If the bString argument is provided it must have the same",
            "  length as model$reacID. ", "The model has ", length(model$reacID),
            "  edges whereas the bitstring has a length of ", length(bString), sep=""))
        }
        # namesSpecies == names(model) try to check if model resemble the output of readSIF ?
        # ideally we should have a type.

        # build the unique vertices from the nameSpecies
        vertices = model$namesSpecies

        # now, we need to split the reaction to get back the different edges
        mysplit = function(x){strsplit(x, "=")}
        reacs = mysplit(model$reacID) # separate the reactions into left and right parts
        tmp <- unlist(mysplit(model$reacID))
        reacs = t(matrix(unlist(mysplit(model$reacID)), ncol=length(tmp)/2)) # reordering

        # Use the bString and indexIntegr input arguments to build up
        if (is.null(bString)){# default is only 1 so all edges are accepted
            optimBStimes<-rep(1,dim(reacs)[1])
        }else{
            optimBStimes<-bString
        }

        optIntegr<-rep(0,length(optimBStimes))
        if (is.null(indexIntegr)==FALSE){
            optIntegr[indexIntegr]<-1
        }

        # finally, build the v1, v2 and edges
        BStimes<-vector()
        Integr<-vector()

        CountReac<-1
        CountAnds<-1
        mysplit2 = function(x){strsplit(x, "+", fixed=TRUE)}
        v1<-vector()
        v2<-vector()
        edges<-vector()
        for (i in 1:dim(reacs)[1]){
          inputs<-unlist(strsplit(reacs[i,1],"+", fixed=TRUE))
          if (length(inputs)==1){
            v1[CountReac] = reacs[i,1]
            edges[CountReac] = 1
            v2[CountReac] = reacs[i,2]
            if (length(grep("!", v1))){
                v1[CountReac] = sub("!", "", v1[CountReac])
                edges[CountReac] = -1
            }
            BStimes[CountReac]<-optimBStimes[i]
            Integr[CountReac]<-optIntegr[i]
            CountReac<-CountReac+1
          }else{
            for (j in seq_along(inputs)){
              v1[CountReac]<-inputs[j]
              edges[CountReac] = 1
              v2[CountReac]<-paste("and",CountAnds,sep="")
              if (length(grep("!", v1[CountReac]))){
                v1[CountReac] = sub("!", "", v1[CountReac])
                edges[CountReac] = -1
              }
              BStimes[CountReac]<-optimBStimes[i]
              Integr[CountReac]<-optIntegr[i]
              CountReac<-CountReac+1
            }
            v1[CountReac]<-paste("and",CountAnds,sep="")
            edges[CountReac] = 1
            v2[CountReac] = reacs[i,2]
            BStimes[CountReac]<-optimBStimes[i]
            Integr[CountReac]<-optIntegr[i]
            CountReac<-CountReac+1
            vertices<-c(vertices,paste("and",CountAnds,sep=""))
            CountAnds<-CountAnds+1
          }
        }
    }

    if (is.null(cnolist) == FALSE){ # if a cnolist is provided, fill
                                    # signals/stimulis/inhitors
        stimuli <- colnames(cnolist@stimuli)
        signals <- colnames(cnolist@signals[[1]])
        inhibitors <- colnames(cnolist@inhibitors)
    }


    # check that the signal and stimuli are indeed in the list of vertices
    # otherwise they will be failures later.
    if (all(signals %in% vertices)==FALSE){
        msg = signals[signals %in% vertices == FALSE]
        print("Those signals were not found in the vertices: ")
        print(msg)
    }

    # build the edges. IGraph does not use names for the vertices but ids
    # that starts at zero. Let us build a data.frame to store the correspondence
    # between the ids and names.
    l = length(vertices) - 1

    # build the graph
    g <- new("graphNEL", nodes=vertices, edgemode="directed")
    weights = rep(1, l) # for now, the weights are to 1.
    for (i in 1:length(v1)){
        g <- addEdge(as.character(v1[i]), as.character(v2[i]), g, weights[i])
    }

    # The graph is now built. We can proceed with node and edge attributes for
    # the output files

    recipEdges="distinct" # an edge A->B does not overlap with B->A

    # --------------------------------------- Build the node and edges attributes list
    nodeAttrs = createNodeAttrs(g, vertices, stimuli, signals, inhibitors, NCNO,
        compressed, graphvizParams)

    res = createEdgeAttrs(v1, v2, edges, BStimes, Integr,
        user_edgecolor=graphvizParams$edgecolor,
        view_empty_edge=graphvizParams$viewEmptyEdges)
    # an alias
    edgeAttrs = res$edgeAttrs

    # createEdge returns the edgeAttrs and a list of edges that correspond to a
    # bistring element that is zero. In such case, the edge is useless and can
    # be removed.
    toremove = res$toremove
    for (x in toremove){
        y = unlist(strsplit(x, "~"))
        g = removeEdge(y[1], y[2], g)
    }
    # Some nodes are now connected to no other nodes. These nodes can be
    # removed. In principle, this is only "and" nodes.
    orphans = nodes(g)[(degree(g)$inDegree  + degree(g)$outDegree) ==0]
    for (x in orphans){
        if (x %in% stimuli == FALSE & x %in% inhibitors == FALSE & x %in% signals == FALSE){
            g = removeNode(x, g)
            # What was returned by createEdgeAttrs in now obsolet and need some cleanup
            f<-function(y){return(x %in% strsplit(y, "~", fixed=TRUE)[[1]])}
            indices = lapply(names(edgeAttrs$color),  f) == TRUE
            edgeAttrs$color[indices] <-NULL
            edgeAttrs$label[indices] <-NULL
            edgeAttrs$lty[indices] <-NULL
            edgeAttrs$arrowhead[indices] <-NULL
            edgeAttrs$penwidth[indices] <-NULL
        }
    }
    # we must rebuild the edges attributes

    # --------------------------- the ranks computation for the layout
    clusters = create_layout(g, signals, stimuli)

    # ------------------------------ general attributes
    # for some obscure reasons, must set style="filled" right her even though it
    # is then overwritten by nodesAttrs$style later on otherwise the output.dot
    # does not contain any style option

    # size does not seem to work in Rgraphviz version 1.32 wait and see for new version.
    fontsize=graphvizParams$fontsize
    attrs <- list(
        node=list(fontsize=fontsize,fontname="Helvetica",style="filled,bold"),
        edge=list(style="solid",penwidth=1,weight="1.0",arrowsize=graphvizParams$arrowsize,minlen=3),
        graph=list(splines=TRUE,size=graphvizParams$size,bgcolor="white",ratio="fill",pad="0.5,0.5",dpi=72)
        )

    copyg <- g

    # current version of Rgraphviz (1.32 feb 2012) does not handle edgewidth
    # properly. A
    if (installed.packages()[,"Version"]["Rgraphviz"] <= "1.33.0"){
        #nodeAttrs$lty = "solid"
        print("plotModel: please upgrade to Rgraphviz >1.33.0 for best output")
        nodelty = "solid"
        savedEdgeAttrs = edgeAttrs$color
        edgeAttrs$color = NULL
    }
    else{
        nodelty=nodeAttrs$lty
    }

    # Set the node Rendering in Rgraphviz
    nodeRenderAttrs  <- setNodeRenderInfo(nodeAttrs, list(lwd=2, lty=nodelty, 
       cex=0.4, fontsize=fontsize, fixedsize=FALSE))
    nodeRenderInfo(g) <- nodeRenderAttrs

   # the arrowhead "normal" is buggy in Rgraphviz version 1.32 so switch to
   # "open" for now. However, the dot output keeps using the normal arrow.
   arrowhead2 = edgeAttrs$arrowhead
   #arrowhead2[arrowhead2=="normal"] = "open"

    # this statement must set recipEdges before calling edgeRenderInfo
    # otherwise feedback loops are not shown properly.
    graphRenderInfo(g) <-  list(recipEdges=recipEdges)

    edgeRenderAttrs  <- setEdgeRenderInfo(edgeAttrs,
        list(arrowhead=arrowhead2, head=v2, tail=v1, 
        lwd=3, lty="solid"))
    edgeRenderInfo(g) <- edgeRenderAttrs

   # Set the edge Rendering in Rgraphviz
#   edgeRenderInfo(g) <- list(
#        col=edgeAttrs$color,
#        arrowhead=arrowhead2,
#        head=v2,
#        tail=v1,
#        label=edgeAttrs$label,
#        lwd=3,
#        lty="solid"    #this fails in some cases even with version >=1.33.1
#    )
    # hack for version of Rgraphviz 1.32. Set back the edgeAttr for the dot
    # output
    if (installed.packages()[,"Version"]["Rgraphviz"] <= "1.33.0"){
        edgeAttrs$color = savedEdgeAttrs
    }

    if (is.null(clusters)==TRUE){
        # finally, the layout for a R plot
        if (show==TRUE){
            x <- layoutGraph(g,layoutType="dot",recipEdges=recipEdges,attrs=attrs)
            renderGraph(x)
            nodeRenderInfo(x) <- nodeRenderAttrs
        }
        #edgeRenderInfo(x) <- edgeRenderAttrs
        edgeAttrs$lty=NULL    # why ?
        toDot(copyg, output_dot, nodeAttrs=nodeAttrs,edgeAttrs=edgeAttrs,attrs=attrs, recipEdges=recipEdges)

        # bug introduced in Rgraphviz 1.34 that set node attributes border.lwd
        # and border.color that are not understood by dot. Best solution is to
        # change Rgraphviz but large latency so we can change the written files
        # afterwards to change the dot file itself
        clean_dot(output_dot)
    }
    else{
        # finally, the layout for a R plot
        #attrs$node = nodeRenderInfo(g)
        #attrs$edge = edgeRenderInfo(g)
        if (show==TRUE){
            x <- layoutGraph(g,layoutType="dot",subGList=clusters,recipEdges=recipEdges,attrs=attrs)


            # known bug in Rgraphviz: renderInfo should be called again after the
            # layout otherwise some attributes are lost in the renderGraph (e.g.,
            # color)
            # note that rendering is now made on x variable (not g)
            nodeRenderInfo(x) <- nodeRenderAttrs
            edgeRenderInfo(x) <- edgeRenderAttrs
    
            renderGraph(x)
        }
        # and save into dot file.
        toDot(copyg, output_dot, nodeAttrs=nodeAttrs, subGList=clusters,
             attrs=attrs, recipEdges=recipEdges, edgeAttrs=edgeAttrs)
        # bug introduced in Rgraphviz 1.34 that set node attributes border.lwd
        # and border.color that are not understood by dot. Best solution is to
        # change Rgraphviz but large latency so we can change the written files
        # afterwards to change the dot file itself
        clean_dot(output_dot)


    }

    if (output != "STDOUT"){dev.off()}
    if (remove_dot==TRUE){file.remove(output_dot)}
    output = list(graph=g, attrs=attrs, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs,clusters=clusters, v1=v1, v2=v2, edges=edges)
}


# if a cnolist, or at least signals/stimuli, then we can create ranking for the
# layout
create_layout <- function(g, signals, stimuli){

    # this algo will tell us the distance between vertices
    # distMatrix columns contains the distance from a vertex to the others
    distMatrix <- floyd.warshall.all.pairs.sp(g)
    #distMatrix <- johnson.all.pairs.sp(g)
    distMatrix[is.infinite(distMatrix) == TRUE] <- -Inf # convert -Inf to be able to grep numbers

    # if signals provided by the user is empty, let us find ourself the nodes
    # with input degrees set to zero:
    if (is.null(stimuli)==TRUE){
        stimuli = nodes(g)[degree(g)$inDegree==0]
    }

    # we will need to know the sinks, which are defined by an outDegree equal to
    # zero. Yet, if signals are already provided, we want to restrict the sinks
    # to this subsets.
    if (is.null(signals)==TRUE){
        sinks  <- nodes(g)[degree(g)$outDegree == 0]
    }
    else{
        sinks  <- signals[degree(g, signals)$outDegree == 0]
    }

    # compute max rank for each column
    ranks <- apply(distMatrix, 2, max)
    mrank = max(ranks, na.rm=TRUE)-1  # -1 because we already know the sinks
    if (mrank == -Inf){
        print("Issue while computing max rank. Skipping the clustering step");
        return (NULL)
    }

    clusters <- list()
    if (mrank >= 1){ # otherwise, nothing to do.
        # for each different rank select the names of the column to create a
        # cluster
        for (rank in 1:mrank){  # starts at 1 although ranks starts at 0.
                                # However, 0 corresponds to the source that are
                                # known already
            # nodes found for a particular rank may contain a sink, that should
            # be removed
            nodes2keep = NULL
            nodes <- names(ranks[which(ranks==rank)])
            for (n in nodes){
                if (any(n==sinks) == FALSE){ nodes2keep <- c(nodes2keep, n)}
            }

            # may fail sometimes
            tryCatch({
                thisCluster <- subGraph(nodes2keep, g)
                thisGraph <-list(graph=thisCluster,cluster=FALSE,attrs=c(rank="same"))
                clusters[[length(clusters)+1]] =  thisGraph},
             error=function(e){"warning: clustering in the layout failed. carry on..."})

        }
    }
    # first the stimulis
    tryCatch(
        {
            tryCatch({
            clusterSource <- subGraph(stimuli, g);
            clusterSource<-list(graph=clusterSource,cluster=FALSE,attrs=c(rank="source"))},
                 error=function(e){print("error during clustering in
subGraph(stimuli, g)? Does the stimuli from your MIDAS are present in the model
?"); print(stimuli)})
            tryCatch(
                {clusters[[length(clusters)+1]] = clusterSource},
                 error=function(e){print("error in clusters2. should never be here")}
            )

        },
        error=function(e){print("warning: clustering the source/stimuli failed. carry on...")}
    )

    # then the signals keeping only those with outDegree==0
    tryCatch(
        {
            clusterSink <- subGraph(sinks, g)
            clusterSink <- list(graph=clusterSink, cluster=FALSE,  attrs=c(rank="sink"))
            tryCatch(
                {clusters[[length(clusters)+1]] = clusterSink},
                error=function(e){print("error in clusters3. should never be here")}
            )
        }, error=function(e){print("warning: clustering the sink failed. carry on...")}
    )


    return(clusters)
}


# Create the node attributes and save in a list to be used either by the
# plot function of the nodeRenderInfo function.
createNodeAttrs <- function(g, vertices, stimuli, signals, inhibitors, NCNO,
compressed, graphvizParams){


    nodeLabels = graphvizParams$nodeLabels
    nodeWidth = graphvizParams$nodeWidth
    nodeHeight = graphvizParams$nodeHeight
    andHeight = graphvizParams$andHeight
    andWidth = graphvizParams$andWidth
    mode = graphvizParams$mode

    # ----------------------------------------------- Build the node attributes list
    fillcolor <- list()
    color <- list()
    style <- list()  # the style of what is inside the node.
    lty <- list()    #style of the contour node
    height <-list()
    label <- list()
    width <- list()
    fixedsize <- list()
    shape <- list()

    # default. Must be filled in case no signals/stimulis/cnolist are provided
    # Default node attributes
    for (s in vertices){
        color[s] <- "black"         # color edge
        fillcolor[s] <- "white"     # fill color
        style[s] <- "filled,bold"
        lty[s] <- "solid"
        label[s] <- s
        height[s] <- nodeHeight
        width[s] <- nodeWidth
        fixedsize[s] <- FALSE
        shape[s] <- "rectangle"
    }

    # user can provide a list of labels for the nodes. Usuful to show attributes
    # of a node.
    if (is.null(nodeLabels)==FALSE){
        names(nodeLabels)<-vertices
        if (length(nodeLabels)!=length(vertices)){
            stop("if nodeLabels provided, it must be of same size as vertices")
        }
        else{
            for (s in vertices){
                label[s] = nodeLabels[s]
            }
        }
    }

    # The stimulis node
    for (s in stimuli){
        fillcolor[s] <- "olivedrab3";
        color[s] <- "olivedrab3";
        style[s] <- "filled"
        color[s] <- "black";
    }
    # The signal nodes
    for (s in signals){ #must be before the inhibitors to allow bicolors
        fillcolor[s] <- "lightblue";
        color[s] <-"lightblue";
        color[s] <- "black";
    }
    # The inhibitor node, that may also belong to the signal category.
    for (s in inhibitors){
        if (length(grep(s, signals))>=1){
            #fillcolor[s] <- "SkyBlue2"
            shape[s]="ellipse"
            style[s] <- "filled,bold"
            color[s] <-"orangered"
        } else {
            if (length(grep(s, stimuli))>=1){
                #fillcolor[s] <- "SkyBlue2"
                shape[s]="ellipse"
                style[s] <- "filled,bold"
                color[s] <-"orangered"
            }
            else{
                fillcolor[s] <- "orangered";
                color[s] <-"orangered";
                color[s] <- "black";
            }
        }
    }
    # The compressed nodes
    for (s in compressed){
        fillcolor[s] <- "white";
        color[s] <- "black";
        style[s] <- "dashed,bold";
        lty[s]="dashed";
    }
    # The NCNO node
    for (s in NCNO){
        fillcolor[s] <- "white";
        color[s] <- "grey";
        fillcolor[s] = "grey"
    }
    # the and gate nodes
    for (s in vertices){
        if (length(grep("and", s))>=1){
            color[s] = "black"
            if (mode=="sbgn"){
                fillcolor[s] = "white"
                label[s] = "and"
            } else{
                fillcolor[s] = "black"
                label[s] = ""
            }
            width[s]=andWidth
            height[s]=andHeight
            fixedsize[s]=FALSE
            shape[s]="circle"

        if (degree(g)$inDegree[s]==3){
            #fillcolor[s] = "blue"
            if (mode=="classic"){
                shape[s]="triangle"
            }
            width[s]=andWidth
            height[s]=andHeight
        }
        if (degree(g)$inDegree[s]==4){
            #fillcolor[s] = "blue"
            if (mode=="classic"){
                shape[s]="rectangle"
            }
            width[s]=andWidth
            height[s]=andHeight
        }

        }
    }
    nodeAttrs <- list(fillcolor=fillcolor, color=color, label=label, width=width, height=height,
        style=style, lty=lty, fixedsize=fixedsize, shape=shape)
    return(nodeAttrs)
}


# Create the node attributes and save in a list to be used either by the
# plot function of the edgeRenderInfo function.
createEdgeAttrs <- function(v1, v2, edges, BStimes ,Integr, user_edgecolor,
    view_empty_edge=TRUE){

    edgewidth_c = 3 # default edge width

    # The edge attributes
    arrowhead <- list()
    edgecolor <- list()
    edgewidth <- list()
    label <- list()
    toremove <- list()
    lty <- list() # not used yet.


    for (i in 1:length(edges)){
        edgename = paste(v1[i], "~", v2[i], sep="")
        edgewidth[edgename] = edgewidth_c    # default edgewidth
        label[edgename] = ""                 # no label on the edge by default
        lty[edgename] = "solid"
        edgecolor[edgename] = "black"        # set a default (useless but safe)

        if (edges[i] == 1){
           arrowhead[edgename] <- "normal"
           edgecolor[edgename] <- user_edgecolor
           if (Integr[i]==1){edgecolor[edgename] <- "purple"}
        }
        else if (edges[i] == -1){
           arrowhead[edgename] <- "tee"
           edgecolor[edgename] <- "red"
           if (Integr[i]==1){edgecolor[edgename] <- "purple"}
        }
        else{
           arrowhead[edgename] <- "normal"
           edgecolor[edgename] <- "blue"
        }

        # BStimes contains the bitstring. color the edges according to its value
        v = (BStimes[i]*100)%/%1

        # width of the edges
        if (v != 100){
            # first, let us build the color
            if (edgecolor[edgename] == 'red'){
                # if red, go from red to light pink color according to v value
                if (view_empty_edge==TRUE){
                    print("showing red link even if empty")
                    color = rgb(1,1-(max(20,v)/100),1-(max(20,v)/100))
                }
                else{
                    color = rgb(1,1-v/100,1-v/100)
                }
            }
            else if (edgecolor[edgename] == 'black'){
                # if black, go from grey dark to grey light color according to v value
                if (view_empty_edge==TRUE){
                    color = paste("grey", as.character(100.-max(20,v)), sep="")
                }
                else{
                    color = paste("grey", as.character(100.-v), sep="")
                }
            }
            else{
                # otherwise, just keep the color identical and only add label
                color = edgecolor[edgename]
            }
            # now, we fill the fields with color and width
            if (v == 0){
                edgewidth[edgename]  = edgewidth_c*(10./100)
                edgecolor[edgename] = color
                # and gates that have no link are removed.
                if (length(grep("and", edgename))>=1){
                    toremove <- append(toremove, edgename)
                }
            }
            else{
              edgecolor[edgename] = color
              edgewidth[edgename] = edgewidth_c*(v/100)
              if (v!=0){
                label[edgename] = as.character(v)
                }else{label[edgename]="0"}
            }
        }
        else {
              #already set at the beginning of the loop: width = edgewidth_c
        }
    }

    edgeAttrs <- list(color=edgecolor,arrowhead=arrowhead,
        penwidth=edgewidth,label=label, lty=lty)

    return(list(toremove=toremove, edgeAttrs=edgeAttrs))
}



clean_dot <- function(filename)
{
    savedata = readLines(filename)
    data = savedata
    for (i in 1:length(data)){
        data[i] = sub("border.lwd=1,", "", data[i])
        data[i] = sub("border.color=black", "", data[i])
    }

    # if died, need to save back the "savedata"
    # otherwise, we can overwrite model.dot
    tryCatch(write(data, file=filename), error=function(e){print("Could not scan the dot file for cleaning (border.lwd and border.color). "); })
}

# Create the list of attributes for the nodes
setNodeRenderInfo <- function(nodeAttrs, extraAttrs)
{
    attrs  <- list(
        fill=nodeAttrs$fillcolor,
        col=nodeAttrs$color,
        style=nodeAttrs$style,
        lty=extraAttrs$lty,
        lwd=extraAttrs$lwd,   
        label=nodeAttrs$label,
        shape=nodeAttrs$shape,
        cex=extraAttrs$cex,
        fontsize=extraAttrs$fontsize,
        iwidth=nodeAttrs$width,
        iheight=nodeAttrs$height,
        fixedsize=extraAttrs$fixedsize)
    return(attrs)
}

# Create the list of attributes for the edges
setEdgeRenderInfo <- function(edgeAttrs, extraAttrs)
{
   attrs <- list(
        col=edgeAttrs$color,
        arrowhead=extraAttrs$arrowhead,
        head=extraAttrs$v2,
        tail=extraAttrs$v1,
        label=edgeAttrs$label,
        lwd=extraAttrs$lwd,
        lty=extraAttrs$lty    #this fails in some cases even with version >=1.33.1
    )
    return(attrs)
}
