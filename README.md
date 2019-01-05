# "Textual Spanning" 

## Reproduction Guide

Dustin S. Stoltz and Marshall A. Taylor

This is the code and data necessary to reproduce the graphs and plots for Stoltz and Taylor (2019) "Textual Spanning: Finding Discursive Holes in Text Networks" in _Socius_.

The measure of textual spanning we propose works on a document by document similarity matrix. The basic data structure in text analysis is an MxN matrices of documents by terms, n-grams, parts of speech, topics and so on. The rows, therefore, are vector representations of the text and can be easily compared for similarity, usually with cosine similarity. In our empirical demonstration, we use a topic model solution to generate a similarity matrix. The result is a one-mode document by document matrix, which can easily be interpreted as a weighted adjacency matrix amendable to network metrics. 

The `textSpan` function takes this document by document similarity matrix and outputs a document specific measure which increases when a document is similar to documents which are not also similar to each other. This is defined by the following equations:

<img src="https://latex.codecogs.com/gif.latex?S_i%20%3D%20%5Csum_j%20%5Cleft%20%28%20p_%7Bij%7D%20&plus;%20%5Csum_q%20%5Cfrac%7Bp_%7Bqj%7D%7D%7Bp_%7Biq%7D%7D%20%5Cright%20%29%5E2"/>

## textSpan Function
----

``` r
    textSpan <- function(A, alpha=1){
        # zero the diagonal of the similarity matrix
        diag(A) <- 0 
        # get denominator (i.e. weighted degree), adjustable by alpha
        den <- (rowSums(A != 0)) * ((rowSums(A)/
                (rowSums(A != 0)))^alpha)
        # divide A by den to get proportional similarities,
        # equation (2) in the paper
        PS <- A/den
        # sum paths of length two using dot product
        PS2 <- PS%*%PS
        # cannot divide matrices, so find 
        # inverse of PS to multiple in next step
        iPS <- solve(PS)
        # remove zero edges and calculate the dyadic spanning 
        # score, which is the term in the parentheses in equation (1)
        SP <- (PS + (PS2*(as.numeric(iPS>0))))^2  
        # calculate cumulative spanning score for each vertex
        cSP <- rowSums(SP) 
        # standardize and invert the scores, equation (3) in the paper
        cSP <- ((cSP-mean(cSP))/sd(cSP))*-1
      }
 ```
 
## Simulated and Empirical Examples
-----
### Load and Prepare Data
Download CSVs, set your working directory, and load these additional packages used in the analysis and visualization:
``` r
        #install.packages("pacman")
        library(pacman)
        pacman::p_load(ggnetwork, ggplot2, ggpubr, igraph,
                        corrplot, tnet, PerformanceAnalytics,
                        stm, reshape2, qgraph, intergraph, 
                        install = T)

```
Load simulated datasets:
``` r 
      sim0 <- read.csv("2_sim_ring_0.csv", stringsAsFactors=FALSE, row.names=1)
      sim1 <- read.csv("2_sim_ring_b.csv", stringsAsFactors=FALSE, row.names=1)
      sim2 <- read.csv("2_sim_clique.csv", stringsAsFactors=FALSE, row.names=1)
      sim3 <- read.csv("2_sim_core_periph.csv", stringsAsFactors=FALSE, row.names=1
```
We use the CMU political blogs dataset which are used with the `stm` package. Because n=13,000 in this dataset, for simplicity we randomly sample 100 blog posts which we use in the following. We offer a subset of the document by term matrix and the document by topic probability matrix, the latter of which is based on the pre-processed topic model solution provided by CMU team (the RData file can be downloaded here: http://goo.gl/VPdxlS).

Load pre-fitted topic model solution on subset of 100 randomly sampled blogposts:
``` r
    tms  <- as.matrix(read.csv("2_subset_topic_solution.csv", 
                                    stringsAsFactors=FALSE, row.names=1))
```
Calculate cosine similarities between documents based on topic model solution
``` r
    cos.tms <- tcrossprod(tms / sqrt(rowSums(tms * tms)))
```

Create tnet object to calculate Opsahl et al's weighted measures
``` r
    # simulated ring graphs
    tn.sim0 <- as.tnet(sim0, type="weighted one-mode tnet")
    tn.sim1 <- as.tnet(sim1, type="weighted one-mode tnet")
    # topic model graph
    tn.tms <- as.tnet(cos.tms, type="weighted one-mode tnet")
```
Create iGraph object for visualization
``` r
    # simulated graphs
    sim.net0 <- graph.adjacency(as.matrix(sim0), diag=F, mode="lower", weighted=T)
    sim.net1 <- graph.adjacency(as.matrix(sim1), diag=F, mode="lower", weighted=T)
    sim.net2 <- graph.adjacency(as.matrix(sim2), diag=F, mode="lower", weighted=T)
    sim.net3 <- graph.adjacency(as.matrix(sim3), diag=F, mode="lower", weighted=T)
    # topic model graph
    g.tms <- graph_from_adjacency_matrix(cos.tms, 
                                        weighted = TRUE, 
                                        diag = FALSE,
                                        add.colnames = NULL)
```

### Calculating Measures
``` r
    # Textual SPANNING --------------------------------------------------------
        # simulated graphs
        # alpha set to 1.0
        V(sim.net0)$span.1.0 <- textSpan(as.matrix(sim0), alpha = 1.0)
        V(sim.net1)$span.1.0 <- textSpan(as.matrix(sim1), alpha = 1.0)
        V(sim.net2)$span.1.0 <- textSpan(as.matrix(sim2), alpha = 1.0)
        V(sim.net3)$span.1.0 <- textSpan(as.matrix(sim3), alpha = 1.0)

        # topic model graph
        # alpha = 1.0
        V(g.tms)$span.1.0 <- textSpan(cos.tms, alpha = 1.0)
    
    # Opsahl et al's WEIGHTED DEGREE
        # simulated ring graphs
        # alpha = 1.0
        V(sim.net0)$deg.1.0 <- as.data.frame(degree_w(tn.sim0, alpha=1.0))$output
        V(sim.net1)$deg.1.0 <- as.data.frame(degree_w(tn.sim1, alpha=1.0))$output
        # topic model graph
        # alpha = 1.0
        V(g.tms)$deg.1.0 <- as.data.frame(degree_w(tn.tms, alpha=1.0))$output
    
    # Opsahl et al's WEIGHTED BETWEENNESS 
        # alpha = 1.0
        V(g.tms)$bet.1.0 <- as.data.frame(betweenness_w(tn.tms, alpha=1.0))$betweenness

    # Opsahl et al's WEIGHTED CLOSENESS
        # alpha = 1.0
        V(g.tms)$clo.1.0 <- as.data.frame(closeness_w(tn.tms, alpha=1.0))$n.closeness

    # Add blog post lengths as vertex attribute
        length <- read.csv("2_doc_lengths.csv")
        V(g.tms)$length <- length$length
```
### Generating Graphs and Plots

#### Simulated Graphs
``` r
    # Prepare simulated graphs for ggplot:
        ## disconnected ring
        l0 <- layout_in_circle(sim.net0)
        sim.net0 <- ggnetwork(sim.net0, layout=l0, weight="weight")
        ## connected ring
        l1 <- layout_in_circle(sim.net1)
        sim.net1 <- ggnetwork(sim.net1, layout=l1, weight="weight")
        ## bifurcated clique
        l2 <- layout_with_fr(sim.net2)
        sim.net2 <- ggnetwork(sim.net2, layout=l2, weight="weight")
        ## core-periphery 
        l3 <- layout_with_fr(sim.net3)
        sim.net3 <- ggnetwork(sim.net3, layout=l3, weight="weight")
        
    # GRAPHS --------------------------------------------------------------
    # Disconnected Ring Graph
    # text SPANNING with alpha = 1.0 
    set.seed(123)
    net0.1.0 <- ggplot(sim.net0, aes(x = x, y = y, xend = xend, yend = yend)) +
                geom_edges(data=sim.net0[which(sim.net0$weight>=.8),], alpha=.5, color="black") +
                geom_edges(data=sim.net0[which(sim.net0$weight<=.5),], alpha=.2, color="gray50") +
                geom_nodes(data=sim.net0, aes(x, y, color=span.1.0), size=3) +
                theme_blank() +
                theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                        legend.title=element_text(face="bold"),
                        axis.text = element_blank(),
                        axis.title = element_blank()) +
                labs(title=expression(bold(paste("Disconnected, Spanning ", alpha, " = 1.0")))) +
                scale_color_gradient(name="", 
                                    limits=c(min(sim.net0$span.1.0), max(sim.net0$span.1.0)),
                                    breaks=c(min(sim.net0$span.1.0), max(sim.net0$span.1.0)), 
                                    labels=c("Low","High"),
                                    low = "#fdbf11", high = "#1696d2") +
                geom_nodelabel_repel(aes(label = ifelse(sim.net0$span.1.0>=quantile(span.1.0, c(.0)),
                                        as.character(round(span.1.0, digits=3)),'')),
                                        box.padding = unit(1, "lines")) 
    # Connected Ring Graph
    # text SPANNING with alpha = 1.0 
    set.seed(123)
    net1.1.0 <- ggplot(sim.net1, aes(x = x, y = y, xend = xend, yend = yend)) +
                geom_edges(data=sim.net1[which(sim.net1$weight>=.8),], alpha=.5, color="black") +
                geom_edges(data=sim.net1[which(sim.net1$weight<=.5),], alpha=.2, color="gray50") +
                geom_nodes(data=sim.net1, aes(x, y, color=span.1.0), size=3) +
                theme_blank() +
                theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                        legend.title=element_text(face="bold"),
                        axis.text = element_blank(),
                        axis.title = element_blank()) +
                labs(title=expression(bold(paste("Connected, Spanning ", alpha, " = 1.0")))) +
                scale_color_gradient(name="", 
                                    limits=c(min(sim.net1$span.1.0), max(sim.net1$span.1.0)),
                                    breaks=c(min(sim.net1$span.1.0), max(sim.net1$span.1.0)), 
                                    labels=c("Low","High"),
                                    low = "#fdbf11", high = "#1696d2") +
                geom_nodelabel_repel(aes(label = ifelse(sim.net1$span.1.0>=quantile(span.1.0, c(.0)),
                                        as.character(round(span.1.0, digits=3)),'')),
                                        box.padding = unit(1, "lines")) 
    # -----------------------
    # Disconnected Ring Graph
    # weighted DEGREE with alpha = 1.0 
    set.seed(123)
    net0.1.0d <- ggplot(sim.net0, aes(x = x, y = y, xend = xend, yend = yend)) +
                geom_edges(data=sim.net0[which(sim.net0$weight>=.8),], alpha=.5, color="black") +
                geom_edges(data=sim.net0[which(sim.net0$weight<=.5),], alpha=.2, color="gray50") +
                geom_nodes(data=sim.net0, aes(x, y, color=deg.1.0), size=3) +
                theme_blank() +
                theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                        legend.title=element_text(face="bold"),
                        axis.text = element_blank(),
                        axis.title = element_blank()) +
                labs(title=expression(bold(paste("Disconnected, Degree ", alpha, " = 1.0")))) +
                scale_color_gradient(name="", 
                                    limits=c(min(sim.net0$deg.1.0), max(sim.net0$deg.1.0)),
                                    breaks=c(min(sim.net0$deg.1.0), max(sim.net0$deg.1.0)), 
                                    labels=c("Low","High"),
                                    low = "#fdbf11", high = "#1696d2") +
                geom_nodelabel_repel(aes(label = ifelse(sim.net0$deg.1.0>=quantile(deg.1.0, c(.01)),
                                        as.character(round(deg.1.0, digits=3)),'')),
                                        box.padding = unit(1, "lines")) 
    # Connected Ring Graph
    # weighted DEGREE with alpha = 1.0 
    set.seed(123)
    net1.1.0d <- ggplot(sim.net1, aes(x = x, y = y, xend = xend, yend = yend)) +
                geom_edges(data=sim.net1[which(sim.net1$weight>=.8),], alpha=.5, color="black") +
                geom_edges(data=sim.net1[which(sim.net1$weight<=.5),], alpha=.2, color="gray50") +
                geom_nodes(data=sim.net1, aes(x, y, color=deg.1.0), size=3) +
                theme_blank() +
                theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                        legend.title=element_text(face="bold"),
                        axis.text = element_blank(),
                        axis.title = element_blank()) +
                labs(title=expression(bold(paste("Connected, Degree ", alpha, " = 1.0")))) +
                scale_color_gradient(name="", 
                                    limits=c(min(sim.net1$deg.1.0), max(sim.net1$deg.1.0)),
                                    breaks=c(min(sim.net1$deg.1.0), max(sim.net1$deg.1.0)), 
                                    labels=c("Low","High"),
                                    low = "#fdbf11", high = "#1696d2") +
                geom_nodelabel_repel(aes(label = ifelse(sim.net1$deg.1.0>=quantile(deg.1.0, c(.01)),
                                        as.character(round(deg.1.0, digits=3)),'')),
                                        box.padding = unit(1, "lines")) 

    ## Arrange Ring Graphs for comparison
        pdf("Figure_2_Ring_Graphs.pdf")
        ggarrange(net0.1.0, net1.1.0, net0.1.0d, net1.1.0d,
            ncol=2, nrow=2, common.legend=T, legend="left")
        dev.off()

    ## --------------------------------------------------------------
    # Bifurcated Clique Graph (Spanning)
    # alpha = 1.0
    set.seed(123)
    net2.1.0 <- ggplot(sim.net2, aes(x = x, y = y, xend = xend, yend = yend)) +
                    geom_edges(data=sim.net2[which(sim.net2$weight==.8),], alpha=.5, color="black") +
                    geom_edges(data=sim.net2[which(sim.net2$weight==.1),], alpha=.1, color="gray50") +
                    geom_nodes(data=sim.net2, aes(x, y, color=span.1.0), size=3) +
                    theme_blank() +
                    theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                            legend.title=element_text(face="bold"),
                            axis.text = element_blank(),
                            axis.title = element_blank()) +
                    labs(title=expression(bold(paste("Bifurcated Clique, Spanning ",alpha, " = 1.0")))) +
                    scale_color_gradient(name="", 
                                        limits=c(min(sim.net2$span.1.0), max(sim.net2$span.1.0)),
                                        breaks=c(min(sim.net2$span.1.0), max(sim.net2$span.1.0)), 
                                        labels=c("Low","High"),
                                        low = "#fdbf11", high = "#1696d2") +
                    geom_nodelabel_repel(aes(label = ifelse(sim.net2$span.1.0>=quantile(span.1.0, c(.90)),
                                            as.character(round(span.1.0, digits=3)),'')),
                                            box.padding = unit(1, "lines"))
    
    pdf("Figure_3_Clique_Graph.pdf", width=8, height=6)
    net2.1.0
    dev.off()
 

    ## --------------------------------------------------------------
    # Core-Periphery Graph (Spanning)
    # alpha = 1.0
    set.seed(123)
    net3.1.0 <- ggplot(sim.net3, aes(x = x, y = y, xend = xend, yend = yend)) +
                    geom_edges(data=sim.net3[which(sim.net3$weight==.8),], alpha=.5, color="black") +
                    geom_edges(data=sim.net3[which(sim.net3$weight==.1),], alpha=.1, color="gray50") +
                    geom_nodes(data=sim.net3, aes(x, y, color=span.1.0), size=3) +
                    theme_blank() +
                    theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                            legend.title=element_text(face="bold"),
                            axis.text = element_blank(),
                            axis.title = element_blank()) +
                    labs(title=expression(bold(paste("Core-Periphery, ", alpha, " = 1.0")))) +
                    scale_color_gradient(name="", 
                                        limits=c(min(sim.net3$span.1.0), max(sim.net3$span.1.0)),
                                        breaks=c(min(sim.net3$span.1.0), max(sim.net3$span.1.0)), 
                                        labels=c("Low","High"),
                                        low = "#fdbf11", high = "#1696d2") +
                    geom_nodelabel_repel(aes(label = ifelse(sim.net3$span.1.0>=quantile(span.1.0, c(.80)),
                                            as.character(round(span.1.0, digits=3)),'')),
                                            box.padding = unit(1, "lines"))
    
    pdf("Figure_4_Core_Graph.pdf", width=8, height=6)
    net3.1.0
    dev.off()
```
#### Topic Model Solution Graphs
``` r
    ## Prepare Topic Model Graph Layout
        # Removing edges for visualization
        p.tms <- igraph::delete.edges(g.tms, E(g.tms)[weight<.6])
        # 
        l5 <- layout_with_fr(p.tms, niter=2000)
        p.tms <- ggnetwork(p.tms, layout=l5, weight="weight")

    # Graph colored by TEXTUAL SPANNING Score, alpha = 1.0
    set.seed(786)
    g.span.1.0 <- ggplot(p.tms, aes(x = x, y = y, xend = xend, yend = yend)) +
                    geom_edges(alpha=.4, color="gray75") +
                    geom_nodes(aes(x, y, size=length, color=span.1.0)) +
                    theme_blank() +
                    theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                            legend.title=element_text(face="bold"),
                            axis.text = element_blank(),
                            axis.title = element_blank()) +
                    labs(title=expression(bold(paste("Textual Spanning, ", alpha, " = 1.0")))) +
                    scale_color_gradient(name="", limits=c(min(p.tms$span.1.0), 
                                                                            max(p.tms$span.1.0)),
                                        breaks=c(min(p.tms$span.1.0), max(p.tms$span.1.0)), 
                                                    labels=c("Low","High"),
                                        low = "#fdbf11", high = "#1696d2") +
                    guides(size = guide_legend(reverse=T, title="Length")) +
                    geom_nodelabel_repel(aes(label = ifelse(p.tms$span.1.0>=quantile(span.1.0, c(.95)),
                                                            as.character(round(span.1.0, digits=3)),'')),
                                        box.padding = unit(1, "lines")) +
                    geom_nodelabel_repel(aes(label = ifelse(p.tms$span.1.0<=quantile(span.1.0, c(.01)),
                                                            as.character(round(span.1.0, digits=3)),'')),
                                        box.padding = unit(1, "lines")) 
    
        pdf("Figure_6_TMS_Span.pdf", width=11, height=8.5)
        g.span.1.0
        dev.off()
        
    # Same as above, but without labels
    set.seed(786)
    g.span.1.0b <- ggplot(p.tms, aes(x = x, y = y, xend = xend, yend = yend)) +
                    geom_edges(alpha=.4, color="gray75") +
                    geom_nodes(aes(x, y, size=length, color=span.1.0)) +
                    theme_blank() +
                    theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                            legend.title=element_text(face="bold"),
                            axis.text = element_blank(),
                            axis.title = element_blank()) +
                    labs(title=expression(bold(paste("Textual Spanning, ", alpha, " = 1.0")))) +
                    scale_color_gradient(name="", limits=c(min(p.tms$span.1.0), 
                                                                            max(p.tms$span.1.0)),
                                        breaks=c(min(p.tms$span.1.0), max(p.tms$span.1.0)), 
                                        labels=c("Low","High"),
                                        low = "#fdbf11", high = "#1696d2") +
                    guides(size = guide_legend(reverse=T, title="Length"))

    ## --------------------------------------------------------------
    # Opsahl et al's WEIGHTED DEGREE, alpha = 1.0
    set.seed(786)
    g.deg.1.0 <- ggplot(p.tms, aes(x = x, y = y, xend = xend, yend = yend)) +
                    geom_edges(alpha=.4, color="gray75") +
                    geom_nodes(aes(x, y, size=length, color=deg.1.0)) +
                    theme_blank() +
                    theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                            legend.title=element_text(face="bold"),
                            axis.text = element_blank(),
                            axis.title = element_blank()) +
                    labs(title=expression(bold(paste("Weighted Degree, ", alpha, " = 1.0")))) +
                    scale_color_gradient(name=" ", limits=c(min(p.tms$deg.1.0), 
                                                                            max(p.tms$deg.1.0)),
                                        breaks=c(min(p.tms$deg.1.0), max(p.tms$deg.1.0)), 
                                                    labels=c("Low","High"),
                                        low = "#fdbf11", high = "#1696d2") +
                    guides(size = guide_legend(reverse=T, title="Length"))

    ## --------------------------------------------------------------
    # Opsahl et al's WEIGHTED BETWEENNESS 
    set.seed(786)
    g.bet.1.0 <- ggplot(p.tms, aes(x = x, y = y, xend = xend, yend = yend)) +
                    geom_edges(alpha=.4, color="gray75") +
                    geom_nodes(aes(x, y, size=length, color=bet.1.0)) +
                    theme_blank() +
                    theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                            legend.title=element_text(face="bold"),
                            axis.text = element_blank(),
                            axis.title = element_blank()) +
                    labs(title=expression(bold(paste("Weighted Betweeness, ", alpha, " = 1.0")))) +
                    scale_color_gradient(name=" ", limits=c(min(p.tms$bet.1.0), 
                                                                            max(p.tms$bet.1.0)),
                                        breaks=c(min(p.tms$bet.1.0), max(p.tms$bet.1.0)), 
                                                    labels=c("Low","High"),
                                        low = "#fdbf11", high = "#1696d2") +
                    guides(size = guide_legend(reverse=T, title="Length"))

    ## --------------------------------------------------------------
    # Opsahl et al's WEIGHTED CLOSENESS
    set.seed(786)
    g.clo.1.0 <- ggplot(p.tms, aes(x = x, y = y, xend = xend, yend = yend)) +
                    geom_edges(alpha=.4, color="gray75") +
                    geom_nodes(aes(x, y, size=length, color=clo.1.0)) +
                    theme_blank() +
                    theme(legend.position="left", plot.title=element_text(face="bold", hjust=.5),
                            legend.title=element_text(face="bold"),
                            axis.text = element_blank(),
                            axis.title = element_blank()) +
                    labs(title=expression(bold(paste("Weighted Closeness, ", alpha, " = 1.0")))) +
                    scale_color_gradient(name="  ", limits=c(min(p.tms$clo.1.0), 
                                                                            max(p.tms$clo.1.0)),
                                        breaks=c(min(p.tms$clo.1.0), max(p.tms$clo.1.0)), 
                                                    labels=c("Low","High"),
                                        low = "#fdbf11", high = "#1696d2") +
                    guides(size = guide_legend(reverse=T, title="Length"))
    
    ## Arrange Graphs w/ Opsahl et al's measures side by side
        pdf("Figure_7_TMS_Opsahl.pdf", width=11, height=10)
        ggarrange(g.span.1.0b, g.deg.1.0, g.bet.1.0, g.clo.1.0,
                    ncol=2, nrow=2, common.legend=T, legend="bottom")
        dev.off()
```

#### Correlation Table Comparing Spanning Scores with Centrality Measures
``` r
        cor <- data.frame()[1:100, ]
        head(cor)
        cor$tms.span.1.0 <- V(g.tms)$span.1.0
        cor$tms.deg.1.0 <- V(g.tms)$deg.1.0
        cor$tms.bet.1.0 <- V(g.tms)$bet.1.0
        cor$tms.clo.1.0 <- V(g.tms)$clo.1.0 
        cor$length <- length$length
        colnames(cor) <- c("Spanning","Degree", "Betweenness","Closeness", "Length")

    ## Calculate correlation using Performance Analytics package
        pdf("Figure_8_Correlations.pdf", width=7, height=5)
        chart.Correlation(cor, histogram=T, pch=19, method="pearson")
        dev.off()
```
#### Generate Box Plot Comparing Top 5 Spanning Posts across Centrality Measures
``` r
    ## Create data frame
        names <- V(g.tms)$name
        spanning <- V(g.tms)$span.1.0
        degree <- V(g.tms)$deg.1.0
        between <- V(g.tms)$bet.1.0
        close <- V(g.tms)$clo.1.0
        span.mat <- cbind(names, spanning, degree, between, close)
        span.mat <- as.data.frame(span.mat)
        span.mat <- span.mat[order(-spanning),]
        # head(span.mat, n=5)

    ## Box plots
        span.mat$degree <- as.numeric(as.character(span.mat$degree))
        span.mat$between <- as.numeric(as.character(span.mat$between))
        span.mat$close <- as.numeric(as.character(span.mat$close))

        span.mat$sddegree <- ((span.mat$degree - mean(span.mat$degree))/sd(span.mat$degree))
        span.mat$sdbetween <- ((span.mat$between - mean(span.mat$between))/sd(span.mat$between))
        span.mat$sdclose <- ((span.mat$close - mean(span.mat$close))/sd(span.mat$close))

        span.mat2 <- melt(span.mat, id.vars="names")
        span.mat2$color[which(span.mat2$variable=="spanning")] <- "#1696d2" 
        span.mat2$color[which(span.mat2$variable=="sddegree")] <- "#ec008b" 
        span.mat2$color[which(span.mat2$variable=="sdbetween")] <- "#fdbf11"
        span.mat2$color[which(span.mat2$variable=="sdclose")] <- "#5c5859"
        span.mat2$color[which(span.mat2$names=="5781" | span.mat2$names=="1926" |
                                span.mat2$names=="9850" | span.mat2$names=="13017" |
                                span.mat2$names=="12698")] <- "black"
        span.mat2$shape <- NA
        span.mat2$shape[which(span.mat2$color=="black")] <- 1
        span.mat2$shape[which(span.mat2$color!="black")] <- 0
    ##

        box <- ggplot(span.mat2[which(span.mat2$variable=="spanning" | span.mat2$variable=="sddegree" |
                                span.mat2$variable=="sdbetween" | span.mat2$variable=="sdclose"),], 
            aes(x=as.factor(variable), y=as.numeric(value), fill=as.factor(variable))) +
        geom_boxplot(aes(fill=variable), alpha=.5, outlier.shape=NA) +
        geom_jitter(aes(color=as.factor(color), shape=as.factor(shape),
                        size=shape)) +
        ylab("Centrality (Standardized)") + xlab("") +
        theme_bw() +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                legend.title=element_text(face="bold")) +
        scale_fill_manual(name="Centrality\nMeasures",
                            breaks=c("spanning","sddegree","sdbetween","sdclose"),
                            labels=c("Spanning","W Degree","W Between","W Closeness"),
                            values=c("#1696d2","#ec008b","#fdbf11","#55b748")) +
        scale_color_manual(values=c("#1696d2","#55b748","#ec008b","#fdbf11","black")) +
        scale_size_continuous(range = c(1,2)) +
        guides(color=F) +
        guides(shape=F) +
        guides(size=F)
    ##
        pdf("Figure_9_box_plot.pdf", width=8, height=6)
        box
        dev.off()
 ```
 #### Generate Topic Correlation Network 
 ``` r
        cor.tms <- cor(tms)
        cor.tms <- ifelse(cor.tms<0, 0, cor.tms)
        new.names <- paste("T", seq(1:20), sep="")
        ##
        qgraph(cor.tms, graph="cor", layout="spring", posCol="black",
                        threshold=.1, labels=new.names, 
                        color=ifelse(rownames(cor.tms)=="V18" | 
                                    rownames(cor.tms)=="V14", "#fdbf11","#d2d2d2"), 
                        vsize=5, 
                        filename="Figure_10_Correlation_Network",
                        filetype="pdf")

```
# END # -------------------------------------------------------------------


