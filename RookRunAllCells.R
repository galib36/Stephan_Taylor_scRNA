#!/usr/bin/env Rscript

library(Rook)
library("dplyr")
library("edgeR")
library("ggplot2")
library("cowplot")
library("BASiCS")
library("scde")
library(repr)
library(biomaRt)
library(GO.db)
library(extRemes)
library(Lmoments)
library(distillery)
library(car)



load("/home/baker/Rna-seq_Data-Analysis/Louisa_Nelson_Single_Cell_Analysis/Heterogeneity_AnalysisReAssigned.RData")

myPort <- 1491
myInterface <- "0.0.0.0"
status <- -1

# R 2.15.1 uses .Internal, but the next release of R will use a .Call.
# Either way it starts the web server.
if (as.integer(R.version[["svn rev"]]) > 59600) {
  status <- .Call(tools:::startHTTPD, myInterface, myPort)
} else {
  status <- .Internal(startHTTPD(myInterface, myPort))
}

if (status == 0) {
  unlockBinding("httpdPort", environment(tools:::startDynamicHelp))
  assign("httpdPort", myPort, environment(tools:::startDynamicHelp))
  
  s <- Rhttpd$new()
  s$listenAddr <- myInterface
  s$listenPort <- myPort
  
  # Change this line to your own application. You can add more than one
  # application if you like
 # s$add(name = "Louisa_AllCell_Cluster", app = make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols = rbind(groups = cutree(hc, 2)), cell.clustering = hc, title = "AllCell cluster", embedding = tSNE.pagoda$Y))
 s$add(name = "Louisa_AllCell_Reassigned_Cluster", app = make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols = rbind(groups = cutree(hc, 4)), cell.clustering = hc, title = "AllCell Reassigned cluster"))
  
  # Now make the console go to sleep. Of course the web server will still be
  # running.
  while (TRUE) Sys.sleep(24 * 60 * 60)
}

# If we get here then the web server didn't start up properly
warning("Oops! Couldn't start Rook app")


