taxa.sig.dend.surv <- function(out, tax.tab, layout.type = "twopi", species.include = "FALSE") {
  
  print("come in success_dend")
  
  names(out) <-c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value", "Q-value"), 1, 5)
  ci.tab.all <- matrix(NA, 1, 1)
  for (i in 1:6) {
      ind.sig <- which(out[[i]]$Q.value < 0.05)
      if (length(ind.sig) >= 1) {
        sig.out <- out[[i]][ind.sig,]
        taxa <- rownames(sig.out)
        text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,1], 3), nsmall = 3), 
                          format(round(sig.out[,"P.value"], 3), nsmall = 3), 
                          format(round(sig.out[,"Q.value"], 3), nsmall = 3))
        ci.tab <- cbind(sig.out[,1])
        text.tab.all <- rbind(text.tab.all, text.tab)
        ci.tab.all <- rbind(ci.tab.all, ci.tab)
      }
  }
  
  text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
  tax.tab.num <- as.matrix(as.data.frame(tax.tab))
  tax.tab.numNA <- as.matrix(as.data.frame(tax.tab))
  
  for (i in 1:6) {
    ind <- which(tax.tab.num[,i+1] == "NANANA")
    if (length(ind) >= 1) {
      tax.tab.num[ind,i+1] <- NA
    }
  }
  rank <- 1
  for (tax.rank in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    ind.sig <- which(text.tab.all[,2] == tax.rank)
    rank.list <- c("p_", "c_", "o_", "f_", "g_", "s_")
    if (length(ind.sig) >= 1) {
      sig.taxa <- text.tab.all[ind.sig,3]
      ind.non.sig <- which(tax.tab.num[,tax.rank] %in% sig.taxa)
      
      if(length(ind.non.sig) >0) {
        non.sig.taxa <- names(table(tax.tab.num[-ind.non.sig,tax.rank]))
        for (i in 1:length(non.sig.taxa)) {
          ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
          tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
        }
      }
      
      for (i in 1:length(sig.taxa)) {
        ind <- which(text.tab.all[,3] == sig.taxa[i])
        ind.rank <- which(text.tab.all[,2] == tax.rank)
        ind <- ind[ind %in% ind.rank]
        sig.taxa.num <- text.tab.all[ind ,1]
        ind.num <- which(tax.tab.num[,tax.rank] == sig.taxa[i])
        tax.tab.num[ind.num,tax.rank] <- sig.taxa.num
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        
      }
      rank <- rank+1
      
    } else {
      non.sig.taxa <- names(table(tax.tab.num[,tax.rank]))
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
      }
      rank <- rank+1
    }
  }
  all.phylum <- tax.tab.num[,"Phylum"]
  all.class <- tax.tab.num[,"Class"]
  all.order <- tax.tab.num[,"Order"]
  all.family <- tax.tab.num[,"Family"]
  all.genus <- tax.tab.num[,"Genus"]
  all.species <- tax.tab.num[,"Species"]
  
  missing <- list()
  connection.with.upper <- list()
  connection <- character()
  
  tax.tab.num[,"Phylum"] <- as.character(lapply(tax.tab.num[,"Phylum"], function(x){ str_replace_all(x, "[[:punct:]]| |=", "")}))
  tax.tab.num[,"Class"] <- as.character(lapply(tax.tab.num[,"Class"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Order"] <- as.character(lapply(tax.tab.num[,"Order"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Family"] <- as.character(lapply(tax.tab.num[,"Family"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Genus"] <- as.character(lapply(tax.tab.num[,"Genus"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Species"] <- as.character(lapply(tax.tab.num[,"Species"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  
  Phylum <- names(table(tax.tab.num[,"Phylum"]))
  Phylum.desc <- c()
  Phylum.desc.dum <- c()
  
  sig.Phylum <- Phylum[!grepl("Phylum", Phylum)]
  ind <- as.numeric(sig.Phylum)+1
  
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)
  
  pos.sig.Phylum <- sig.Phylum[ind.pos.sig] 
  neg.sig.Phylum <- sig.Phylum[ind.neg.sig]
  na.Phylum <- Phylum[grepl("NA", Phylum)]
  non.sig.Phylum <- Phylum[!is.element(Phylum,c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum))]
  
  dummy <- c()
  for(i in 1:sum(!is.na(Phylum))){
    dummy <- c(dummy, paste0("dum",i))
  }
  
  phylum.wdum1 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  phylum.wdum2 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  
  for(i in seq(1,length(Phylum),2)){
    phylum.wdum1[i] <- dummy[i]
  }
  for(i in seq(2,length(Phylum),2)){
    phylum.wdum2[i] <- dummy[i]
  }
  
  phylum.wdum <- rbind(phylum.wdum1, phylum.wdum2)
  k <- 1
  for (i in 1:length(Phylum)) {                                 # just include the following for loops inside this for loop and save each graph as pdf
    ind <- which(tax.tab.num[,"Phylum"] == Phylum[i])           # for (i in 1: length(phylum)) {~~~~, for (i in 1:length(Class))~~~~, for (i in 1: length(Order)~~....)}
    if (sum(!is.na(tax.tab.num[ind,"Class"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Class"]))
      desc <- str_replace_all(desc, " ", "")
      Phylum.desc[i] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")
      
      if(Phylum[i] %in% phylum.wdum2){
        Phylum.desc.dum[k] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }else{
        Phylum.desc.dum[k] <- paste(phylum.wdum2[which(phylum.wdum1==Phylum[i])], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }
      connection <- c(connection, desc)
    }
    connection.with.upper[[1]] <- connection
  }
  
  Phylum.desc <- Phylum.desc[!is.na(Phylum.desc)]
  
  Bacteria.desc.dum <- paste("Bacteria -> {", paste(phylum.wdum1, collapse = " "), "}", sep = "", collapse = "")
  phy.to.phy <- paste(phylum.wdum1, " -> {", phylum.wdum2, "} ", sep = "")
  
  Phylum.desc <- paste(Phylum.desc, collapse = "; ")
  
  connection <- character()
  Class <- names(table(tax.tab.num[,"Class"]))
  Class.desc <- c()
  for (i in 1:length(Class)) {
    ind <- which(tax.tab.num[,"Class"] == Class[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Order"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Order"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Class.desc[i] <- paste(Class[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[2]] <- connection
  }
  
  Class.desc <- Class.desc[!is.na(Class.desc)]
  Class.desc <- paste(Class.desc, collapse = "; ")
  
  connection <- character()
  Order <- names(table(tax.tab.num[,"Order"]))
  Order <- Order[!Order %in% "NA"]
  Order.desc <- c()
  for (i in 1:length(Order)) {
    ind <- which(tax.tab.num[,"Order"] == Order[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Family"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Family"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Order.desc[i] <- paste(Order[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[3]] <- connection
  }
  Order.desc <- Order.desc[!is.na(Order.desc)]
  Order.desc <- paste(Order.desc, collapse = "; ")
  
  connection <- character()
  Family <- names(table(tax.tab.num[,"Family"]))
  Family <- Family[!Family %in% "NA"]
  Family.desc <- c()
  for (i in 1:length(Family)) {
    ind <- which(tax.tab.num[,"Family"] == Family[i])
    if (sum(!is.na(tax.tab.num[ind,"Genus"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Genus"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Family.desc[i] <- paste(Family[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[4]] <- connection
  }
  Family.desc <- Family.desc[!Family.desc %in% "NA"]
  Family.desc <- paste(Family.desc, collapse = "; ")
  
  connection <- character()
  Genus <- names(table(tax.tab.num[,"Genus"]))
  Genus <- Genus[!Genus %in% "NA"]
  Genus.desc <- c()
  
  if(length(setdiff(Genus,connection.with.upper[[4]])) > 0){
    for (i in 1:length(Genus)) {
      if(Genus[i] != setdiff(Genus, connection.with.upper[[4]])){
        ind <- which(tax.tab.num[,"Genus"] == Genus[i])
        
        if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
          desc <- names(table(tax.tab.num[ind,"Species"]))
          desc <- str_replace_all(desc, " ", "")
          desc <- desc[!(desc %in% "NA")]
          Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
          connection <- c(connection, desc)
        }
        connection.with.upper[[5]] <- connection 
      }
    }
  } else if(length(setdiff(Genus,connection.with.upper[[4]])) == 0){
    for (i in 1:length(Genus)) {
      ind <- which(tax.tab.num[,"Genus"] == Genus[i])
      
      if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
        desc <- names(table(tax.tab.num[ind,"Species"]))
        desc <- str_replace_all(desc, " ", "")
        desc <- desc[!(desc %in% "NA")]
        Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
        connection <- c(connection, desc)
      }
      connection.with.upper[[5]] <- connection 
    }
  }
  Genus.desc <- Genus.desc[!is.na(Genus.desc)]
  Genus.desc <- paste(Genus.desc, collapse = "; ")
  
  #Species <- names(table(tax.tab.num[,"Species"]))
  Class <- connection.with.upper[[1]]
  Order <- connection.with.upper[[2]]
  Family <- connection.with.upper[[3]]
  Genus <- connection.with.upper[[4]]
  Species <- connection.with.upper[[5]]
  
  missing[[1]] <- NA
  missing[[2]] <- setdiff(all.class, connection.with.upper[[1]])
  missing[[3]] <- setdiff(all.order, connection.with.upper[[2]])
  missing[[4]] <- setdiff(all.family, connection.with.upper[[3]])
  missing[[5]] <- setdiff(all.genus, connection.with.upper[[4]])
  missing[[6]] <- setdiff(all.species, connection.with.upper[[5]])
  
  Class <- Class[!(Class %in% "NA")]
  sig.Class <- Class[!grepl("Class", Class)]
  #ind <- text.tab.all[,1] %in% Class[!grepl("Class", Class)]
  ind <- as.numeric(sig.Class)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Class <- sig.Class[ind.pos.sig] 
  neg.sig.Class <- sig.Class[ind.neg.sig]
  
  Order <- Order[!(Order %in% "NA")]
  sig.Order <- Order[!grepl("Order", Order)]
  #ind <- text.tab.all[,1] %in% Order[!grepl("Order", Order)]
  ind <- as.numeric(sig.Order)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Order <- sig.Order[ind.pos.sig] 
  neg.sig.Order <- sig.Order[ind.neg.sig]
  
  Family <- Family[!(Family %in% "NA")]
  sig.Family <- Family[!grepl("Family", Family)]
  na.Family <- Family[grepl("NA", Family)]
  #ind <- text.tab.all[,1] %in% Family[!grepl("Family", Family)]
  ind <- as.numeric(sig.Family)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Family <- as.character(sig.Family[ind.pos.sig])
  neg.sig.Family <- as.character(sig.Family[ind.neg.sig])
  
  Genus <- Genus[!(Genus %in% "NA")]
  sig.Genus <- Genus[!grepl("Genus", Genus)]
  ind <- as.numeric(sig.Genus)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Genus <- sig.Genus[ind.pos.sig] 
  neg.sig.Genus <- sig.Genus[ind.neg.sig]
  
  Species <- Species[!(Species %in% "NA")]
  sig.Species <- Species[!grepl("Species", Species)]
  ind <- as.numeric(sig.Species)+1
  ind.pos.sig <- which(ci.tab.all[ind[!is.na(ind)],1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind[!is.na(ind)],1] < 0)  
  
  pos.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  neg.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  
  pos.sig.Species <- pos.sig.Species[ind.pos.sig] 
  neg.sig.Species <- neg.sig.Species[ind.neg.sig]
  
  if (species.include) {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, Genus.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus, pos.sig.Species)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus, neg.sig.Species)
    total.group <- c(Class, Order, Family, Genus, Species)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
    
  } else {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus)
    total.group <- c(Class, Order, Family, Genus)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
  }
  Taxa.desc <- str_remove_all(Taxa.desc, " NA;")
  
  flow.text.dum <- list()
  length(flow.text.dum) <- 6
  
  neg <- 1
  non <- 1
  na <- 1
  for(i in 1:length(phylum.wdum1)){
    if(phylum.wdum1[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum2[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum1[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum2[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum1[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }else if(phylum.wdum2[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }
    else if(phylum.wdum1[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }else if(phylum.wdum2[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }
  }
  
  if(length(pos.sig.taxa) >0){
    for(i in 1:length(pos.sig.taxa)){
      flow.text.dum[[2]][i] <- paste(pos.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = indianred, label = ",pos.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(neg.sig.taxa) >0){
    for(i in 1:length(neg.sig.taxa)){
      flow.text.dum[[3]][i] <- paste(neg.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = steelblue, label = ",neg.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(non.sig.taxa) >0){
    for(i in 1:length(non.sig.taxa)){
      flow.text.dum[[4]][i] <- paste(non.sig.taxa[i], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
    }
  }
  
  if(length(na.taxa) >0){
    for(i in 1:length(na.taxa)){
      flow.text.dum[[5]][i] <- paste(na.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
    }
  }
  
  if(length(dummy) >0){
    for(i in 1:length(dummy)){
      flow.text.dum[[6]][i] <- paste(dummy[i], paste("[fontname = Helvetica, shape = point, width = 0, style = invi, label = '']", sep = " "))
    }
  }
  
  flow.text.dum.complete <- c()
  for(i in 1:6){
    flow.text.dum.complete <- paste(flow.text.dum.complete, flow.text.dum[[i]], sep = "\n ")
  }
  
  flow.text <- paste("digraph ", layout.type, " {",
                     paste("graph [layout = ",layout.type,", root = Bacteria, ranksep = \"2 0.3 1 1 1 \"]", collapse = ""),
                     "Bacteria [fontname = Helvetica, fontsize = 12, fixedsize = true, width = 1, style = filled, fillcolor = lightblue, label = Bacteria]",
                     paste(flow.text.dum.complete, collapse = ""),
                     "edge [color = black, arrowhead = none]",
                     Taxa.desc,
                     "}", sep = "\n")
  
  taxon.tab <- data.frame("ID" = text.tab.all[-1,1],
                          #"Status" = ifelse( ci.tab.all[-1] < 0, "neg", "pos"), # "blue", "red"),
                          "Taxon" = sub(".*;", "", text.tab.all[-1,3]))
  
  return(list(flow.text = grViz(flow.text), taxon.tab = taxon.tab, ci.tab.all = ci.tab.all) )
}
  
  
taxa.cox.test <- function(survtime, censoring, taxa, cov.dat = NULL, multi.method = "BH") {
  
  censoring <- censoring
  cens.var.ord <- unique(censoring[[1]])[order(unique(censoring[[1]]))]
  print(cens.var.ord)
  
  if( !is.numeric(cens.var.ord) ){
    censoring[ censoring == cens.var.ord[1] ] <- 0
    censoring[ censoring == cens.var.ord[2] ] <- 1
    censoring <- as.numeric(censoring)
  }
  else if( !(cens.var.ord[1] == 0 || cens.var.ord[1] == 1)  ) {
    if( !(cens.var.ord[2] == 0 || cens.var.ord[2] == 1)  ){
      censoring[ censoring == cens.var.ord[1] ] <- 0
      censoring[ censoring == cens.var.ord[2] ] <- 1
      censoring <- as.numeric(censoring)
    }
  }
  
  cox.test <- list()
  for(i in 1:6) {
    
    taxa[[i]] <- taxa[[i]][!is.na(match(rownames(taxa[[i]]), rownames(survtime))),]
    survtime1 <- survtime[match(rownames(survtime), rownames(taxa[[i]]))]
    censoring1 <- censoring[match(rownames(survtime), rownames(taxa[[i]]))]
    
    n.tax <- ncol(taxa[[i]])
    cox.out <- matrix(NA, n.tax, 6)
    for (j in 1:n.tax) {
      taxon <- taxa[[i]][,j]
      if (is.null(cov.dat)){
        fit <- coxph(Surv(survtime1[[1]], censoring1[[1]]) ~  scale(taxon))#, silent = TRUE)
      } else {
        dat <- as.data.frame(cbind(survtime, censoring[[1]], taxon, cov.dat))
        f <- formula(paste("Surv(survtime1[[1]], censoring1[[1]]) ~  scale(taxon) +",  paste(colnames(dat)[4:ncol(dat)], collapse = "+")))
        fit <- try(coxph(f, data = dat), silent = TRUE)
      }
      
      
      if (class(fit) == "try-error"){
        cox.out[j,] <- c(rep(NA, 6))
      } else {
        fit.info <- summary(fit)
        
        cox.out[j,] <- c( round(fit.info$coefficients[1,1], 3), 
                          round(fit.info$coefficients[1,3], 3), 
                          round(fit.info$conf.int[1,1]    , 3),            #round(exp(fit$coefficients), 3)[1],
                          round(fit.info$conf.int[1,3]    , 3),
                          round(fit.info$conf.int[1,4]    , 3),
                          round(fit.info$coefficients[1,5], 3))
      }
    }
    
    cox.out <- as.data.frame(cox.out)
    rownames(cox.out) <- colnames(taxa[[i]])
    colnames(cox.out) <- c("Coef", "SE", "HR", "Lower", "Upper", "P.value")
    cox.out <- bin.q.func(cox.out, method = multi.method)
    cox.test[[i]] <- cox.out
  }
  names(cox.test) <- names(taxa)
  return(cox.test)
}

surv.taxa.forest.plot.pages <- function(all.taxa.q.out, species.include, mult.test.cor = "TRUE") {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out$Q.value < 0.05)  
    } else {
      ind.sig <- which(out$P.value < 0.05)
    }
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    plot.per.page <- total
    num.pages <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
    }
  } else if(total<=80 & total>=0) {
    num.pages = 1
    plot.per.page <-total
  }
  
  return(num.pages)
}

surv.taxa.forest.plot.pages1 <- function(all.taxa.q.out, taxa.names.out, species.include, report.type = "Coef", mult.test.cor = "TRUE") {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  if(report.type == "Est"){
    report.txt <- "Est."
  } else if(report.type == "OR") {
    report.txt <- "OR"
  } else if(report.type == "HR") {
    report.txt <- c("Coef", "HR")
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out$Q.value < 0.05)  
    } else {
      ind.sig <- which(out$P.value < 0.05)
    }
    
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    capacity <- c(40:80)
    plot.per.page <- total
    num.pages <- 0
    num.mod <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
      num.mod <- mod(total, num.pages)
    }
  } else if(total<=80 & total>=1) {
    num.pages = 1
    plot.per.page <-total
    num.mod <- 0
  }
  
  all.text.tab <- list()
  all.ci.tab <- list()
  if(total >0) {
    plot.taxa <- list()
    rank <- 1
    tax.rank <- c("Phylum","Class","Order","Family","Genus","Species")
    current.page <- 1
    text.tab.all <- 0
    text.tab.save <- numeric()
    ci.tab.save <- numeric()
    info.taxa <- numeric()
    info.data <- numeric()
    info.rank <- numeric()
    info.ci <- numeric()
    sig.out <- numeric()
    actual.plot.per.page <- plot.per.page
    
    if (text.tab.all == 0) {
      if (mult.test.cor){
        text.tab.all <- matrix(c("ID", "Rank", "Taxon", report.txt, "P-value", "Q-value"), 1, (5+length(report.txt)))  
      } else {
        text.tab.all <- matrix(c("ID", "Rank", "Taxon", report.txt, "P-value"), 1, (4+length(report.txt)))
      }
      ci.tab.all <- matrix(c(NA, NA, NA), 1, 3)
    }
    
    for (rank in 1:range) {
      if(!is.null(all.taxa.q.out[[rank]])){
        ind.sig <- sig.by.rank[[rank]]
        sig.out <- all.taxa.q.out[[rank]][sig.by.rank[[rank]],]
        
        info.taxa <- c(info.taxa, taxa.names.out$names[[rank]][sig.by.rank[[rank]]])
        info.data <- rbind(info.data, sig.out)
        info.ci <- rbind(info.ci, cbind(sig.out[,report.type], sig.out[,"Lower"], sig.out[,"Upper"]))
        info.rank <- c(info.rank, rep(tax.rank[rank], length(sig.by.rank[[rank]])))
      }
    }
    
    if(mult.test.cor){
      info <- as.matrix(cbind(c(1:nrow(info.data)), info.rank, info.taxa, format(round(info.data[,c(report.txt, "P.value", "Q.value")],3), nsmall = 3)))
    } else {
      info <- as.matrix(cbind(c(1:nrow(info.data)), info.rank, info.taxa, format(round(info.data[,c(report.txt, "P.value")],3), nsmall = 3)))
    }
    
    
    for(p in 1:num.pages) {
      if(p <= num.mod){
        actual.plot.per.page <- plot.per.page+1
        initial <- 1+actual.plot.per.page*(p-1)
      } else {
        initial <- 1+actual.plot.per.page*(p-1)
        actual.plot.per.page <- plot.per.page
      }
      names(info) <- names(text.tab.all)
      all.text.tab[[p]] <- rbind(as.matrix(text.tab.all), info[c(initial:(initial+actual.plot.per.page-1)),])
      all.ci.tab[[p]] <- rbind(ci.tab.all, info.ci[c(initial:(initial+actual.plot.per.page-1)),])
    }
  } else {
    all.text.tab <- NULL
    all.ci.tab <- NULL
  }
  return(list(all.text.tab = all.text.tab, all.ci.tab = all.ci.tab))
}


duplicate.surv <- function(duplicate.taxa, taxon.inplot, duplicate.full.list){
  print("surv:success in duplicates if first here")
  if(length(duplicate.taxa %in% taxon.inplot)>0) {
    duplicate.taxa <- unlist(duplicate.full.list)[duplicate.taxa %in% taxon.inplot]
    par(mar=c(0, 0.5, 0, 0.5))
    text(x=0, y=0.5, paste(duplicate.taxa, collapse = "\n"), cex = 0.75, adj = c(0, NA))
  } else {
    text(x=0.5, y=0.5, "")
  }
}


survival.pages2 <- function(page.taxa.q.out, page) {
  
  print("inside survival forest plot: success")
  text.tab.all <- page.taxa.q.out$all.text.tab[[page]]
  ci.tab.all <- page.taxa.q.out$all.ci.tab[[page]]
  
  if(is.null(text.tab.all) & is.null(ci.tab.all)){
    plot.new()
    text(x = 0.5, y = 0.5, "No significant taxa are found.", 
         cex = 1.2, col = "black")
  }else{
    str.max <- lapply(page.taxa.q.out$all.text.tab, nchar)
    for(i in 1:length(page.taxa.q.out$all.text.tab)){
      str.max[[i]] <- str.max[[i]][,3]
    }
    maxStr <- max(unlist(str.max))
    if(!is.numeric(maxStr)){
      maxStr <- 0
    }
    
    text.tab.all[,3] <- substr(text.tab.all[,3], 1, 55)
    text.tab.all[,6] <- p.value.0.1_char(text.tab.all[,6]) 
    text.tab.all[,7] <- p.value.0.1_char(text.tab.all[,7])
    
    #par(mar=c(0, 0.2, 0, 0.2))
    if(text.tab.all[1,5] == "Coef"){
      if(nrow(ci.tab.all) <= 5) {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), lineheight = unit(1, "cm"), #line.margin = unit(0.2, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
      else if(nrow(ci.tab.all) <= 45){
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      } else {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), 
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
      
    } else {
      if(nrow(ci.tab.all) <= 5){
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), lineheight = unit(1, "cm"), #line.margin = unit(0.08, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      } else if(nrow(ci.tab.all) <= 45){
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      } else {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), 
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
    } 
  }
}


# 
# 
# taxa.surv.forest.plot.pages2 <- function(page.taxa.q.out, page) {
#   
#   text.tab.all <- page.taxa.q.out$all.text.tab[[page]]
#   ci.tab.all <- page.taxa.q.out$all.ci.tab[[page]]
#   
#   if(is.null(text.tab.all) & is.null(ci.tab.all)){
#     plot.new()
#     text(x = 0.5, y = 0.5, "No significant taxa are found.", 
#          cex = 1.2, col = "black")
#   }else{
#     str.max <- lapply(page.taxa.q.out$all.text.tab, nchar)
#     for(i in 1:length(page.taxa.q.out$all.text.tab)){
#       str.max[[i]] <- str.max[[i]][,3]
#     }
#     maxStr <- max(unlist(str.max))
#     if(!is.numeric(maxStr)){
#       maxStr <- 0
#     }
#     
#     text.tab.all[,3] <- substr(text.tab.all[,3], 1, 55)
# 
#     if(text.tab.all[1,5] == "Coef"){
#       if(nrow(ci.tab.all) <= 5) {
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = unit(1, "cm"), #line.margin = unit(0.2, "cm"),
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       }
#       else if(nrow(ci.tab.all) <= 45){
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       } else {
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), 
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       }
#       
#     } else {
#       if(nrow(ci.tab.all) <= 5){
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = unit(1, "cm"), #line.margin = unit(0.08, "cm"),
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       } else if(nrow(ci.tab.all) <= 45){
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       } else {
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), 
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       }
#     } 
#   }
# }
# 
# 
