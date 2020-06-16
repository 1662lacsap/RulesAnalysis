#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("genefilter", version = "3.8")

library (data.table)
library(tableHTML)
library(DT)
library(dplyr)
library(shinyjs)
library(dmdScheme)
library(shinydashboard)
library(formattable)
library(xlsx)
library(readxl)
library(DDoutlier)

library(formattable)

library(stringr)


library(DMwR)
library(lattice)
library(grid)
library(Rlof)
library(doParallel)
library(foreach)
library(iterators)
library(parallel)
library(cluster)

library(StatMatch)
library(clusterSim)
library(clValid)


library(stats)
library(factoextra)

server <- function(input, output, session) {
  
  # overwriting R's version (not advised)
  `<=` <- function(x, y) (isTRUE(all.equal(x, y)) | (x < y))
  `>=` <- function(x, y) (isTRUE(all.equal(x, y)) | (x > y))
 
  Parametry<- c("Indeks sylwetki","Indeks Dunna","Indeks Daviesa-Bouldina","Pseudo F", "Indeks Huberta", "Miara Gamma", "Liczba grup", "k-sąsiadów lub set.seed", "% reguł nietypowych jako liczba", "Współczynnik CPCC", "Metoda łączenia", "Dataset / Reguły nietypowe" )
  Z_Regułami_Nietypowymi <- 1:12
  Bez_Reguł_Nietypowych <- 1:12
  wyniki<-data.frame(Parametry,Z_Regułami_Nietypowymi,Bez_Reguł_Nietypowych)
  
  
  output$contents  <- renderTable({
     
    req(input$file1)
    
    tryCatch(
      {
        collection<- fread(input$file1$datapath,
                            header = input$header,
                            sep = input$sep)
        file1 = input$file1
        saveRDS(file1 , file = "file1.rds")
        quote2 = input$quote2                   
        saveRDS(quote2, file = "quote2.rds")
        
        ATTRIBUTES <- str_sub(string=collection[2], start=12)
        ATTRIBUTES1<- as.numeric(ATTRIBUTES)
        DECISION_VALUES<- str_sub(string=collection[ATTRIBUTES1+3], start=17)
        DECISION_VALUES1<-as.numeric(DECISION_VALUES)
        saveRDS(DECISION_VALUES1, file = "RULES1.rds")
        RULES<- str_sub(string=collection[ATTRIBUTES1+3+DECISION_VALUES1+1], start=7)
        RULES1<-as.numeric(RULES)
        REGULY_START<-(ATTRIBUTES1+3+DECISION_VALUES1+1+1)
        REGULY_KONIEC<-(REGULY_START+RULES1-1)
        collection1<-collection[(REGULY_START:REGULY_KONIEC),]
        collection2<-str_split(string=collection1$V1, pattern = "=>")
        collection3<-data.table(Reduce(rbind, collection2))
        collection4<-collection3[,1]
        m <- as.data.frame(matrix(0, ncol = ATTRIBUTES1, nrow = RULES1))
        q <- ATTRIBUTES1+2
        k <- c(3:q)
        l<-c(1:RULES1)
        
        for (j in l){
          for (i in k){
            
            collection_3 = unlist(gregexpr(pattern =' ',collection[i]))
            collection_4<-str_sub(string=collection[i], start = 1, end = collection_3[1]-1)
            collection5 <- str_extract(string=collection4[j], pattern = collection_4)
            
            if(is.na(collection5) == TRUE ||  length(unlist(gregexpr(pattern=collection5, collection4[j])))>1) {
              m[j,i-2] = input$quote } else {   
                collection6=unlist(gregexpr(pattern=collection5, collection4[j]))
                collection7 = nchar(collection5)+1
                collection_start<- str_sub(string = collection4[j], start = collection6)
                collection7a<-gregexpr(pattern =')', collection_start)
                ostatni_skladnik = collection7a[[1]][1]
                ostatni_skladnik1=ostatni_skladnik-2
                collection8<-str_sub(string = collection4[j], start = (collection6 + collection7), end =  collection6 + ostatni_skladnik1)
                collection9=unlist(gregexpr(pattern="=",  collection8))
                if (collection9==-1)
                {m[j,i-2] = collection8} else { 
                  m[j,i-2] = input$quote
                }
              }
          }
        }
        
    
        p<-(ATTRIBUTES1-1)
        d <- m[ ,1:p]
        
       write.csv(d,"wczytanyZbior.csv", row.names = FALSE)
               
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    if(input$disp == "head") {
      return(head(d))
    }
    else {
      return(d)
    }
    
  }) 
    
  output$wyjscieTabela <-render_tableHTML({ 
      
   qoute2 <- readRDS(file = "quote2.rds")
    
    g <- fread ("wczytanyZbior.csv", sep=',', header = TRUE, stringsAsFactors = FALSE, quote = qoute2)
    g[g=="NA"] <- NA  
    gower_mat1 <- gower.dist (g ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL)
    gower_mat1[is.nan(gower_mat1)] <- 0
    gower_mat <-as.dist(gower_mat1) 
    clusters <- hclust(gower_mat, method = input$metoda)    
    wyniki[11,2] <- input$metoda
    
    #Cophenetic Distances for a Hierarchical Clustering
    d1 <- gower_mat
    hc <- clusters
    d2 <- cophenetic(hc)
    wyniki[10,2] <- cor(d1, d2)
    
    # Ocena jakosci grupowania - cluster validity 
     cutree_k <- input$liczbaGrup[1]
     wyniki[7,2]<-cutree_k
     wyniki[7,3]<-cutree_k
     cutree_k_LOF <- input$liczbaGrup[1]
     cutree_k_COF <- input$liczbaGrup[1]
     cutree_k_kmeans <- input$liczbaGrup[1]
     cutree_k_SMALLCLUSTER <- input$liczbaGrup[1]
    
    #cutree(clusters, k = cutree_k) - wektor liczb calkowitych dodatnich informujacy 
    #o przynaleznosci obiektow do klasy
    cutree_wektor <- cutree(clusters, k = cutree_k) # k liczba klas
    
    #1 Oparte o odleglosc miedzy jednostkami w skupieniu i miedzy skupieniami 
    
     si3 <- silhouette(cutree_wektor, gower_mat)
     icq_sylwetki <- mean(si3[,"sil_width"])
      wyniki[1,2]<-icq_sylwetki
    
     wyniki[2,2]<-dunn(gower_mat, cutree_wektor)
    
    #2 Oparte o rozproszenie jednostek w skupieniu i odleglosci miedzy skupieniami
    
    index_DB <- index.DB(gower_mat, cutree_wektor, centrotypes="centroids", p = 2, q = 2)
    wyniki[3,2]<-index_DB[1]

    #3 Oparte na sumie kwadratow wewnatrz skupien i miedzy skupieniami
    
    icq_index_G1 <- index.G1(gower_mat1, cutree_wektor, d = NULL, centrotypes = "centroids")
    wyniki[4,2]<-icq_index_G1
   
    icq_index_G3 <- index.G3(gower_mat, cutree_wektor)
    wyniki[5,2]<-icq_index_G3
    
    icq_index_G2 <- index.G2(gower_mat, cutree_wektor)
    wyniki[6,2]<-icq_index_G2
   
    wyniki[12,2] <-readRDS(file = "file1.rds") 
    

    #LOF
    
    if (input$algorytm == "LOF")
    {
   
      qoute2 <- readRDS(file = "quote2.rds")
      
      e <- fread ("wczytanyZbior.csv", sep=',', header = TRUE,stringsAsFactors = FALSE, quote = qoute2)
      e[e=="NA"] <- NA
     
     
      RULES1<-nrow(e)
      howmany <-(input$procentRegulNietypowych[1]/100)*RULES1
      howmany<-ceiling(howmany)
      wyniki[9,2]<-0
      wyniki[9,3]<-howmany
      set_k <- input$kSasiadow[1] 
      wyniki[8,2]<-NA
      wyniki[8,3]<-set_k
      
      outlier.scores <- lofactor(gower_mat1, k = set_k)
      outliers_LOF <- order(outlier.scores, decreasing = T) [1:howmany]
      wyniki[12,3] <- toString(outliers_LOF)
    
      p <- outliers_LOF
      data_inf1 <- e[-p, ]
      
      # dendrogram AHC klastrowanie hierarchiczne po usunieciu odchylen - outliers_LOF
      gower_mat1 <- gower.dist ( data_inf1 ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL )
      gower_mat1[is.nan(gower_mat1)] <- 0
      
      gower_mat <-as.dist(gower_mat1) 
       
      clusters <- hclust(gower_mat,method = input$metoda) 
     
      wyniki[11,3] <- input$metoda
      
      #Cophenetic Distances for a Hierarchical Clustering
      #require(graphics)
      d1 <- gower_mat
      hc <- clusters
      d2 <- cophenetic(hc)
      wyniki[10,3] <- cor(d1, d2)
      cutree_wektor <- cutree(clusters, k = cutree_k_LOF)
      si3 <- silhouette(cutree_wektor, gower_mat) 
      icq1 <- mean(si3[,"sil_width"])
      wyniki[1,3]<-icq1
      wyniki[2,3]<-dunn(gower_mat, cutree_wektor)
      index_DB <- index.DB(gower_mat, cutree_wektor, centrotypes="centroids", p = 2, q = 2)
      wyniki[3,3]<-index_DB[1]
      icq_index_G1 <- index.G1(gower_mat1, cutree_wektor, d = NULL, centrotypes = "centroids")
      wyniki[4,3]<-icq_index_G1
      icq_index_G3 <- index.G3(gower_mat, cutree_wektor)
      wyniki[5,3]<-icq_index_G3
      icq_index_G2 <- index.G2(gower_mat, cutree_wektor)
      wyniki[6,3]<-icq_index_G2
      
    }
    
    #COF

    if (input$algorytm == "COF")
    {
      
      qoute2 <- readRDS(file = "quote2.rds")
      f <- fread ("wczytanyZbior.csv", sep=',', header = TRUE,stringsAsFactors = FALSE, quote = qoute2)
      f[f=="NA"] <- NA
      RULES1<-nrow(f)
      howmany <-(input$procentRegulNietypowych[1]/100)*RULES1
      howmany<-ceiling(howmany)
      wyniki[9,2]<-0
      wyniki[9,3]<-howmany
      set_k <- input$kSasiadow[1] 
      wyniki[8,2]<-NA
      wyniki[8,3]<-set_k
      outlier.scores <- COF(gower_mat1, k = set_k)
      outliers_COF <- order(outlier.scores, decreasing = T) [1:howmany]
      wyniki[12,3] <- toString(outliers_COF)
      q <- outliers_COF
      data_inf2 <- f[-q, ]
      gower_mat1 <- gower.dist ( data_inf2 ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL )
      gower_mat1[is.nan(gower_mat1)] <- 0
      gower_mat <-as.dist(gower_mat1)  
      clusters <- hclust(gower_mat,method = input$metoda) 
      wyniki[11,3] <- input$metoda
      d1 <- gower_mat
      hc <- clusters
      d2 <- cophenetic(hc)
      wyniki[10,3] <- cor(d1, d2)
      cutree_wektor <- cutree(clusters, k = cutree_k_COF)
      si3 <- silhouette(cutree_wektor, gower_mat) 
      icq2 <- mean(si3[,"sil_width"])
      wyniki[1,3]<-icq2
      wyniki[2,3]<-dunn(gower_mat, cutree_wektor)
      index_DB <- index.DB(gower_mat, cutree_wektor, centrotypes="centroids", p = 2, q = 2)
      wyniki[3,3]<-index_DB[1]
      icq_index_G1 <- index.G1(gower_mat1, cutree_wektor, d = NULL, centrotypes = "centroids")
      wyniki[4,3]<-icq_index_G1
      icq_index_G3 <- index.G3(gower_mat, cutree_wektor)
      wyniki[5,3]<-icq_index_G3
      icq_index_G2 <- index.G2(gower_mat, cutree_wektor)
      wyniki[6,3]<-icq_index_G2
      

      # KMEANS
      
    }
      
    if (input$algorytm == "KMEANS")
      
    {
      qoute2 <- readRDS(file = "quote2.rds")
      h <- fread ("wczytanyZbior.csv", sep=',', header = TRUE,stringsAsFactors = FALSE, quote = qoute2)
      h[h=="NA"] <- NA
      RULES1<-nrow(h)
      howmany <-(input$procentRegulNietypowych[1]/100)*RULES1
      howmany<-ceiling(howmany)
      wyniki[9,2]<-0
      wyniki[9,3]<-howmany
      set_k <- input$setSeed[1] 
      wyniki[8,2]<-NA
      wyniki[8,3]<-set_k
      RNGkind(sample.kind = "Rounding")
      set.seed(set_k)
      DECISION_VALUES1<-readRDS(file = "RULES1.rds")
      kmeans.result <- kmeans(gower_mat1 , centers = DECISION_VALUES1)
      centers <- kmeans.result$centers[kmeans.result$cluster, ]
      distances <- sqrt(rowSums((as.numeric(unlist(gower_mat1)) - centers)^2))
      outliers_kmeans <- order( distances, decreasing = T ) [1:howmany]
      wyniki[12,3] <- toString(outliers_kmeans)
      r <- outliers_kmeans
      data_inf3 <- h[-r, ]
      gower_mat1 <- gower.dist ( data_inf3 ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL )
      gower_mat1[is.nan(gower_mat1)] <- 0
      gower_mat <-as.dist(gower_mat1) 
      clusters <- hclust(gower_mat, method =  input$metoda) 
      wyniki[11,3] <-  input$metoda
      d1 <- gower_mat
      hc <- clusters
      d2 <- cophenetic(hc) 
      wyniki[10,3] <- cor(d1, d2)
      cutree_wektor <- cutree(clusters, k = cutree_k_kmeans)
      si3 <- silhouette(cutree_wektor, gower_mat) 
      icq3 <- mean(si3[,"sil_width"])
      wyniki[1,3]<-icq3
      wyniki[2,3]<-dunn(gower_mat, cutree_wektor)
      index_DB<-index.DB(gower_mat, cutree_wektor, centrotypes="centroids", p = 2, q = 2)
      wyniki[3,3]<-index_DB[1]
      icq_index_G1 <- index.G1(gower_mat1, cutree_wektor, d = NULL, centrotypes = "centroids")
      wyniki[4,3]<-icq_index_G1
      icq_index_G3 <- index.G3(gower_mat, cutree_wektor)
      wyniki[5,3]<-icq_index_G3 
      icq_index_G2 <- index.G2(gower_mat, cutree_wektor)
      wyniki[6,3]<-icq_index_G2
      
    }


    # SMALLCLUSTER

    if (input$algorytm == "SMALLCLUSTER")
    {
     
      qoute2 <- readRDS(file = "quote2.rds")
      j <- fread ("wczytanyZbior.csv", sep=',', header = TRUE,stringsAsFactors = FALSE, quote = qoute2)
      j[j=="NA"] <- NA
      cutree_k_SMALLCLUSTER <- input$liczbaGrup[1]
      cutree_k_SMALLCLUSTER_outliers <- input$smallCluster[1]
      
      if (cutree_k_SMALLCLUSTER_outliers > nrow(j))
      { return("Wybrana liczba klastrów jest większa niż liczba obserwacji") }
      
      wyniki[9,2]<-0
      wyniki[9,3]<-"smallest clusters" 
      wyniki[8,2]<-NA
      wyniki[8,3]<-NA
      gower_mat1 <- gower.dist ( j ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL )
      gower_mat1[is.nan(gower_mat1)] <- 0
      gower_mat <-as.dist(gower_mat1) 
      clusters <- hclust(gower_mat,method = input$metoda)
      cutree_wektor <- cutree(clusters, k = cutree_k_SMALLCLUSTER_outliers)
      
      ilosc_obserwacji <-length(cutree_wektor) 
      ilosc_obserwacji_1 <-length(which(cutree_wektor == 1))
      ilosc_obserwacji_2 <-length(which(cutree_wektor == 2))
      ilosc_obserwacji_3 <-length(which(cutree_wektor == 3))
      ilosc_obserwacji_4 <-length(which(cutree_wektor == 4))
      ilosc_obserwacji_5 <-length(which(cutree_wektor == 5))
      ilosc_obserwacji_6 <-length(which(cutree_wektor == 6))
      ilosc_obserwacji_7 <-length(which(cutree_wektor == 7))
      ilosc_obserwacji_8 <-length(which(cutree_wektor == 8))
      ilosc_obserwacji_9 <-length(which(cutree_wektor == 9))
      ilosc_obserwacji_10 <-length(which(cutree_wektor == 10))
      
      
      if ((ilosc_obserwacji_1/ilosc_obserwacji) > 0 & (ilosc_obserwacji_1/ilosc_obserwacji) <= input$procentZbioru)
      {
        a1<-TRUE
      }
      else
      {
        a1<-FALSE
      }
      
      
      if ((ilosc_obserwacji_2/ilosc_obserwacji) > 0 & (ilosc_obserwacji_2/ilosc_obserwacji) <= input$procentZbioru)
      {
        b2<-TRUE
      }
      else
      {
        b2<-FALSE
      }
      
      if ((ilosc_obserwacji_3/ilosc_obserwacji) > 0 & (ilosc_obserwacji_3/ilosc_obserwacji) <= input$procentZbioru)
      {
        c3<-TRUE
      }
      else
      {
        c3<-FALSE
      }
      
      if ((ilosc_obserwacji_4/ilosc_obserwacji) > 0 & (ilosc_obserwacji_4/ilosc_obserwacji) <= input$procentZbioru)
      {
        d4<-TRUE
      }
      else
      {
        d4<-FALSE
      }
      
      if ((ilosc_obserwacji_5/ilosc_obserwacji) > 0 & (ilosc_obserwacji_5/ilosc_obserwacji) <= input$procentZbioru)
      {
        e5<-TRUE
      }
      else
      {
        e5<-FALSE
      }
      
      if ((ilosc_obserwacji_6/ilosc_obserwacji) > 0 & (ilosc_obserwacji_6/ilosc_obserwacji) <= input$procentZbioru)
      {
        f6<-TRUE
      }
      else
      {
        f6<-FALSE
      }
      
      if ((ilosc_obserwacji_7/ilosc_obserwacji) > 0 & (ilosc_obserwacji_7/ilosc_obserwacji) <= input$procentZbioru)
      {
        g7<-TRUE
      }
      else
      {
        g7<-FALSE
      }
      
      if ((ilosc_obserwacji_8/ilosc_obserwacji) > 0 & (ilosc_obserwacji_8/ilosc_obserwacji) <= input$procentZbioru)
      {
        h8<-TRUE
      }
      else
      {
        h8<-FALSE
      }
      
      if ((ilosc_obserwacji_9/ilosc_obserwacji) > 0 & (ilosc_obserwacji_9/ilosc_obserwacji) <= input$procentZbioru)
      {
        i9<-TRUE
      }
      else
      {
        i9<-FALSE
      }
      
      if ((ilosc_obserwacji_10/ilosc_obserwacji) > 0 & (ilosc_obserwacji_10/ilosc_obserwacji) <= input$procentZbioru)
      {
        j10<-TRUE
      }
      else
      {
        j10<-FALSE
      }
      
      
###########################################################################################################      
      
      #a1 b2 c3 d4
      
      if (a1 == TRUE & b2 == TRUE & c3 == TRUE & d4 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4)
      }
      
      #a1 b2 c3 e5
      
      else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & e5 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5)
      }
      
      
      #a1 b2 c3 f6
      
      else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6)
      }
      
      
      #a1 b2 c3 g7
      
      else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7)
      }
      
      
      #a1 b2 c3 h8
      
      else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8)
      }
      
      
      #a1 b2 c3 i9
      
      else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9)
      }
      
      
      #a1 b2 c3 j10
      
      else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_10)
      }
      
      #a1 b2 d4 e5
      
      else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & e5 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
      }
      
      #a1 b2 d4 f6
      
      else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
      }
      
      #a1 b2 d4 g7
      
      else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
      }
      
      #a1 b2 d4 h8
      
      else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
      }
      
      #a1 b2 d4 i9
      
      else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
      }
      
      #a1 b2 d4 j10
      
      else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
      }
      
      #a1 b2 e5 f6
      
      else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      #a1 b2 e5 g7
      
      else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      #a1 b2 e5 h8
      
      else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      #a1 b2 e5 i9
      
      else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      #a1 b2 e5 j10
      
      else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      #a1 b2 f6 g7
      
      else if (a1 == TRUE & b2 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #a1 b2 f6 h8
      
      else if (a1 == TRUE & b2 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #a1 b2 f6 i9
      
      else if (a1 == TRUE & b2 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #a1 b2 f6 j10
      
      else if (a1 == TRUE & b2 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #a1 b2 g7 h8
      
      else if (a1 == TRUE & b2 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #a1 b2 g7 i9
      
      else if (a1 == TRUE & b2 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #a1 b2 g7 j10
      
      else if (a1 == TRUE & b2 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #a1 b2 h8 i9
      
      else if (a1 == TRUE & b2 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #a1 b2 h8 j10
      
      else if (a1 == TRUE & b2 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #a1 b2 i9 j10
      
      else if (a1 == TRUE & b2 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #a1 c3 d4 e5
      
      else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & e5 == TRUE)
      {
        outliers_SMALLCLUSTER_1 <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
      }
      
      #a1 c3 d4 f6
      
      else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
      }
      
      #a1 c3 d4 g7
      
      else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
      }
      
      #a1 c3 d4 h8
      
      else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
      }
      
      #a1 c3 d4 i9
      
      else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
      }
      
      #a1 c3 d4 j10
      
      else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
      }
      
      #a1 c3 e5 f6
      
      else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      #a1 c3 e5 g7
      
      else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      #a1 c3 e5 h8
      
      else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      #a1 c3 e5 i9
      
      else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      #a1 c3 e5 j10
      
      else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      #a1 c3 f6 g7
      
      else if (a1 == TRUE & c3 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #a1 c3 f6 h8
      
      else if (a1 == TRUE & c3 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #a1 c3 f6 i9
      
      else if (a1 == TRUE & c3 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #a1 c3 f6 j10
      
      else if (a1 == TRUE & c3 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #a1 c3 g7 h8
      
      else if (a1 == TRUE & c3 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #a1 c3 g7 i9
      
      else if (a1 == TRUE & c3 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #a1 c3 g7 j10
      
      else if (a1 == TRUE & c3 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #a1 c3 h8 i9
      
      else if (a1 == TRUE & c3== TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #a1 c3 h8 j10
      
      else if (a1 == TRUE & c3 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #a1 c3 i9 j10
      
      else if (a1 == TRUE & c3 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #a1 d4 e5 f6
      
      else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      #a1 d4 e5 g7
      
      else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      #a1 d4 e5 h8
      
      else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      #a1 d4 e5 i9
      
      else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      #a1 d4 e5 j10
      
      else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      #a1 d4 f6 g7
      
      else if (a1 == TRUE & d4 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #a1 d4 f6 h8
      
      else if (a1 == TRUE & d4 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #a1 d4 f6 i9
      
      else if (a1 == TRUE & d4 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #a1 d4 f6 j10
      
      else if (a1 == TRUE & d4 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #a1 d4 g7 h8
      
      else if (a1 == TRUE & d4 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #a1 d4 g7 i9
      
      else if (a1 == TRUE & d4 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #a1 d4 g7 j10
      
      else if (a1 == TRUE & d4 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #a1 d4 h8 i9
      
      else if (a1 == TRUE & d4 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #a1 d4 h8 j10
      
      else if (a1 == TRUE & d4 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #a1 d4 i9 j10
      
      else if (a1 == TRUE & d4 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #a1 e5 f6 g7
      
      else if (a1 == TRUE & e5 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #a1 e5 f6 h8
      
      else if (a1 == TRUE & e5 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #a1 e5 f6 i9
      
      else if (a1 == TRUE & e5 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #a1 e5 f6 j10
      
      else if (a1 == TRUE & e5 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #a1 e5 g7 h8
      
      else if (a1 == TRUE & e5 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #a1 e5 g7 i9
      
      else if (a1 == TRUE & e5 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #a1 e5 g7 j10
      
      else if (a1 == TRUE & e5 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #a1 e5 h8 i9
      
      else if (a1 == TRUE & e5 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #a1 e5 h8 j10
      
      else if (a1 == TRUE & e5 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #a1 e5 i9 j10
      
      else if (a1 == TRUE & e5 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #a1 f6 g7 h8
      
      else if (a1 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #a1 f6 g7 i9
      
      else if (a1 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #a1 f6 g7 j10
      
      else if (a1 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #a1 f6 h8 i9
      
      else if (a1 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #a1 f6 h8 j10
      
      else if (a1 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #a1 f6 i9 j10
      
      else if (a1 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #a1 g7 h8 i9
      
      else if (a1 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #a1 g7 h8 j10
      
      else if (a1 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #a1 g7 i9 j10
      
      else if (a1 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #a1 h8 i9 j10
      
      else if (a1 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #b2 c3 d4 e5
      
      else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & e5 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
      }
      
      #b2 c3 d4 f6
      
      else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
      }
      
      #b2 c3 d4 g7
      
      else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
      }
      
      #b2 c3 d4 h8
      
      else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
      }
      
      #b2 c3 d4 i9
      
      else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
      }
      
      #b2 c3 d4 j10
      
      else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
      }
      
      #b2 c3 e5 f6
      
      else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      #b2 c3 e5 g7
      
      else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      #b2 c3 e5 h8
      
      else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      #b2 c3 e5 i9
      
      else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      #b2 c3 e5 j10
      
      else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      #b2 c3 f6 g7
      
      else if (b2 == TRUE & c3 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #b2 c3 f6 h8
      
      else if (b2 == TRUE & c3 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #b2 c3 f6 i9
      
      else if (b2 == TRUE & c3 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #b2 c3 f6 j10
      
      else if (b2 == TRUE & c3 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #b2 c3 g7 h8
      
      else if (b2 == TRUE & c3 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #b2 c3 g7 i9
      
      else if (b2 == TRUE & c3 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #b2 c3 g7 j10
      
      else if (b2 == TRUE & c3 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #b2 c3 h8 i9
      
      else if (b2 == TRUE & c3 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #b2 c3 h8 j10
      
      else if (b2 == TRUE & c3 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #b2 c3 i9 j10
      
      else if (b2 == TRUE & c3 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #b2 d4 e5 f6
      
      else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      #b2 d4 e5 g7
      
      else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      #b2 d4 e5 h8
      
      else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      #b2 d4 e5 i9
      
      else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      #b2 d4 e5 j10
      
      else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      #b2 d4 f6 g7
      
      else if (b2 == TRUE & d4 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #b2 d4 f6 h8
      
      else if (b2 == TRUE & d4 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #b2 d4 f6 i9
      
      else if (b2 == TRUE & d4 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #b2 d4 f6 j10
      
      else if (b2 == TRUE & d4 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #b2 d4 g7 h8
      
      else if (b2 == TRUE & d4 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #b2 d4 g7 i9
      
      else if (b2 == TRUE & d4 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #b2 d4 g7 j10
      
      else if (b2 == TRUE & d4 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #b2 d4 h8 i9
      
      else if (b2 == TRUE & d4 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #b2 d4 h8 j10
      
      else if (b2 == TRUE & d4 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #b2 d4 i9 j10
      
      else if (b2 == TRUE & d4 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #b2 e5 f6 g7
      
      else if (b2 == TRUE & e5 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #b2 e5 f6 h8
      
      else if (b2 == TRUE & e5 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #b2 e5 f6 i9
      
      else if (b2 == TRUE & e5 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #b2 e5 f6 j10
      
      else if (b2 == TRUE & e5 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #b2 e5 g7 h8
      
      else if (b2 == TRUE & e5 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #b2 e5 g7 i9
      
      else if (b2 == TRUE & e5 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #b2 e5 g7 j10
      
      else if (b2 == TRUE & e5 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #b2 e5 h8 i9
      
      else if (b2 == TRUE & e5 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #b2 e5 h8 j10
      
      else if (b2 == TRUE & e5 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #b2 e5 i9 j10
      
      else if (b2 == TRUE & e5 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #b2 f6 g7 h8
      
      else if (b2 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #b2 f6 g7 i9
      
      else if (b2 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #b2 f6 g7 j10
      
      else if (b2 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #b2 f6 h8 i9
      
      else if (b2 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #b2 f6 h8 j10
      
      else if (b2 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #b2 f6 i9 j10
      
      else if (b2 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #b2 g7 h8 i9
      
      else if (b2 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #b2 g7 h8 j10
      
      else if (b2 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #b2 g7 i9 j10
      
      else if (b2 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #b2 h8 i9 j10
      
      else if (b2 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #c3 d4 e5 f6
      
      else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      #c3 d4 e5 g7
      
      else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      #c3 d4 e5 h8
      
      else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      #c3 d4 e5 i9
      
      else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      #c3 d4 e5 j10
      
      else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      #c3 d4 f6 g7
      
      else if (c3 == TRUE & d4 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #c3 d4 f6 h8
      
      else if (c3 == TRUE & d4 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #c3 d4 f6 i9
      
      else if (c3 == TRUE & d4 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #c3 d4 f6 j10
      
      else if (c3 == TRUE & d4 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #c3 d4 g7 h8
      
      else if (c3 == TRUE & d4 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #c3 d4 g7 i9
      
      else if (c3 == TRUE & d4 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #c3 d4 g7 j10
      
      else if (c3 == TRUE & d4 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #c3 d4 h8 i9
      
      else if (c3 == TRUE & d4 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #c3 d4 h8 j10
      
      else if (c3 == TRUE & d4 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #c3 d4 i9 j10
      
      else if (c3 == TRUE & d4 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #c3 e5 f6 g7
      
      else if (c3 == TRUE & e5 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #c3 e5 f6 h8
      
      else if (c3 == TRUE & e5 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #c3 e5 f6 i9
      
      else if (c3 == TRUE & e5 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #c3 e5 f6 j10
      
      else if (c3 == TRUE & e5 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #c3 e5 g7 h8
      
      else if (c3 == TRUE & e5 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #c3 e5 g7 i9
      
      else if (c3 == TRUE & e5 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #c3 e5 g7 j10
      
      else if (c3 == TRUE & e5 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #c3 e5 h8 i9
      
      else if (c3 == TRUE & e5 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #c3 e5 h8 j10
      
      else if (c3 == TRUE & e5 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #c3 e5 i9 j10
      
      else if (c3 == TRUE & e5 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #c3 f6 g7 h8
      
      else if (c3 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #c3 f6 g7 i9
      
      else if (c3 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #c3 f6 g7 j10
      
      else if (c3 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #c3 f6 h8 i9
      
      else if (c3 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #c3 f6 h8 j10
      
      else if (c3 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #c3 f6 i9 j10
      
      else if (c3 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #c3 g7 h8 i9
      
      else if (c3 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #c3 g7 h8 j10
      
      else if (c3 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #c3 g7 i9 j10
      
      else if (c3 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #c3 h8 i9 j10
      
      else if (c3 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #d4 e5 f6 g7
      
      else if (d4 == TRUE & e5 == TRUE & f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      #d4 e5 f6 h8
      
      else if (d4 == TRUE & e5 == TRUE & f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      #d4 e5 f6 i9
      
      else if (d4 == TRUE & e5 == TRUE & f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      #d4 e5 f6 j10
      
      else if (d4 == TRUE & e5 == TRUE & f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      #d4 e5 g7 h8
      
      else if (d4 == TRUE & e5 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #d4 e5 g7 i9
      
      else if (d4 == TRUE & e5 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #d4 e5 g7 j10
      
      else if (d4 == TRUE & e5 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #d4 e5 h8 i9
      
      else if (d4 == TRUE & e5 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #d4 e5 h8 j10
      
      else if (d4 == TRUE & e5 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #d4 e5 i9 j10
      
      else if (d4 == TRUE & e5 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #d4 f6 g7 h8
      
      else if (d4 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #d4 f6 g7 i9
      
      else if (d4 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #d4 f6 g7 j10
      
      else if (d4 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #d4 f6 h8 i9
      
      else if (d4 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #d4 f6 h8 j10
      
      else if (d4 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #d4 f6 i9 j10
      
      else if (d4 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #d4 g7 h8 i9
      
      else if (d4 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #d4 g7 h8 j10
      
      else if (d4 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #d4 g7 i9 j10
      
      else if (d4 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #d4 h8 i9 j10
      
      else if (d4 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #e5 f6 g7 h8
      
      else if (e5 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      #e5 f6 g7 i9
      
      else if (e5 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      #e5 f6 g7 j10
      
      else if (e5 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #e5 f6 h8 i9
      
      else if (e5 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #e5 f6 h8 j10
      
      else if (e5 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #e5 f6 i9 j10
      
      else if (e5 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #e5 g7 h8 i9
      
      else if (e5 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #e5 g7 h8 j10
      
      else if (e5 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #e5 g7 i9 j10
      
      else if (e5 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #e5 h8 i9 j10
      
      else if (e5 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #f6 g7 h8 i9
      
      else if (f6 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      #f6 g7 h8 j10
      
      else if (f6 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      #f6 g7 i9 j10
      
      else if (f6 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #f6 h8 i9 j10
      
      else if (f6 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      #g7 h8 i9 j10
      
      else if (g7 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #a1 b2 c3
      
      
      else if (a1 == TRUE & b2 == TRUE & c3 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3)
      }
      
      
      
      #a1 b2 d4
      
      else if (a1 == TRUE & b2 == TRUE & d4 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4)
      }
      
      
      #a1 b2 e5
      
      else if (a1 == TRUE & b2 == TRUE & e5 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5)
      }
      
      
      #a1 b2 f6
      
      
      else if (a1 == TRUE & b2 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6)
      }
      
      
      #a1 b2 g7
      
      else if (a1 == TRUE & b2 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7)
      }
      
      
      #a1 b2 h8
      
      else if (a1 == TRUE & b2 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8)
      }
      
      
      #a1 b2 i9
      
      else if (a1 == TRUE & b2 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_9)
      }
      
      
      #a1 b2 j10
      
      else if (a1 == TRUE & b2 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_10)
      }
      
      #a1 c3 d4
      
      else if (a1 == TRUE & c3 == TRUE & d4 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4)
      }
      
      
      #a1 c3 e5
      
      else if (a1 == TRUE & c3 == TRUE & e5 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5)
      }
      
      
      #a1 c3 f6
      
      else if (a1 == TRUE & c3 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6)
      }
      
      
      #a1 c3 g7
      
      else if (a1 == TRUE & c3 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7)
      }
      
      #a1 c3 h8
      
      else if (a1 == TRUE & c3 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8)
      }
      
      
      #a1 c3 i9
      
      else if (a1 == TRUE & c3 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9)
      }
      
      
      #a1 c3 j10
      
      else if (a1 == TRUE & c3 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_10)
      }
      
      
      #a1 d4 e5
      
      else if (a1 == TRUE & d4 == TRUE & e5 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
      }
      
      
      #a1 d4 f6
      
      else if (a1 == TRUE & d4 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
      }
      
      
      #a1 d4 g7
      
      else if (a1 == TRUE & d4 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
      }
      
      
      #a1 d4 h8
      
      else if (a1 == TRUE & d4 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
      }
      
      
      #a1 d4 i9
      
      else if (a1 == TRUE & d4 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
      }
      
      
      #a1 d4 j10
      
      else if (a1 == TRUE & d4 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
      }
      
      
      #a1 e5 f6
      
      else if (a1 == TRUE & e5 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      
      #a1 e5 g7
      
      else if (a1 == TRUE & e5 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      
      #a1 e5 h8
      
      else if (a1 == TRUE & e5 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      
      #a1 e5 i9
      
      else if (a1 == TRUE & e5 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      
      #a1 e5 j10
      
      else if (a1 == TRUE & e5 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      #a1 f6 g7
      
      else if (a1 == TRUE & f6 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      
      #a1 f6 h8
      
      else if (a1 == TRUE & f6 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      
      #a1 f6 i9
      
      else if (a1 == TRUE & f6 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      
      #a1 f6 j10
      
      else if (a1 == TRUE & f6 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      
      #a1 g7 h8
      
      else if (a1 == TRUE & g7 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      
      #a1 g7 i9
      
      else if (a1 == TRUE & g7 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      
      #a1 g7 j10
      
      else if (a1 == TRUE & g7 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      #a1 h8 i9
      
      else if (a1 == TRUE & h8 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      
      #a1 h8 j10
      
      else if (a1 == TRUE & h8 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      
      #a1 i9 j10
      
      else if (a1 == TRUE & i9 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #b2 c3 d4
      
      else if (b2 == TRUE & c3 == TRUE & d4 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4)
      }
      
      
      #b2 c3 e5
      
      else if (b2 == TRUE & c3 == TRUE & e5 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5)
      }
      
      
      #b2 c3 f6
      
      else if (b2 == TRUE & c3 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6)
      }
      
      
      #b2 c3 g7
      
      else if (b2 == TRUE & c3 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7)
      }
      
      
      #b2 c3 h8
      
      else if (b2 == TRUE & c3 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8)
      }
      
      
      #b2 c3 i9
      
      else if (b2 == TRUE & c3 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9)
      }
      
      
      #b2 c3 j10
      
      else if (b2 == TRUE & c3 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_10)
      }
      
      
      #b2 d4 e5
      
      else if (b2 == TRUE & d4 == TRUE & e5 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
      }
      
      
      #b2 d4 f6
      
      else if (b2 == TRUE & d4 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
      }
      
      
      #b2 d4 g7
      
      else if (b2 == TRUE & d4 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
      }
      
      
      #b2 d4 h8
      
      else if (b2 == TRUE & d4 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
      }
      
      
      #b2 d4 i9
      
      else if (b2 == TRUE & d4 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
      }
      
      
      #b2 d4 j10
      
      else if (b2 == TRUE & d4 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
      }
      
      
      #b2 e5 f6
      
      else if (b2 == TRUE & e5 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      
      #b2 e5 g7
      
      else if (b2 == TRUE & e5 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      
      #b2 e5 h8
      
      else if (b2 == TRUE & e5 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      
      #b2 e5 i9
      
      else if (b2 == TRUE & e5 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      
      #b2 e5 j10
      
      else if (b2 == TRUE & e5 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      
      #b2 f6 g7
      
      else if (b2 == TRUE & f6 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      
      #b2 f6 h8
      
      else if (b2 == TRUE & f6 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      
      #b2 f6 i9
      
      else if (b2 == TRUE & f6 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      
      #b2 f6 j10
      
      else if (b2 == TRUE & f6 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      
      #b2 g7 h8
      
      else if (b2 == TRUE & g7 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      
      #b2 g7 i9
      
      else if (b2 == TRUE & g7 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      
      #b2 g7 j10
      
      else if (b2 == TRUE & g7 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      
      #b2 h8 i9
      
      else if (b2 == TRUE & h8 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      
      #b2 h8 j10
      
      else if (b2 == TRUE & h8 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      
      #b2 i9 j10
      
      else if (b2 == TRUE & i9 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #c3 d4 e5
      
      else if (c3 == TRUE & d4 == TRUE & e5 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
      }
      
      
      #c3 d4 f6
      
      else if (c3 == TRUE & d4 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
      }
      
      
      #c3 d4 g7
      
      else if (c3 == TRUE & d4 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
      }
      
      
      #c3 d4 h8
      
      else if (c3 == TRUE & d4 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
      }
      
      
      #c3 d4 i9
      
      else if (c3 == TRUE & d4 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
      }
      
      
      #c3 d4 j10
      
      else if (c3 == TRUE & d4 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
      }
      
      
      #c3 e5 f6
      
      else if (c3 == TRUE & e5 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      
      #c3 e5 g7
      
      else if (c3 == TRUE & e5 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      
      #c3 e5 h8
      
      else if (c3 == TRUE & e5 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      
      #c3 e5 i9
      
      else if (c3 == TRUE & e5 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      
      #c3 e5 j10
      
      else if (c3 == TRUE & e5 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      
      #c3 f6 g7
      
      else if (c3 == TRUE & f6 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      
      #c3 f6 h8
      
      else if (c3 == TRUE & f6 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      
      #c3 f6 i9
      
      else if (c3 == TRUE & f6 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      
      #c3 f6 j10
      
      else if (c3 == TRUE & f6 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      
      #c3 g7 h8
      
      else if (c3 == TRUE & g7 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      
      #c3 g7 i9
      
      else if (c3 == TRUE & g7 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      
      #c3 g7 j10
      
      else if (c3 == TRUE & g7 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      
      #c3 h8 i9
      
      else if (c3 == TRUE & h8 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      
      #c3 h8 j10
      
      else if (c3 == TRUE & h8 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      
      #c3 i9 j10
      
      else if (c3 == TRUE & i9 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #d4 e5 f6
      
      else if (d4 == TRUE & e5 == TRUE & f6 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      
      #d4 e5 g7
      
      else if (d4 == TRUE & e5 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      
      #d4 e5 h8
      
      else if (d4 == TRUE & e5 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      
      #d4 e5 i9
      
      else if (d4 == TRUE & e5 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      
      #d4 e5 j10
      
      else if (d4 == TRUE & e5 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      
      #d4 f6 g7
      
      else if (d4 == TRUE & f6 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      
      #d4 f6 h8
      
      else if (d4 == TRUE & f6 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      
      #d4 f6 i9
      
      else if (d4 == TRUE & f6 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      
      #d4 f6 j10
      
      else if (d4 == TRUE & f6 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      
      #d4 g7 h8
      
      else if (d4 == TRUE & g7 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      
      #d4 g7 i9
      
      else if (d4 == TRUE & g7 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      
      #d4 g7 j10
      
      else if (d4 == TRUE & g7 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      
      #d4 h8 i9
      
      else if (d4 == TRUE & h8 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      
      #d4 h8 j10
      
      else if (d4 == TRUE & h8 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      
      #d4 i9 j10
      
      else if (d4 == TRUE & i9 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #e5 f6 g7
      
      else if (e5 == TRUE & f6 == TRUE & g7 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      
      #e5 f6 h8
      
      else if (e5 == TRUE & f6 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      
      #e5 f6 i9
      
      else if (e5 == TRUE & f6 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      
      #e5 f6 j10
      
      else if (e5 == TRUE & f6 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      
      #e5 g7 h8
      
      else if (e5 == TRUE & g7 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      
      #e5 g7 i9
      
      else if (e5 == TRUE & g7 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      
      #e5 g7 j10
      
      else if (e5 == TRUE & g7 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      
      #e5 h8 i9
      
      else if (e5 == TRUE & h8 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      
      #e5 h8 j10
      
      else if (e5 == TRUE & h8 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      
      #e5 i9 j10
      
      else if (e5 == TRUE & i9 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #f6 g7 h8
      
      else if (f6 == TRUE & g7 == TRUE & h8 == TRUE )
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      
      #f6 g7 i9
      
      else if (f6 == TRUE & g7 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      
      #f6 g7 j10
      
      else if (f6 == TRUE & g7 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      
      #f6 h8 i9
      
      else if (f6 == TRUE & h8 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      
      #f6 h8 j10
      
      else if (f6 == TRUE & h8 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      
      #f6 i9 j10
      
      else if (f6 == TRUE & i9 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #g7 h8 i9
      
      else if (g7 == TRUE & h8 == TRUE & i9 == TRUE )
      {
        outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      
      #g7 h8 j10
      
      else if (g7 == TRUE & h8 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      
      #g7 i9 j10
      
      else if (g7 == TRUE & i9 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #h8 i9 j10
      
      else if (h8 == TRUE & i9 == TRUE & j10 == TRUE )
      {
        outliers_SMALLCLUSTER_8  <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #a1 b2
      
      
      else if (a1 == TRUE & b2 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2)
      }
      
      
      #a1 c3
      
      
      else if (a1 == TRUE & c3 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3)
      }
      
      
      #a1 d4
      
      else if (a1 == TRUE & d4 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4)
      }
      
      
      #a1 e5
      
      else if (a1 == TRUE & e5 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5)
      }
      
      
      #a1 f6
      
      else if (a1 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6)
      }
      
      
      #a1 g7
      
      else if (a1 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7)
      }
      
      
      #a1 h8
      
      else if (a1 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_8)
      }
      
      
      #a1 i9
      
      else if (a1 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_9)
      }
      
      
      #a1 j10
      
      else if (a1 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_10)
      }
      
      
      #b2 c3
      
      else if (b2 == TRUE & c3 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3)
      }
      
      
      #b2 d4
      
      else if (b2 == TRUE & d4 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4)
      }
      
      
      #b2 e5
      
      else if (b2 == TRUE & e5 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5)
      }
      
      
      #b2 f6
      
      else if (b2 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6)
      }
      
      
      #b2 g7
      
      else if (b2 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7)
      }
      
      
      #b2 h8
      
      else if (b2 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8)
      }
      
      
      #b2 i9
      
      else if (b2 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_9)
      }
      
      
      #b2 j10
      
      else if (b2 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_10)
      }
      
      
      #c3 d4
      
      else if (c3 == TRUE & d4 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4)
      }
      
      
      #c3 e5
      
      else if (c3 == TRUE & e5 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5)
      }
      
      
      #c3 f6
      
      else if (c3 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6)
      }
      
      
      #c3 g7
      
      else if (c3 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7)
      }
      
      
      #c3 h8
      
      else if (c3 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8)
      }
      
      
      #c3 i9
      
      else if (c3 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9)
      }
      
      
      #c3 j10
      
      else if (c3 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_10)
      }
      
      
      #d4 e5
      
      else if (d4 == TRUE & e5 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
      }
      
      
      #d4 f6
      
      else if (d4 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
      }
      
      
      #d4 g7
      
      else if (d4 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
      }
      
      
      #d4 h8
      
      else if (d4 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
      }
      
      
      #d4 i9
      
      else if (d4 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
      }
      
      
      #d4 j10
      
      else if (d4 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
      }
      
      
      #e5 f6
      
      else if (e5 == TRUE & f6 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
      }
      
      
      #e5 g7
      
      else if (e5 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
      }
      
      
      #e5 h8
      
      else if (e5 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
      }
      
      
      #e5 i9
      
      else if (e5 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
      }
      
      
      #e5 j10
      
      else if (e5 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
      }
      
      
      #f6 g7
      
      else if (f6 == TRUE & g7 == TRUE)
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
      }
      
      
      #f6 h8
      
      else if (f6 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
      }
      
      
      #f6 i9
      
      else if (f6 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
      }
      
      
      #f6 j10
      
      else if (f6 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
      }
      
      
      #g7 h8
      
      else if (g7 == TRUE & h8 == TRUE)
      {
        outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
      }
      
      
      #g7 i9
      
      else if (g7 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
      }
      
      
      #g7 j10
      
      else if (g7 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
      }
      
      
      #h8 i9
      
      else if (h8 == TRUE & i9 == TRUE)
      {
        outliers_SMALLCLUSTER_8  <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
      }
      
      
      #h8 j10
      
      else if (h8 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_8  <-which(cutree_wektor == 8)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
      }
      
      
      #i9 j10
      
      else if (i9 == TRUE & j10 == TRUE)
      {
        outliers_SMALLCLUSTER_9  <-which(cutree_wektor == 9)
        outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
      }
      
      
      #a1
      
      else if (a1 == TRUE)
      {
        outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1)
      }
      
      #b2
      
      else if (b2 == TRUE)
      {
        outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2)
      }
      
      #c3
      
      else if (c3 == TRUE)
      {
        outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3)
      }
      
      #d4
      
      else if (d4 == TRUE)
      {
        outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4)
      }
      
      #e5
      
      else if (e5 == TRUE)
      {
        outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5)
      }
      
      #f6
      
      else if (f6 == TRUE)
      {
        outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6)
      }
      
      #g7
      
      else if (g7 == TRUE)
      {
        outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7)
      }
      
      #h8
      
      else if (h8 == TRUE)
      {
        outliers_SMALLCLUSTER_8  <-which(cutree_wektor == 8)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_8)
      }
      
      #i9
      
      else if (i9 == TRUE)
      {
        outliers_SMALLCLUSTER_9  <-which(cutree_wektor == 9)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_9)
      }
      
      
      #j10
      
      else if (j10 == TRUE)
      {
        outliers_SMALLCLUSTER_10  <-which(cutree_wektor == 10)
        
        outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_10)
      }
      
        else
        {
          outliers_SMALLCLUSTER <- "t"
        }
        
       
###############################################################################
  
  if ( outliers_SMALLCLUSTER == "t")
  {    wyniki[12,3] <- "not small clusters"}
else
       {wyniki[12,3] <- toString(outliers_SMALLCLUSTER)}
       
        if (outliers_SMALLCLUSTER == "t")
        { data_inf4 <- j }
     else
     {
    x <- outliers_SMALLCLUSTER
      data_inf4 <- j[-x, ]
   }   
      
      
      gower_mat1 <- gower.dist ( data_inf4 ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL )
      gower_mat1[is.nan(gower_mat1)] <- 0
      gower_mat <-as.dist(gower_mat1) 
     
      clusters <- hclust(gower_mat,method = input$metoda) 
      wyniki[11,3] <- input$metoda
      d1 <- gower_mat
      hc <- clusters
      d2 <- cophenetic(hc) 
      wyniki[10,3] <- cor(d1, d2)
      
      cutree_wektor <- cutree(clusters, k = cutree_k_SMALLCLUSTER)  
      si3 <- silhouette(cutree_wektor, gower_mat) 
      icq1 <- mean(si3[,"sil_width"])
      wyniki[1,3]<-icq1
      
      wyniki[2,3]<-dunn(gower_mat, cutree_wektor)
      index_DB <- index.DB(gower_mat, cutree_wektor, centrotypes="centroids", p = 2, q = 2)
      wyniki[3,3]<-index_DB[1]
      icq_index_G1 <- index.G1(gower_mat1, cutree_wektor, d = NULL, centrotypes = "centroids")
      wyniki[4,3]<-icq_index_G1
      icq_index_G3 <- index.G3(gower_mat, cutree_wektor)
      wyniki[5,3]<-icq_index_G3
      icq_index_G2 <- index.G2(gower_mat, cutree_wektor)
      wyniki[6,3]<-icq_index_G2
      
    }
    
  
       
       if ( ((as.numeric(wyniki[1,2])) < (as.numeric((wyniki[1,3])))) &
            ((as.numeric(wyniki[2,2])) <= (as.numeric((wyniki[2,3])))) & 
            ((as.numeric(wyniki[3,2])) > (as.numeric((wyniki[3,3])))) & 
            ((as.numeric(wyniki[4,2])) < (as.numeric((wyniki[4,3])))) & 
            ((as.numeric(wyniki[5,2])) >= (as.numeric((wyniki[5,3])))) & 
            ((as.numeric(wyniki[6,2])) <= (as.numeric((wyniki[6,3])))) & 
           ((as.numeric(wyniki[10,2])) < (as.numeric((wyniki[10,3])))) )
  {
      tableHTML(wyniki, caption = 'Analiza wpływu reguł nietypowych na jakość grupowania:', footer = '*Spróbuj dobrać pola wyboru i paski wyszukiwania tak, aby niebieski nagłówek tabeli zmienił się na zielony',  
                           collapse = 'separate_shiny', spacing = '6px 5px', escape = TRUE, rownames = FALSE) %>%
        add_css_row(css = list('background-color', '#e6f0ff'),
                    rows = odd(2:13)) %>%
         add_css_row(css = list('background-color', '#ffffff'),
                                  rows = even(2:13))  %>%
        
        add_css_row(css = list('background-color', '#90EE90'),  
                    rows = odd(1:1))
      }
    
    else 
    {
      tableHTML(wyniki, caption = 'Analiza wpływu reguł nietypowych na jakość grupowania:', footer = '*Spróbuj dobrać pola wyboru i paski wyszukiwania tak, aby niebieski nagłówek tabeli zmienił się na zielony',  
                collapse = 'separate_shiny', spacing = '6px 5px', escape = TRUE, rownames = FALSE) %>%
        add_css_row(css = list('background-color', '#e6f0ff'),
                    rows = odd(2:13)) %>%
        add_css_row(css = list('background-color', '#ffffff'),
                    rows = even(2:13))  %>%
        add_css_row(css = list('background-color', '#C6D9F6'),
                    rows = odd(1:1))
    }  
      
                             
  })
  
  
  output$wyjscieDendrogram <- renderPlot({
  
    observeEvent(input$show, {
      showNotification("Message text",
                       action = session$reload(), "Reload page")
      
    })
    
    fn <- "wczytanyZbior.csv"
    #Check its existence
    if (!file.exists(fn)) 
     return(NULL)
    qoute2 <- readRDS(file = "quote2.rds")
    my_data <- fread ("wczytanyZbior.csv", sep=',', header = TRUE, stringsAsFactors = FALSE, quote = qoute2)
    my_data[my_data=="NA"] <- NA
    
    gower_mat1 <- gower.dist ( my_data ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL )
    gower_mat1[is.nan(gower_mat1)] <- 0
    gower_mat <-as.dist(gower_mat1) 
    clusters <- hclust(gower_mat, method = input$metoda)
    
    J<-readRDS(file = "file1.rds")
    
    collection_30 = unlist(gregexpr(pattern =",",toString(J)))
    collection_40<-str_sub(string=toString(J), start = 1, end = collection_30[1]-1)
    
    
  plot(fviz_dend(clusters, k = input$liczbaGrup[1],font.main = 20, main = paste("Dendrogram przed usunięciem reguł nietypowych, dataset = ",collection_40), 
                   xlab = paste("Metoda łączenia grup:",input$metoda,"|", "Miara odległości między regułami: Gower " ),cex = 0.5,  # label size
                 k_colors = c("#2E9FDF", "#FC4E07", "#07fc4e", "#3339FF", "#00AFBB", "#FF33F0", "#E7B800", "#C0C0C0","#800080", "#008080", "#808000", "#2F4F4F", "#8B0000", "#FF1493","#BDB76B", "#FFD700", "#483D8B", "#00FF7F","#00CED1","#0000CD" ), 
           color_labels_by_k = TRUE, # color labels by groups
          rect = FALSE # Add rectangle around groups 
    ))
    
    
    output$plot2 <- renderPlot({
      
      
      if (input$algorytm == "LOF")
        
      {
        
        qoute2 <- readRDS(file = "quote2.rds")
        
        g <- fread ("wczytanyZbior.csv", sep=',', header = TRUE, stringsAsFactors = FALSE, quote = qoute2)
        g[g=="NA"] <- NA
        
        RULES1<-nrow(g)
        howmany <-(input$procentRegulNietypowych[1]/100)*RULES1
        howmany<-ceiling(howmany)
     
        set_k <- input$kSasiadow[1] 
     
        outlier.scores <- lofactor(gower_mat1, k = set_k)
        outliers_LOF <- order(outlier.scores, decreasing = T) [1:howmany]
        
        p <- outliers_LOF
        data_inf1 <- g[-p, ]
           
        gower_mat1 <- gower.dist (data_inf1 ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL)
        gower_mat1[is.nan(gower_mat1)] <- 0
        gower_mat <-as.dist(gower_mat1) 
        clusters <- hclust(gower_mat, method = input$metoda)   
        
        plot(fviz_dend(clusters, k = input$liczbaGrup[1],font.main = 20, main = paste("Dendrogram po usunięciu reguł nietypowych - metoda LOF, dataset = ",collection_40), 
                       xlab = paste("Metoda łączenia grup:",input$metoda,"|", "Miara odległości między regułami: Gower " ),cex = 0.5,   # label size
                       k_colors = c("#2E9FDF", "#FC4E07", "#07fc4e", "#3339FF", "#00AFBB", "#FF33F0", "#E7B800", "#C0C0C0","#800080", "#008080", "#808000", "#2F4F4F", "#8B0000", "#FF1493","#BDB76B", "#FFD700", "#483D8B", "#00FF7F","#00CED1","#0000CD" ),
                       color_labels_by_k = TRUE, # color labels by groups
                       rect = FALSE # Add rectangle around groups 
        ))
        
      }
      
      
      
      if (input$algorytm == "COF")
        
      {
        
        qoute2 <- readRDS(file = "quote2.rds")
        f <- fread ("wczytanyZbior.csv", sep=',', header = TRUE,stringsAsFactors = FALSE, quote = qoute2)
        f[f=="NA"] <- NA
        RULES1<-nrow(f)
        howmany <-(input$procentRegulNietypowych[1]/100)*RULES1
        howmany<-ceiling(howmany)
        set_k <- input$kSasiadow[1] 
        outlier.scores <- COF(gower_mat1, k = set_k)
        outliers_COF <- order(outlier.scores, decreasing = T) [1:howmany]
        q <- outliers_COF
        data_inf2 <- f[-q, ]
        
        gower_mat1 <- gower.dist ( data_inf2 ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL )
        gower_mat1[is.nan(gower_mat1)] <- 0
        gower_mat <-as.dist(gower_mat1) 
        
        clusters <- hclust(gower_mat,method = input$metoda) 
        
        
        plot(fviz_dend(clusters, k = input$liczbaGrup[1],font.main = 20, main = paste("Dendrogram po usunięciu reguł nietypowych - metoda COF, dataset = ",collection_40), 
                       xlab = paste("Metoda łączenia grup:",input$metoda,"|", "Miara odległości między regułami: Gower " ),cex = 0.5,   # label size
                       k_colors = c("#2E9FDF", "#FC4E07", "#07fc4e", "#3339FF", "#00AFBB", "#FF33F0", "#E7B800", "#C0C0C0","#800080", "#008080", "#808000", "#2F4F4F", "#8B0000", "#FF1493","#BDB76B", "#FFD700", "#483D8B", "#00FF7F","#00CED1","#0000CD" ),
                       color_labels_by_k = TRUE, # color labels by groups
                       rect = FALSE # Add rectangle around groups 
        ))
        
      }
      
      if (input$algorytm == "KMEANS")
        
      {
        
        qoute2 <- readRDS(file = "quote2.rds")
        h <- fread ("wczytanyZbior.csv", sep=',', header = TRUE,stringsAsFactors = FALSE, quote = qoute2)
        h[h=="NA"] <- NA
        RULES1<-nrow(h)
        howmany <-(input$procentRegulNietypowych[1]/100)*RULES1
        howmany<-ceiling(howmany)
      
        set_k <- input$setSeed[1] 
        RNGkind(sample.kind = "Rounding")
        set.seed(set_k) 
        DECISION_VALUES1<-readRDS(file = "RULES1.rds")
        kmeans.result <- kmeans(gower_mat1 , centers = DECISION_VALUES1)
        centers <- kmeans.result$centers[kmeans.result$cluster, ]
        distances <- sqrt(rowSums((as.numeric(unlist(gower_mat1)) - centers)^2))
        outliers_kmeans <- order( distances, decreasing = T ) [1:howmany]
        
        r <- outliers_kmeans
        data_inf3 <- h[-r, ]
        
        gower_mat1 <- gower.dist ( data_inf3 ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL )
        gower_mat1[is.nan(gower_mat1)] <- 0
        gower_mat <-as.dist(gower_mat1) 
        
        
        clusters <- hclust(gower_mat, method =  input$metoda) 
        
        
        plot(fviz_dend(clusters, k = input$liczbaGrup[1],font.main = 20, main = paste("Dendrogram po usunięciu reguł nietypowych - metoda KMEANS, dataset = ",collection_40), 
                       xlab = paste("Metoda łączenia grup:",input$metoda,"|", "Miara odległości między regułami: Gower " ),cex = 0.5,   # label size
                       k_colors = c("#2E9FDF", "#FC4E07", "#07fc4e", "#3339FF", "#00AFBB", "#FF33F0", "#E7B800", "#C0C0C0","#800080", "#008080", "#808000", "#2F4F4F", "#8B0000", "#FF1493","#BDB76B", "#FFD700", "#483D8B", "#00FF7F","#00CED1","#0000CD" ),
                       color_labels_by_k = TRUE, # color labels by groups
                       rect = FALSE # Add rectangle around groups 
        ))
        
        
      }
      
      
      if (input$algorytm == "SMALLCLUSTER")
        
      {
        
        qoute2 <- readRDS(file = "quote2.rds")
        
        t <- fread ("wczytanyZbior.csv", sep=',', header = TRUE, stringsAsFactors = FALSE, quote = qoute2)
        t[t=="NA"] <- NA
        cutree_k_SMALLCLUSTER <- input$liczbaGrup[1]
        cutree_k_SMALLCLUSTER_outliers <- input$smallCluster[1]
        
        if (cutree_k_SMALLCLUSTER_outliers > nrow(t))
          { return("Wybrana liczba klastrów jest większa niż liczba obserwacji") }
        
        gower_mat1 <- gower.dist ( t ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL )
        gower_mat1[is.nan(gower_mat1)] <- 0
        gower_mat <-as.dist(gower_mat1) 
        clusters <- hclust(gower_mat,method = input$metoda)
        cutree_wektor <- cutree(clusters, k = cutree_k_SMALLCLUSTER_outliers)
        
        
        ilosc_obserwacji <-length(cutree_wektor) 
        ilosc_obserwacji_1 <-length(which(cutree_wektor == 1))
        ilosc_obserwacji_2 <-length(which(cutree_wektor == 2))
        ilosc_obserwacji_3 <-length(which(cutree_wektor == 3))
        ilosc_obserwacji_4 <-length(which(cutree_wektor == 4))
        ilosc_obserwacji_5 <-length(which(cutree_wektor == 5))
        ilosc_obserwacji_6 <-length(which(cutree_wektor == 6))
        ilosc_obserwacji_7 <-length(which(cutree_wektor == 7))
        ilosc_obserwacji_8 <-length(which(cutree_wektor == 8))
        ilosc_obserwacji_9 <-length(which(cutree_wektor == 9))
        ilosc_obserwacji_10 <-length(which(cutree_wektor == 10))
        
        
        
        
        if ((ilosc_obserwacji_1/ilosc_obserwacji) > 0 & (ilosc_obserwacji_1/ilosc_obserwacji) <= input$procentZbioru)
        {
          a1<-TRUE
        }
        else
        {
          a1<-FALSE
        }
        
        
        if ((ilosc_obserwacji_2/ilosc_obserwacji) > 0 & (ilosc_obserwacji_2/ilosc_obserwacji) <= input$procentZbioru)
        {
          b2<-TRUE
        }
        else
        {
          b2<-FALSE
        }
        
        if ((ilosc_obserwacji_3/ilosc_obserwacji) > 0 & (ilosc_obserwacji_3/ilosc_obserwacji) <= input$procentZbioru)
        {
          c3<-TRUE
        }
        else
        {
          c3<-FALSE
        }
        
        if ((ilosc_obserwacji_4/ilosc_obserwacji) > 0 & (ilosc_obserwacji_4/ilosc_obserwacji) <= input$procentZbioru)
        {
          d4<-TRUE
        }
        else
        {
          d4<-FALSE
        }
        
        if ((ilosc_obserwacji_5/ilosc_obserwacji) > 0 & (ilosc_obserwacji_5/ilosc_obserwacji) <= input$procentZbioru)
        {
          e5<-TRUE
        }
        else
        {
          e5<-FALSE
        }
        
        if ((ilosc_obserwacji_6/ilosc_obserwacji) > 0 & (ilosc_obserwacji_6/ilosc_obserwacji) <= input$procentZbioru)
        {
          f6<-TRUE
        }
        else
        {
          f6<-FALSE
        }
        
        if ((ilosc_obserwacji_7/ilosc_obserwacji) > 0 & (ilosc_obserwacji_7/ilosc_obserwacji) <= input$procentZbioru)
        {
          g7<-TRUE
        }
        else
        {
          g7<-FALSE
        }
        
        if ((ilosc_obserwacji_8/ilosc_obserwacji) > 0 & (ilosc_obserwacji_8/ilosc_obserwacji) <= input$procentZbioru)
        {
          h8<-TRUE
        }
        else
        {
          h8<-FALSE
        }
        
        if ((ilosc_obserwacji_9/ilosc_obserwacji) > 0 & (ilosc_obserwacji_9/ilosc_obserwacji) <= input$procentZbioru)
        {
          i9<-TRUE
        }
        else
        {
          i9<-FALSE
        }
        
        if ((ilosc_obserwacji_10/ilosc_obserwacji) > 0 & (ilosc_obserwacji_10/ilosc_obserwacji) <= input$procentZbioru)
        {
          j10<-TRUE
        }
        else
        {
          j10<-FALSE
        }
        
        
        ################################################################################################################
        
        #a1 b2 c3 d4
        
        if (a1 == TRUE & b2 == TRUE & c3 == TRUE & d4 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4)
        }
        
        #a1 b2 c3 e5
        
        else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & e5 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5)
        }
        
        
        #a1 b2 c3 f6
        
        else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6)
        }
        
        
        #a1 b2 c3 g7
        
        else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7)
        }
        
        
        #a1 b2 c3 h8
        
        else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8)
        }
        
        
        #a1 b2 c3 i9
        
        else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9)
        }
        
        
        #a1 b2 c3 j10
        
        else if (a1 == TRUE & b2 == TRUE & c3 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_10)
        }
        
        #a1 b2 d4 e5
        
        else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & e5 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
        }
        
        #a1 b2 d4 f6
        
        else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
        }
        
        #a1 b2 d4 g7
        
        else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
        }
        
        #a1 b2 d4 h8
        
        else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
        }
        
        #a1 b2 d4 i9
        
        else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
        }
        
        #a1 b2 d4 j10
        
        else if (a1 == TRUE & b2 == TRUE & d4 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
        }
        
        #a1 b2 e5 f6
        
        else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        #a1 b2 e5 g7
        
        else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        #a1 b2 e5 h8
        
        else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        #a1 b2 e5 i9
        
        else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        #a1 b2 e5 j10
        
        else if (a1 == TRUE & b2 == TRUE & e5 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        #a1 b2 f6 g7
        
        else if (a1 == TRUE & b2 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #a1 b2 f6 h8
        
        else if (a1 == TRUE & b2 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #a1 b2 f6 i9
        
        else if (a1 == TRUE & b2 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #a1 b2 f6 j10
        
        else if (a1 == TRUE & b2 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #a1 b2 g7 h8
        
        else if (a1 == TRUE & b2 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #a1 b2 g7 i9
        
        else if (a1 == TRUE & b2 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #a1 b2 g7 j10
        
        else if (a1 == TRUE & b2 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #a1 b2 h8 i9
        
        else if (a1 == TRUE & b2 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #a1 b2 h8 j10
        
        else if (a1 == TRUE & b2 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #a1 b2 i9 j10
        
        else if (a1 == TRUE & b2 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #a1 c3 d4 e5
        
        else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & e5 == TRUE)
        {
          outliers_SMALLCLUSTER_1 <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
        }
        
        #a1 c3 d4 f6
        
        else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
        }
        
        #a1 c3 d4 g7
        
        else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
        }
        
        #a1 c3 d4 h8
        
        else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
        }
        
        #a1 c3 d4 i9
        
        else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
        }
        
        #a1 c3 d4 j10
        
        else if (a1 == TRUE & c3 == TRUE & d4 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
        }
        
        #a1 c3 e5 f6
        
        else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        #a1 c3 e5 g7
        
        else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        #a1 c3 e5 h8
        
        else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        #a1 c3 e5 i9
        
        else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        #a1 c3 e5 j10
        
        else if (a1 == TRUE & c3 == TRUE & e5 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        #a1 c3 f6 g7
        
        else if (a1 == TRUE & c3 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #a1 c3 f6 h8
        
        else if (a1 == TRUE & c3 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #a1 c3 f6 i9
        
        else if (a1 == TRUE & c3 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #a1 c3 f6 j10
        
        else if (a1 == TRUE & c3 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #a1 c3 g7 h8
        
        else if (a1 == TRUE & c3 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #a1 c3 g7 i9
        
        else if (a1 == TRUE & c3 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #a1 c3 g7 j10
        
        else if (a1 == TRUE & c3 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #a1 c3 h8 i9
        
        else if (a1 == TRUE & c3== TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #a1 c3 h8 j10
        
        else if (a1 == TRUE & c3 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #a1 c3 i9 j10
        
        else if (a1 == TRUE & c3 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #a1 d4 e5 f6
        
        else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        #a1 d4 e5 g7
        
        else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        #a1 d4 e5 h8
        
        else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        #a1 d4 e5 i9
        
        else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        #a1 d4 e5 j10
        
        else if (a1 == TRUE & d4 == TRUE & e5 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        #a1 d4 f6 g7
        
        else if (a1 == TRUE & d4 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #a1 d4 f6 h8
        
        else if (a1 == TRUE & d4 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #a1 d4 f6 i9
        
        else if (a1 == TRUE & d4 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #a1 d4 f6 j10
        
        else if (a1 == TRUE & d4 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #a1 d4 g7 h8
        
        else if (a1 == TRUE & d4 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #a1 d4 g7 i9
        
        else if (a1 == TRUE & d4 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #a1 d4 g7 j10
        
        else if (a1 == TRUE & d4 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #a1 d4 h8 i9
        
        else if (a1 == TRUE & d4 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #a1 d4 h8 j10
        
        else if (a1 == TRUE & d4 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #a1 d4 i9 j10
        
        else if (a1 == TRUE & d4 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #a1 e5 f6 g7
        
        else if (a1 == TRUE & e5 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #a1 e5 f6 h8
        
        else if (a1 == TRUE & e5 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #a1 e5 f6 i9
        
        else if (a1 == TRUE & e5 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #a1 e5 f6 j10
        
        else if (a1 == TRUE & e5 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #a1 e5 g7 h8
        
        else if (a1 == TRUE & e5 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #a1 e5 g7 i9
        
        else if (a1 == TRUE & e5 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #a1 e5 g7 j10
        
        else if (a1 == TRUE & e5 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #a1 e5 h8 i9
        
        else if (a1 == TRUE & e5 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #a1 e5 h8 j10
        
        else if (a1 == TRUE & e5 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #a1 e5 i9 j10
        
        else if (a1 == TRUE & e5 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #a1 f6 g7 h8
        
        else if (a1 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #a1 f6 g7 i9
        
        else if (a1 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #a1 f6 g7 j10
        
        else if (a1 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #a1 f6 h8 i9
        
        else if (a1 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #a1 f6 h8 j10
        
        else if (a1 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #a1 f6 i9 j10
        
        else if (a1 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #a1 g7 h8 i9
        
        else if (a1 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #a1 g7 h8 j10
        
        else if (a1 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #a1 g7 i9 j10
        
        else if (a1 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #a1 h8 i9 j10
        
        else if (a1 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #b2 c3 d4 e5
        
        else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & e5 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
        }
        
        #b2 c3 d4 f6
        
        else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
        }
        
        #b2 c3 d4 g7
        
        else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
        }
        
        #b2 c3 d4 h8
        
        else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
        }
        
        #b2 c3 d4 i9
        
        else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
        }
        
        #b2 c3 d4 j10
        
        else if (b2 == TRUE & c3 == TRUE & d4 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
        }
        
        #b2 c3 e5 f6
        
        else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        #b2 c3 e5 g7
        
        else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        #b2 c3 e5 h8
        
        else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        #b2 c3 e5 i9
        
        else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        #b2 c3 e5 j10
        
        else if (b2 == TRUE & c3 == TRUE & e5 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        #b2 c3 f6 g7
        
        else if (b2 == TRUE & c3 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #b2 c3 f6 h8
        
        else if (b2 == TRUE & c3 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #b2 c3 f6 i9
        
        else if (b2 == TRUE & c3 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #b2 c3 f6 j10
        
        else if (b2 == TRUE & c3 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #b2 c3 g7 h8
        
        else if (b2 == TRUE & c3 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #b2 c3 g7 i9
        
        else if (b2 == TRUE & c3 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #b2 c3 g7 j10
        
        else if (b2 == TRUE & c3 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #b2 c3 h8 i9
        
        else if (b2 == TRUE & c3 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #b2 c3 h8 j10
        
        else if (b2 == TRUE & c3 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #b2 c3 i9 j10
        
        else if (b2 == TRUE & c3 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #b2 d4 e5 f6
        
        else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        #b2 d4 e5 g7
        
        else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        #b2 d4 e5 h8
        
        else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        #b2 d4 e5 i9
        
        else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        #b2 d4 e5 j10
        
        else if (b2 == TRUE & d4 == TRUE & e5 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        #b2 d4 f6 g7
        
        else if (b2 == TRUE & d4 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #b2 d4 f6 h8
        
        else if (b2 == TRUE & d4 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #b2 d4 f6 i9
        
        else if (b2 == TRUE & d4 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #b2 d4 f6 j10
        
        else if (b2 == TRUE & d4 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #b2 d4 g7 h8
        
        else if (b2 == TRUE & d4 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #b2 d4 g7 i9
        
        else if (b2 == TRUE & d4 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #b2 d4 g7 j10
        
        else if (b2 == TRUE & d4 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #b2 d4 h8 i9
        
        else if (b2 == TRUE & d4 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #b2 d4 h8 j10
        
        else if (b2 == TRUE & d4 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #b2 d4 i9 j10
        
        else if (b2 == TRUE & d4 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #b2 e5 f6 g7
        
        else if (b2 == TRUE & e5 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #b2 e5 f6 h8
        
        else if (b2 == TRUE & e5 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #b2 e5 f6 i9
        
        else if (b2 == TRUE & e5 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #b2 e5 f6 j10
        
        else if (b2 == TRUE & e5 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #b2 e5 g7 h8
        
        else if (b2 == TRUE & e5 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #b2 e5 g7 i9
        
        else if (b2 == TRUE & e5 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #b2 e5 g7 j10
        
        else if (b2 == TRUE & e5 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #b2 e5 h8 i9
        
        else if (b2 == TRUE & e5 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #b2 e5 h8 j10
        
        else if (b2 == TRUE & e5 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #b2 e5 i9 j10
        
        else if (b2 == TRUE & e5 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #b2 f6 g7 h8
        
        else if (b2 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #b2 f6 g7 i9
        
        else if (b2 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #b2 f6 g7 j10
        
        else if (b2 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #b2 f6 h8 i9
        
        else if (b2 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #b2 f6 h8 j10
        
        else if (b2 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #b2 f6 i9 j10
        
        else if (b2 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #b2 g7 h8 i9
        
        else if (b2 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #b2 g7 h8 j10
        
        else if (b2 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #b2 g7 i9 j10
        
        else if (b2 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #b2 h8 i9 j10
        
        else if (b2 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #c3 d4 e5 f6
        
        else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        #c3 d4 e5 g7
        
        else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        #c3 d4 e5 h8
        
        else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        #c3 d4 e5 i9
        
        else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        #c3 d4 e5 j10
        
        else if (c3 == TRUE & d4 == TRUE & e5 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        #c3 d4 f6 g7
        
        else if (c3 == TRUE & d4 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #c3 d4 f6 h8
        
        else if (c3 == TRUE & d4 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #c3 d4 f6 i9
        
        else if (c3 == TRUE & d4 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #c3 d4 f6 j10
        
        else if (c3 == TRUE & d4 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #c3 d4 g7 h8
        
        else if (c3 == TRUE & d4 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #c3 d4 g7 i9
        
        else if (c3 == TRUE & d4 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #c3 d4 g7 j10
        
        else if (c3 == TRUE & d4 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #c3 d4 h8 i9
        
        else if (c3 == TRUE & d4 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #c3 d4 h8 j10
        
        else if (c3 == TRUE & d4 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #c3 d4 i9 j10
        
        else if (c3 == TRUE & d4 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #c3 e5 f6 g7
        
        else if (c3 == TRUE & e5 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #c3 e5 f6 h8
        
        else if (c3 == TRUE & e5 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #c3 e5 f6 i9
        
        else if (c3 == TRUE & e5 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #c3 e5 f6 j10
        
        else if (c3 == TRUE & e5 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #c3 e5 g7 h8
        
        else if (c3 == TRUE & e5 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #c3 e5 g7 i9
        
        else if (c3 == TRUE & e5 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #c3 e5 g7 j10
        
        else if (c3 == TRUE & e5 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #c3 e5 h8 i9
        
        else if (c3 == TRUE & e5 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #c3 e5 h8 j10
        
        else if (c3 == TRUE & e5 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #c3 e5 i9 j10
        
        else if (c3 == TRUE & e5 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #c3 f6 g7 h8
        
        else if (c3 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #c3 f6 g7 i9
        
        else if (c3 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #c3 f6 g7 j10
        
        else if (c3 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #c3 f6 h8 i9
        
        else if (c3 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #c3 f6 h8 j10
        
        else if (c3 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #c3 f6 i9 j10
        
        else if (c3 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #c3 g7 h8 i9
        
        else if (c3 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #c3 g7 h8 j10
        
        else if (c3 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #c3 g7 i9 j10
        
        else if (c3 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #c3 h8 i9 j10
        
        else if (c3 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #d4 e5 f6 g7
        
        else if (d4 == TRUE & e5 == TRUE & f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        #d4 e5 f6 h8
        
        else if (d4 == TRUE & e5 == TRUE & f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        #d4 e5 f6 i9
        
        else if (d4 == TRUE & e5 == TRUE & f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        #d4 e5 f6 j10
        
        else if (d4 == TRUE & e5 == TRUE & f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        #d4 e5 g7 h8
        
        else if (d4 == TRUE & e5 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #d4 e5 g7 i9
        
        else if (d4 == TRUE & e5 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #d4 e5 g7 j10
        
        else if (d4 == TRUE & e5 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #d4 e5 h8 i9
        
        else if (d4 == TRUE & e5 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #d4 e5 h8 j10
        
        else if (d4 == TRUE & e5 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #d4 e5 i9 j10
        
        else if (d4 == TRUE & e5 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #d4 f6 g7 h8
        
        else if (d4 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #d4 f6 g7 i9
        
        else if (d4 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #d4 f6 g7 j10
        
        else if (d4 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #d4 f6 h8 i9
        
        else if (d4 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #d4 f6 h8 j10
        
        else if (d4 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #d4 f6 i9 j10
        
        else if (d4 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #d4 g7 h8 i9
        
        else if (d4 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #d4 g7 h8 j10
        
        else if (d4 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #d4 g7 i9 j10
        
        else if (d4 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #d4 h8 i9 j10
        
        else if (d4 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #e5 f6 g7 h8
        
        else if (e5 == TRUE & f6 == TRUE & g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        #e5 f6 g7 i9
        
        else if (e5 == TRUE & f6 == TRUE & g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        #e5 f6 g7 j10
        
        else if (e5 == TRUE & f6 == TRUE & g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #e5 f6 h8 i9
        
        else if (e5 == TRUE & f6 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #e5 f6 h8 j10
        
        else if (e5 == TRUE & f6 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #e5 f6 i9 j10
        
        else if (e5 == TRUE & f6 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #e5 g7 h8 i9
        
        else if (e5 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #e5 g7 h8 j10
        
        else if (e5 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #e5 g7 i9 j10
        
        else if (e5 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #e5 h8 i9 j10
        
        else if (e5 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #f6 g7 h8 i9
        
        else if (f6 == TRUE & g7 == TRUE & h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        #f6 g7 h8 j10
        
        else if (f6 == TRUE & g7 == TRUE & h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        #f6 g7 i9 j10
        
        else if (f6 == TRUE & g7 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #f6 h8 i9 j10
        
        else if (f6 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        #g7 h8 i9 j10
        
        else if (g7 == TRUE & h8 == TRUE & i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        
        #a1 b2 c3
        
        
        else if (a1 == TRUE & b2 == TRUE & c3 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3)
        }
        
        
        
        #a1 b2 d4
        
        else if (a1 == TRUE & b2 == TRUE & d4 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4)
        }
        
        
        #a1 b2 e5
        
        else if (a1 == TRUE & b2 == TRUE & e5 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5)
        }
        
        
        #a1 b2 f6
        
        
        else if (a1 == TRUE & b2 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6)
        }
        
        
        #a1 b2 g7
        
        else if (a1 == TRUE & b2 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7)
        }
        
        
        #a1 b2 h8
        
        else if (a1 == TRUE & b2 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8)
        }
        
        
        #a1 b2 i9
        
        else if (a1 == TRUE & b2 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_9)
        }
        
        
        #a1 b2 j10
        
        else if (a1 == TRUE & b2 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_10)
        }
        
        #a1 c3 d4
        
        else if (a1 == TRUE & c3 == TRUE & d4 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4)
        }
        
        
        #a1 c3 e5
        
        else if (a1 == TRUE & c3 == TRUE & e5 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5)
        }
        
        
        #a1 c3 f6
        
        else if (a1 == TRUE & c3 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6)
        }
        
        
        #a1 c3 g7
        
        else if (a1 == TRUE & c3 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7)
        }
        
        #a1 c3 h8
        
        else if (a1 == TRUE & c3 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8)
        }
        
        
        #a1 c3 i9
        
        else if (a1 == TRUE & c3 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9)
        }
        
        
        #a1 c3 j10
        
        else if (a1 == TRUE & c3 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_10)
        }
        
        
        #a1 d4 e5
        
        else if (a1 == TRUE & d4 == TRUE & e5 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
        }
        
        
        #a1 d4 f6
        
        else if (a1 == TRUE & d4 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
        }
        
        
        #a1 d4 g7
        
        else if (a1 == TRUE & d4 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
        }
        
        
        #a1 d4 h8
        
        else if (a1 == TRUE & d4 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
        }
        
        
        #a1 d4 i9
        
        else if (a1 == TRUE & d4 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
        }
        
        
        #a1 d4 j10
        
        else if (a1 == TRUE & d4 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
        }
        
        
        #a1 e5 f6
        
        else if (a1 == TRUE & e5 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        
        #a1 e5 g7
        
        else if (a1 == TRUE & e5 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        
        #a1 e5 h8
        
        else if (a1 == TRUE & e5 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        
        #a1 e5 i9
        
        else if (a1 == TRUE & e5 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        
        #a1 e5 j10
        
        else if (a1 == TRUE & e5 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        #a1 f6 g7
        
        else if (a1 == TRUE & f6 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        
        #a1 f6 h8
        
        else if (a1 == TRUE & f6 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        
        #a1 f6 i9
        
        else if (a1 == TRUE & f6 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        
        #a1 f6 j10
        
        else if (a1 == TRUE & f6 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        
        #a1 g7 h8
        
        else if (a1 == TRUE & g7 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        
        #a1 g7 i9
        
        else if (a1 == TRUE & g7 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        
        #a1 g7 j10
        
        else if (a1 == TRUE & g7 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        #a1 h8 i9
        
        else if (a1 == TRUE & h8 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        
        #a1 h8 j10
        
        else if (a1 == TRUE & h8 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        
        #a1 i9 j10
        
        else if (a1 == TRUE & i9 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        #b2 c3 d4
        
        else if (b2 == TRUE & c3 == TRUE & d4 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4)
        }
        
        
        #b2 c3 e5
        
        else if (b2 == TRUE & c3 == TRUE & e5 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5)
        }
        
        
        #b2 c3 f6
        
        else if (b2 == TRUE & c3 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6)
        }
        
        
        #b2 c3 g7
        
        else if (b2 == TRUE & c3 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7)
        }
        
        
        #b2 c3 h8
        
        else if (b2 == TRUE & c3 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8)
        }
        
        
        #b2 c3 i9
        
        else if (b2 == TRUE & c3 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9)
        }
        
        
        #b2 c3 j10
        
        else if (b2 == TRUE & c3 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_10)
        }
        
        
        #b2 d4 e5
        
        else if (b2 == TRUE & d4 == TRUE & e5 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
        }
        
        
        #b2 d4 f6
        
        else if (b2 == TRUE & d4 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
        }
        
        
        #b2 d4 g7
        
        else if (b2 == TRUE & d4 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
        }
        
        
        #b2 d4 h8
        
        else if (b2 == TRUE & d4 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
        }
        
        
        #b2 d4 i9
        
        else if (b2 == TRUE & d4 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
        }
        
        
        #b2 d4 j10
        
        else if (b2 == TRUE & d4 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
        }
        
        
        #b2 e5 f6
        
        else if (b2 == TRUE & e5 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        
        #b2 e5 g7
        
        else if (b2 == TRUE & e5 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        
        #b2 e5 h8
        
        else if (b2 == TRUE & e5 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        
        #b2 e5 i9
        
        else if (b2 == TRUE & e5 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        
        #b2 e5 j10
        
        else if (b2 == TRUE & e5 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        
        #b2 f6 g7
        
        else if (b2 == TRUE & f6 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        
        #b2 f6 h8
        
        else if (b2 == TRUE & f6 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        
        #b2 f6 i9
        
        else if (b2 == TRUE & f6 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        
        #b2 f6 j10
        
        else if (b2 == TRUE & f6 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        
        #b2 g7 h8
        
        else if (b2 == TRUE & g7 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        
        #b2 g7 i9
        
        else if (b2 == TRUE & g7 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        
        #b2 g7 j10
        
        else if (b2 == TRUE & g7 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        
        #b2 h8 i9
        
        else if (b2 == TRUE & h8 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        
        #b2 h8 j10
        
        else if (b2 == TRUE & h8 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        
        #b2 i9 j10
        
        else if (b2 == TRUE & i9 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        #c3 d4 e5
        
        else if (c3 == TRUE & d4 == TRUE & e5 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
        }
        
        
        #c3 d4 f6
        
        else if (c3 == TRUE & d4 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
        }
        
        
        #c3 d4 g7
        
        else if (c3 == TRUE & d4 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
        }
        
        
        #c3 d4 h8
        
        else if (c3 == TRUE & d4 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
        }
        
        
        #c3 d4 i9
        
        else if (c3 == TRUE & d4 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
        }
        
        
        #c3 d4 j10
        
        else if (c3 == TRUE & d4 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
        }
        
        
        #c3 e5 f6
        
        else if (c3 == TRUE & e5 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        
        #c3 e5 g7
        
        else if (c3 == TRUE & e5 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        
        #c3 e5 h8
        
        else if (c3 == TRUE & e5 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        
        #c3 e5 i9
        
        else if (c3 == TRUE & e5 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        
        #c3 e5 j10
        
        else if (c3 == TRUE & e5 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        
        #c3 f6 g7
        
        else if (c3 == TRUE & f6 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        
        #c3 f6 h8
        
        else if (c3 == TRUE & f6 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        
        #c3 f6 i9
        
        else if (c3 == TRUE & f6 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        
        #c3 f6 j10
        
        else if (c3 == TRUE & f6 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        
        #c3 g7 h8
        
        else if (c3 == TRUE & g7 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        
        #c3 g7 i9
        
        else if (c3 == TRUE & g7 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        
        #c3 g7 j10
        
        else if (c3 == TRUE & g7 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        
        #c3 h8 i9
        
        else if (c3 == TRUE & h8 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        
        #c3 h8 j10
        
        else if (c3 == TRUE & h8 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        
        #c3 i9 j10
        
        else if (c3 == TRUE & i9 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        #d4 e5 f6
        
        else if (d4 == TRUE & e5 == TRUE & f6 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        
        #d4 e5 g7
        
        else if (d4 == TRUE & e5 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        
        #d4 e5 h8
        
        else if (d4 == TRUE & e5 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        
        #d4 e5 i9
        
        else if (d4 == TRUE & e5 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        
        #d4 e5 j10
        
        else if (d4 == TRUE & e5 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        
        #d4 f6 g7
        
        else if (d4 == TRUE & f6 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        
        #d4 f6 h8
        
        else if (d4 == TRUE & f6 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        
        #d4 f6 i9
        
        else if (d4 == TRUE & f6 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        
        #d4 f6 j10
        
        else if (d4 == TRUE & f6 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        
        #d4 g7 h8
        
        else if (d4 == TRUE & g7 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        
        #d4 g7 i9
        
        else if (d4 == TRUE & g7 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        
        #d4 g7 j10
        
        else if (d4 == TRUE & g7 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        
        #d4 h8 i9
        
        else if (d4 == TRUE & h8 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        
        #d4 h8 j10
        
        else if (d4 == TRUE & h8 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        
        #d4 i9 j10
        
        else if (d4 == TRUE & i9 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        #e5 f6 g7
        
        else if (e5 == TRUE & f6 == TRUE & g7 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        
        #e5 f6 h8
        
        else if (e5 == TRUE & f6 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        
        #e5 f6 i9
        
        else if (e5 == TRUE & f6 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        
        #e5 f6 j10
        
        else if (e5 == TRUE & f6 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        
        #e5 g7 h8
        
        else if (e5 == TRUE & g7 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        
        #e5 g7 i9
        
        else if (e5 == TRUE & g7 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        
        #e5 g7 j10
        
        else if (e5 == TRUE & g7 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        
        #e5 h8 i9
        
        else if (e5 == TRUE & h8 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        
        #e5 h8 j10
        
        else if (e5 == TRUE & h8 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        
        #e5 i9 j10
        
        else if (e5 == TRUE & i9 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        #f6 g7 h8
        
        else if (f6 == TRUE & g7 == TRUE & h8 == TRUE )
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        
        #f6 g7 i9
        
        else if (f6 == TRUE & g7 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        
        #f6 g7 j10
        
        else if (f6 == TRUE & g7 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        
        #f6 h8 i9
        
        else if (f6 == TRUE & h8 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        
        #f6 h8 j10
        
        else if (f6 == TRUE & h8 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        
        #f6 i9 j10
        
        else if (f6 == TRUE & i9 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        #g7 h8 i9
        
        else if (g7 == TRUE & h8 == TRUE & i9 == TRUE )
        {
          outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        
        #g7 h8 j10
        
        else if (g7 == TRUE & h8 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        
        #g7 i9 j10
        
        else if (g7 == TRUE & i9 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        #h8 i9 j10
        
        else if (h8 == TRUE & i9 == TRUE & j10 == TRUE )
        {
          outliers_SMALLCLUSTER_8  <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        #a1 b2
        
        
        else if (a1 == TRUE & b2 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_2 <-which(cutree_wektor == 2)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_2)
        }
        
        
        #a1 c3
        
        
        else if (a1 == TRUE & c3 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_3)
        }
        
        
        
        #a1 d4
        
        else if (a1 == TRUE & d4 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_4)
        }
        
        
        #a1 e5
        
        else if (a1 == TRUE & e5 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_5)
        }
        
        
        #a1 f6
        
        else if (a1 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_6)
        }
        
        
        #a1 g7
        
        else if (a1 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_7)
        }
        
        
        #a1 h8
        
        else if (a1 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_8)
        }
        
        
        #a1 i9
        
        else if (a1 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_9)
        }
        
        
        #a1 j10
        
        else if (a1 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1,outliers_SMALLCLUSTER_10)
        }
        
        
        #b2 c3
        
        else if (b2 == TRUE & c3 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_3 <-which(cutree_wektor == 3)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_3)
        }
        
        
        #b2 d4
        
        else if (b2 == TRUE & d4 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_4)
        }
        
        
        #b2 e5
        
        else if (b2 == TRUE & e5 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_5)
        }
        
        
        #b2 f6
        
        else if (b2 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_6)
        }
        
        
        #b2 g7
        
        else if (b2 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_7)
        }
        
        
        #b2 h8
        
        else if (b2 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_8)
        }
        
        
        #b2 i9
        
        else if (b2 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_9)
        }
        
        
        #b2 j10
        
        else if (b2 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2,outliers_SMALLCLUSTER_10)
        }
        
        
        #c3 d4
        
        else if (c3 == TRUE & d4 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_4 <-which(cutree_wektor == 4)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_4)
        }
        
        
        #c3 e5
        
        else if (c3 == TRUE & e5 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_5)
        }
        
        
        #c3 f6
        
        else if (c3 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_6)
        }
        
        
        #c3 g7
        
        else if (c3 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_7)
        }
        
        
        #c3 h8
        
        else if (c3 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_8)
        }
        
        
        #c3 i9
        
        else if (c3 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_9)
        }
        
        
        #c3 j10
        
        else if (c3 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3,outliers_SMALLCLUSTER_10)
        }
        
        
        #d4 e5
        
        else if (d4 == TRUE & e5 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_5 <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_5)
        }
        
        
        #d4 f6
        
        else if (d4 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_6)
        }
        
        
        #d4 g7
        
        else if (d4 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_7)
        }
        
        
        #d4 h8
        
        else if (d4 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_8)
        }
        
        
        #d4 i9
        
        else if (d4 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_9)
        }
        
        
        #d4 j10
        
        else if (d4 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4,outliers_SMALLCLUSTER_10)
        }
        
        
        #e5 f6
        
        else if (e5 == TRUE & f6 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_6 <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_6)
        }
        
        
        #e5 g7
        
        else if (e5 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_7)
        }
        
        
        #e5 h8
        
        else if (e5 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_8)
        }
        
        
        #e5 i9
        
        else if (e5 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_9)
        }
        
        
        #e5 j10
        
        else if (e5 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5,outliers_SMALLCLUSTER_10)
        }
        
        
        #f6 g7
        
        else if (f6 == TRUE & g7 == TRUE)
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_7 <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_7)
        }
        
        
        #f6 h8
        
        else if (f6 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_8)
        }
        
        
        #f6 i9
        
        else if (f6 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_9)
        }
        
        
        #f6 j10
        
        else if (f6 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6,outliers_SMALLCLUSTER_10)
        }
        
        
        #g7 h8
        
        else if (g7 == TRUE & h8 == TRUE)
        {
          outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_8 <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_8)
        }
        
        
        #g7 i9
        
        else if (g7 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_9)
        }
        
        
        #g7 j10
        
        else if (g7 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7,outliers_SMALLCLUSTER_10)
        }
        
        
        #h8 i9
        
        else if (h8 == TRUE & i9 == TRUE)
        {
          outliers_SMALLCLUSTER_8  <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_9 <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_9)
        }
        
        
        #h8 j10
        
        else if (h8 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_8  <-which(cutree_wektor == 8)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_8,outliers_SMALLCLUSTER_10)
        }
        
        
        #i9 j10
        
        else if (i9 == TRUE & j10 == TRUE)
        {
          outliers_SMALLCLUSTER_9  <-which(cutree_wektor == 9)
          outliers_SMALLCLUSTER_10 <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_9,outliers_SMALLCLUSTER_10)
        }
        
        
        #a1
        
        else if (a1 == TRUE)
        {
          outliers_SMALLCLUSTER_1  <-which(cutree_wektor == 1)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_1)
        }
        
        #b2
        
        else if (b2 == TRUE)
        {
          outliers_SMALLCLUSTER_2  <-which(cutree_wektor == 2)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_2)
        }
        
        #c3
        
        else if (c3 == TRUE)
        {
          outliers_SMALLCLUSTER_3  <-which(cutree_wektor == 3)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_3)
        }
        
        #d4
        
        else if (d4 == TRUE)
        {
          outliers_SMALLCLUSTER_4  <-which(cutree_wektor == 4)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_4)
        }
        
        #e5
        
        else if (e5 == TRUE)
        {
          outliers_SMALLCLUSTER_5  <-which(cutree_wektor == 5)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_5)
        }
        
        #f6
        
        else if (f6 == TRUE)
        {
          outliers_SMALLCLUSTER_6  <-which(cutree_wektor == 6)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_6)
        }
        
        #g7
        
        else if (g7 == TRUE)
        {
          outliers_SMALLCLUSTER_7  <-which(cutree_wektor == 7)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_7)
        }
        
        #h8
        
        else if (h8 == TRUE)
        {
          outliers_SMALLCLUSTER_8  <-which(cutree_wektor == 8)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_8)
        }
        
        #i9
        
        else if (i9 == TRUE)
        {
          outliers_SMALLCLUSTER_9  <-which(cutree_wektor == 9)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_9)
        }
        
        
        #j10
        
        else if (j10 == TRUE)
        {
          outliers_SMALLCLUSTER_10  <-which(cutree_wektor == 10)
          
          outliers_SMALLCLUSTER <- c(outliers_SMALLCLUSTER_10)
        }
        
        else
        {
          outliers_SMALLCLUSTER <- "t"
        }
        
################################################################################      
        
        if ( outliers_SMALLCLUSTER == "t")
        {    wyniki[12,3] <- "not small clusters"}
        else
        {wyniki[12,3] <- toString(outliers_SMALLCLUSTER)}
        
        if (outliers_SMALLCLUSTER == "t")
        { data_inf4 <- t }
        else
        {
          xx <- outliers_SMALLCLUSTER
          data_inf4 <- t[-xx, ]
        }   
        
        gower_mat1 <- gower.dist (data_inf4 ,  rngs = NA ,  KR.corr = TRUE ,  var.weights  =  NULL)
        gower_mat1[is.nan(gower_mat1)] <- 0
        gower_mat <-as.dist(gower_mat1) 
        clusters <- hclust(gower_mat, method = input$metoda)   
        
        #plot(clusters, hang = -1, main = "Dendrogram prezentujący podział zbioru danych na grupy po usunięciu reguł nietypowych - metoda SMALLCLUSTER")

        plot(fviz_dend(clusters, k = input$liczbaGrup[1],font.main = 20, main = paste("Dendrogram po usunięciu reguł nietypowych - metoda SMALLCLUSTER, dataset = ",collection_40), 
                       xlab = paste("Metoda łączenia grup:",input$metoda,"|", "Miara odległości między regułami: Gower " ),cex = 0.5,   # label size
                       k_colors = c("#2E9FDF", "#FC4E07", "#07fc4e", "#3339FF", "#00AFBB", "#FF33F0", "#E7B800", "#C0C0C0","#800080", "#008080", "#808000", "#2F4F4F", "#8B0000", "#FF1493","#BDB76B", "#FFD700", "#483D8B", "#00FF7F","#00CED1","#0000CD" ),
                       color_labels_by_k = TRUE, # color labels by groups
                       rect = FALSE # Add rectangle around groups 
        ))
        
      }
      
      
      
      
    })
    

  })
  
  
}

