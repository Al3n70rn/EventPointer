#' EventPointer
#' 
#' Performs statistical analysis to obtain  significant alternative splicing events
#' @param Design Design matrix for the experiment
#' @param Contrast Contrast matrix for the experiment
#' @param affy Preprocessed aroma.affymetrix variable
#' @param array Type of array used in the experiment (HTA or Hjay)
#' @param Filter Binary variable to filter events according to expression
#' @param Qn Quantile used to filter the events (Bounded between 0-1, Q1 would be 0.25)
#' @return Data frame with statistical values for each event
#' @examples
#' \dontrun{
#'          library(EventPointer)
#' 
#'          setwd("J:/EventPointer/")
#'          verbose <- Arguments$getVerbose(-8);
#'          timestampOn(verbose);
#'          projectName <- "SRSF1"
#'          chipType <- "HTA_AS"
#'          cdfGFile <- "HTA_AS,r"
#'          cdfG <- AffymetrixCdfFile$byChipType(cdfGFile)
#'          cs <- AffymetrixCelSet$byName(projectName, cdf=cdfG)
#'          bc <- NormExpBackgroundCorrection(cs, method="mle", tag=c("*","r11"));
#'          csBC <- process(bc,verbose=verbose,ram=20);
#'          qn <- QuantileNormalization(csBC, typesToUpdate="pm");
#'          csN <- process(qn,verbose=verbose,ram=20);
#'          plmEx <- ExonRmaPlm(csN, mergeGroups=FALSE)
#'          fit(plmEx, verbose=verbose)
#'          cesEx <- getChipEffectSet(plmEx)
#'
#'          load("J:/EventPointer/DyCmatrix.Rdata")
#'          Events<-EventPointer(Design=Dmatrix,Contrast=Cmatrix,affy=cesEx,array="HTA",Filter=T)
#'        }
#' @export 
EventPointer<-function(Design,Contrast,affy,array,Filter=F,Qn=0.25)
{
  # Statistical analysis for detecting alternative splicing using affymetrix microarrays
  # and limma framework

  # verbose <- Arguments$getVerbose(-8)
  # cat(verbose," Analysis Started")
  # Hora<-Sys.time()
  cat(format(Sys.time(), "%X")," Running EventPointer")

  # Extract the summarized intensities for each of the probesets in the cdf.
  # The probesets of the cdf correspond to each of the paths (Path 1, Path 2 and Reference)
  # for each of the events detectable by the array.
  ExFit <- extractDataFrame(affy, addNames=TRUE)

  # Obtain the number of columns of the ExFit variable
  ncc<-ncol(ExFit)

  # Sanity check
  # All the events and probesets must be ordered in the following way:
  # Event_1_Ref
  # Event_1_P1
  # Event_1_P2
  # Event_2_Ref...

  # Giving each Paths an number indicator then we can order them by using it
  ExFit[ExFit[,2]=="_Ref",4] <- 1
  ExFit[ExFit[,2]=="_P1",4] <- 2
  ExFit[ExFit[,2]=="_P2",4] <- 3

  # The Probesets are ordered first by unit and then by group, with units being
  # the events and groups the paths
  NuevoOrden <- order(ExFit[,3],ExFit[,4])
  ExFit <- ExFit[NuevoOrden,]


  # AuxM is the Auxiliary Matrix used in the statistical analysis,
  # the rows of this matrix are used to indicate to which path each
  # of the probesets belong to.

  #     R   P1  P2
  #     1   0   0
  #     1   1   0
  #     1   1   1

  AuxM<-matrix(c(1,0,0,1,1,0,1,1,1),nrow=3,byrow=T)

  # With a Kronecker product, each of the elements of the Design Matrix
  # are replaced by AuxM.

  D <- kronecker(Design, AuxM)

  # We just keep the intensity values for the probesets in each sample
  A<-ExFit[,6:ncc]

  # The expression matrix (Ymat) must be in a different arrangement than ExFit.
  # ExFit has for rows the different paths for each event and for columns the samples
  # of the experiment. However, Ymat must have for rows each of the samples and the paths
  # it belongs to and for columns each of the events.

  #                       Original Matrix                                               Required Ordered Matrix
  #
  #                 Sample 1    Sample 2    Sample 3                                  Event 1     Event 2     Event 3
  #   Event 1_Ref                                                    Sample 1 _ Ref
  #   Event 1_P1                                                     Sampel 1 _ P1
  #   Event 1_P2                                                     Sample 1 _ P2
  #   Event 2_Ref                                            .       Sample 2 _ Ref
  #   Event 2_P1                                    ............     Sample 2 _ P1
  #   Event 2_P2                                    .............    Sample 2 _ P2
  #   Event 3_Ref                                   ............     Sample 3 _ Ref
  #   Event 3_P1                                             .       Sample 3 _ P1
  #   Event 3_P2                                                     Sample 3 _ P2


  II <- lapply(seq(1,nrow(A),3),function(x) rep(x+0:2,ncol(A)))
  JJ <- rep(rep(1:ncol(A),each=3),length(unique(ExFit[,1])))
  B <- A[cbind(unlist(II),JJ)]
  Ymat <- t(matrix(B,nrow=length(unique(ExFit[,1])),byrow=TRUE))
  colnames(Ymat)<-unique(ExFit[,1])
  rownames(Ymat)<-paste(rep(colnames(ExFit)[6:ncc],each=3),"_",c("Ref","P1","P2"),sep="")

  # Transform the matrix values to log2
  Ymat<-log2(Ymat)

  ### Start Statistical Analysis
  # cat(" \nPerforming Statistical Analysis with Limma...")


  # Linear model using Ymat and Design matrices
  fit <- lmFit(object=t(Ymat),design=D)

  # The contrasts we are interested in are the ones
  # related with each Path, and we apply a kronecker
  # product of the contrast matrix with the corresponding
  # vector for each Path (P1 =  1 1 0 ; P2 = 1 1 1)
  P1<-kronecker(Contrast,matrix(c(1,1,0),nrow=3))
  P2<-kronecker(Contrast,matrix(c(1,1,1),nrow=3))

  # C contains the vectors of the contrasts
  C<-cbind(P1,P2)

  # Compute estimated coefficients and standard errors for the given contrasts
  fit2 <- contrasts.fit(fit, C)

  # Empirical Bayesian Statistics
  fit2 <- eBayes(fit2)

  # Obtain the ranking of events for each of the contrasts
  T2<-topTable(fit2,coef=1,number=Inf)
  T3<-topTable(fit2,coef=2,number=Inf)

  # Both tables (T2 and T3) must be in the same order
  EvsIds<-rownames(T2)
  ii3<-match(EvsIds,rownames(T3))
  T3<-T3[ii3,]

  # One table is created by merging both, as both have the same column names, we rename the columns
  # from one of the tables with letters to avoid problems.
  colnames(T3)<-letters[1:ncol(T3)]
  T34_345 <- cbind(T2,T3)

  # By taking both pvalues, one from each contrast, we perform an Irwin Hall Sumarizaion
  # to obtain 1 pvalue for each of the events
  Values1<-IHsummarization(T34_345[,4],T34_345[,3],T34_345[,10],T34_345[,9])

  # Using the file EventsFound.txt we get all the information for each of the events
  # such as type of event, genomic position, gene, etc..

  if(array=="HTA")
  {
    data(HTA_Events)
    
    # Order the events from the txt in the same way as the tables from each contrast
    HTA_Events[,4]<-gsub(pattern = " ","",HTA_Events[,4])
    Ids<-paste(as.vector(HTA_Events[,2]),"_",as.vector(HTA_Events[,4]),sep="")
    index<-match(EvsIds,Ids)
    HTA_Events<-HTA_Events[index,]
    
    # Put all the information in one data frame
    Result<-data.frame(HTA_Events[,1],HTA_Events[,2],HTA_Events[,3],HTA_Events[,4],HTA_Events[,11],Values1$Tstats,Values1$Pvalues,stringsAsFactors=F)
    
    rm(HTA_Events)
    
  }else if(array=="Hjay")
  {
    data(Hjay_Events)
    
    # Order the events from the txt in the same way as the tables from each contrast
    
    Hjay_Events[,4]<-gsub(pattern = " ","",Hjay_Events[,4])
    Ids<-paste(as.vector(Hjay_Events[,2]),"_",as.vector(Hjay_Events[,4]),sep="")
    index<-match(EvsIds,Ids)
    Hjay_Events<-Hjay_Events[index,]
    
    # Put all the information in one data frame
    Result<-data.frame(Hjay_Events[,1],Hjay_Events[,2],Hjay_Events[,3],Hjay_Events[,4],Hjay_Events[,11],Values1$Tstats,Values1$Pvalues,stringsAsFactors=F)
    
    rm(Hjay_Events)
    
  }
  
  colnames(Result)<-c("HGNC Symbol","Ensembl Gene Id","Event Type","Event Number","Genomic Position",
                      "Splicing Z Value","Splicing Pvalue")
  rownames(Result)<-paste(Result[,2],"_",Result[,4],sep="")

  # Order the final dataframe according to the pvalue
  Result<-Result[order(Result[,7]),]

  # Remove events that are detected in abnormal genomic positions
  RowS<-nrow(Result)
  Pgen<-unlist(strsplit(as.vector(Result[,"Genomic Position"]),":"))[seq(1,RowS*2,by=2)]
  Hs<-grep("H",Pgen)
  GLs<-grep("GL",Pgen)
  Result<-Result[-c(Hs,GLs),]
  
  # Filter Events
  
  if(Filter==T)
  {
    Y<-extractDataFrame(affy,addNames=T)
    Y<-Y[,c(1,2,6:ncc)]
    
    Maxs<-rowMaxs(as.matrix(Y[,3:ncol(Y)]))
    Maxs_P1<-Maxs[seq(1,length(Maxs),3)]
    Maxs_P2<-Maxs[seq(2,length(Maxs),3)]
    Maxs_Ref<-Maxs[seq(3,length(Maxs),3)]
    
    Th<-quantile(Maxs_Ref,Qn)
    Filt<-which(Maxs_P1>Th & Maxs_P2>Th  & Maxs_Ref>Th)
    
    EE<-unique(Y[,1])
    
    FilterMatrix<-matrix(0,nrow=length(EE),ncol=2)
    FilterMatrix[,1]<-EE
    FilterMatrix[Filt,2]<-1
    
    Ixx<-match(rownames(Result),FilterMatrix[,1])
    FilterMatrix<-FilterMatrix[Ixx,]
    Jxx<-which(FilterMatrix[,2]==1)
    
    Result<-Result[Jxx,]

    
  }
  
  # Return the Result to the user
  # cat("Done")
  cat("\n",format(Sys.time(), "%X")," Analysis Completed")
  return(Result)
}
