#' Protein Domains Enrichment Analysis
#'
#' Performs an enrichment analysis for the protein domains affected by the alternative splicing events
#' @param Events Dataframe returned by EventPointer. COntains information for all alternative splicing events
#' @return Result List with 4 elements: 1) Up Regulated Domains 2) Down Regulated Domains 3) Pvalues 4) Zvalues
#' @export
#' @examples
#' \dontrun{
#' DomainEnrichment(Events)}

DomainEnrichment<-function(Events)
{

  data(DomainsEvents)
  RnamesA<-paste(rownames(Events),"_PathA",sep="")
  RnamesB<-paste(rownames(Events),"_PathB",sep="")
  Rnames<-c(RnamesA,RnamesB)
  ix<-match(Rnames,colnames(DomainsEvents))
  DomainsEvents<-DomainsEvents[,ix]
  Zts<-Events[,"Splicing Z Value"]
  dummy <- runif(length(Zts))
  Dejar <- which(as.vector(abs(DomainsEvents %*% c(dummy, -dummy))>1e-12))
  DomainsEvents2<-DomainsEvents[Dejar,]
  DomainsEvents2@x<-rep(1,length(DomainsEvents2@x))
  Domains.z<-Wilcoxon.z(c(Zts,-Zts),Matrix::t(DomainsEvents2))*sqrt(2)
  names(Domains.z)<-rownames(DomainsEvents2)
  Pvalues<-2*pnorm(-abs(Domains.z))
  Values<-data.frame(Domain=rownames(DomainsEvents2),Pvalue=Pvalues,Zvalue=Domains.z)
  Values<-Values[order(abs(Values[,2])),]

  ZPos<-which(sign(Values[,3])>0)[1:20]
  ZNeg<-which(sign(Values[,3])<0)[1:20]

  Positive_Domains<-Values[ZPos,]
  Negative_Domains<-Values[ZNeg,]

  data(Pfam)

  iix<-match(Positive_Domains[,1],Pfam@data$id)
  Info<-Pfam@data$description[iix]
  Positive_Domains<-cbind(Positive_Domains,Info)

  iix<-match(Negative_Domains[,1],Pfam@data$id)
  Info<-Pfam@data$description[iix]
  Negative_Domains<-cbind(Negative_Domains,Info)

  Zvalues<-Values[,3]
  names(Zvalues)<-Values[,1]

  Pvalues<-Values[,2]
  names(Pvalues)<-Values[,1]

  Result<-list(UpRegulated_Domain=Positive_Domains,DownRegulated_Domain=Negative_Domains,Pvalues=Pvalues,Zvalues=Zvalues)

  return(Result)

}

