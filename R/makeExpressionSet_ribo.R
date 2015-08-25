makeExpressionSet_ribo <- function(rnaCountFiles) {

# rnaExonCountFiles = list.files(file.path(dir,"data/rna"),pattern="*.txt",full.name=TRUE)
  options(stringsAsFactors=FALSE)
  
  require(parallel)
  
  # RIBO file(s)
  riboCountTable = read.table(riboCountFiles,sep="\t",header=T)
  colnames(riboCountTable) = toupper(colnames(riboCountTable))
  colnames(riboCountTable) = gsub(colnames(riboCountTable),
                                  pattern=".",replace="",fixed=TRUE)  
  
  fData_ribo = riboCountTable[,-2]
  fData_ribo = fData_ribo[order(fData_ribo$ENSGID),]
  
  
  # GENE LENGTH
  geneLen = data.frame(ENSGID=riboCountTable$ENSGID,len=riboCountTable$GENELENGTH)
  
  # make expression data sets
  require(Biobase)
  
  featureData = new("AnnotatedDataFrame",data=data.frame(ENSGID=geneLen$ENSGID,len=geneLen$len))
  
  rownames(fData_ribo) = fData_ribo[,1] 
  
  fData = fData_ribo[,-1]
  
  pData = data.frame(seqData=rep("ribo",(NCOL(fData_ribo)-1)),
                     species=c(rep(c( rep("chimp",5),rep("human",6) ,
                                    rep("rhesus",5)),1)),
                     celline=c(do.call(rbind,lapply(colnames(fData),
                                          function(x) {xx=gsub(x,pattern=".rna",replace="");
                                                       xx=gsub(xx,pattern=".ribo",replace="");xx}))) )
  rownames(pData) = colnames(fData)
  metaData = data.frame(labelDescription=c("ribo data",
                                           "human/chipmp/rhesus",
                                           "name of cell line"),
                        row.names=c("seqData","species","celline"))
  pData = new("AnnotatedDataFrame",data=pData,varMetadata=metaData)
  
  experimentData <- new("MIAME",title="Ribo-seq cross species comparison")
  eSet = ExpressionSet(assayData = as.matrix(fData),
                       phenoData=pData,
                       experimentData=experimentData)

  featureData(eSet) = featureData                  
  
  list(eSet=eSet,geneLen=geneLen)
}


