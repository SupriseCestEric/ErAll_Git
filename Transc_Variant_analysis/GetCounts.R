

makeMatrix<- function(firstRow,fileList){
  Matrix1<-firstRow
  for (item in fileList) {
    tmp<-read.table(item, header=TRUE)
    Matrix1<-cbind(Matrix1,tmp[,4])
  }
  return(Matrix1)
}

main<-function(){
  #
  #GetCounts is a function to rip the 4th column from kallisto output (est. counts) and generates a matrix with all patients
  #
  #
  #Usage: Rscript GetCounts.R path_to_wdir Output_matrix_name
  #
  #
  #Source code may be modified to grab other columns, and caution must be taken when working on windows
  #Windows users may type paths in unix style "C:/Users/Eric" instead of "C:\Users\Eric".
  #
  #watch out for other .tsv files, make sure none are in the working dir apart from the kallisto output files.
  #This function grabs the .tsv files recursively, so clear any subdirectories of unwanted .tsv files as well.
  
  #get the correct path for wdir and output matrix name. 
  args<-commandArgs(trailingOnly = TRUE)
  
  #setwd
  folder=args[1];
  setwd(folder)
  #get name
  out_name=args[2]
  
  #list files in wdir recursively to get all of kallisto's tsv files
  filenames<-list.files(pattern = "*.tsv", full.names=TRUE, recursive = TRUE)
  
  #get only the first column of transcript IDs in refseq format
  gene_names <-read.table(filenames[1], header=TRUE)
  gene_names <-as.character(gene_names$target_id)
  
  #loop over files and grab 4th column (est. counts), 
  gene_names <- makeMatrix(gene_names,filenames)
  gene_names<-as.data.frame(gene_names)
  #Replace rownames with 1st column, then delete it. change head to filename
  rownames(gene_names)=gene_names[,1] #set first column as row names
  gene_names=gene_names[,-1] #remove first column that is useless now
  colnames(gene_names)<-filenames
  
  #Output matrix
  write.table(gene_names, file=out_name)
}

main()

