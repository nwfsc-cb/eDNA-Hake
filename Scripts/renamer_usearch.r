usearch.rename.back<- function(fasta.in, fasta.out=NULL){
 require(stringr)
  require(seqinr)
   if(!is.null(fasta.out)){}else{
    fasta.out = file.path(dirname(fasta.in),str_replace(basename(fasta.in), ".fasta$", "_renamed.fasta"))
  }
  
  read.fasta(fasta.in) -> temp
str_remove_all(names(temp), ";size=[0-9]+;$") -> names(temp)

write.fasta(temp,names(temp),file.out = fasta.out )
  
  
}
