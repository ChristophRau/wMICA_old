#Inputs  
#input:  Rows are Samples, Columns are Genes
#output_prefix: prefix to affix to outputted files. 
#n.cores:  number of cores.  Defaults to total number minus 1  (detectCores())
#threshold:  lowest MINE value to send to ICMg.  Defaults to .3.  Lower equals 'better' but takes a LOT longer
#C:  Number of Modules in ICMg
#Seeds:  Which genes to seed.  Names need to be the same as in the input file
#Weight:  Weight to give each Seed.
#alpha:  for ICMg
#beta:  for ICMg
#B.num:  Number of Burnin rounds
#B.size:  Size of one burnin round
#S.num:  Number of Sample Rounds
#S.size:  Size of one sample round

'wMICA.Run'<-function(input,output_prefix="MICA",n.cores=detectCores()-1,threshold=.4,C=20,Seeds=c(),Weight=0,
                      alpha=10,beta=.01,B.num=10,B.size=10,S.num=10,S.size=10,MIC_ARRAY_FILE=NULL){
  
  gene_names=colnames(input)
  if(length(Seeds)>0){
  temp=match(Seeds[,1],gene_names)
  Seeds[,1]=temp
  Seeds=apply(Seeds,2,as.numeric)
  }
  input=apply(input,2,as.numeric)
  if(length(MIC_ARRAY_FILE)>0){
    print("Reading in MIC File.  Skipping MINERVA")
    MINE=read.csv(MIC_ARRAY_FILE,row.names=1)
  } else {
  print(paste0("Starting MINERVA with ",n.cores, " cores.  Please Wait."))
  MINE=mine(x=input,n.cores=n.cores)$MIC
  outname=paste0(output_prefix,"_MICARRAY.csv")
  print(paste0("saving MIC Array as ",outname))
  write.csv(MINE,file=outname)
  }
  outname=paste0(output_prefix,"_Links.csv")
  print(paste0("Creating Links File ",outname))
  outfile=file(outname,"w")
  firstRow="Gene1,Gene2,Strength"
  cat(firstRow,file=outfile,sep="\n")
  MINE_THRESH=MINE>threshold
  edges=0
  for(i in 1:(nrow(MINE)-1)){
    for(j in (i+1):nrow(MINE)){
      if(MINE_THRESH[i,j]){
        edges=edges+1
        outrow=paste(i,j,MINE[i,j],sep=",")
        cat(outrow,file=outfile,sep="\n")
      }
    }
  }
close(outfile)  
LINKS=read.csv(outname)
LINKS=apply(LINKS,2,as.numeric)
nodes=c(LINKS[,1],LINKS[,2])
nodes=length(table(nodes))

print(paste0("There are ",nodes," genes and ", edges," edges.  Starting Weighted ICMg."))
wICMg=ICMg.links.sampler.SeedsAndWeights(L = LINKS,C = C,Seeds = Seeds,Weight = Weight,
                                         alpha = alpha,beta = beta,B.num = B.num,
                                         B.size = B.size,S.num = S.num,S.size = S.size,
                                         C.boost = 1)
print("Calculating Memberships...")
wICMg$comp.memb=ICMg.get.comp.memberships(LINKS,wICMg)
outdata=wICMg$comp.memb
rownames(outdata)=c(1:C)
colnames(outdata)=gene_names      
outname=paste0(output_prefix,"_MODULEMEMBERSHIPS.csv")
print(paste0("Writing Output File ",outname))
write.csv(outdata,file=outname)
print("Done!")
}

