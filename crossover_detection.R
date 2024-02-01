#!/usr/bin/Rscript

# This script detect crossovers for each chromosome for a given parent (maternal or paternal) and outputs a bed file with chromosome, individual_id, and CO boundaries. A .txt file for each chromosome with SNP filtering information is also produced.

# The script requires scaffold boundary coordinates, duohmm corrected haplotype file, genotyping error file, and duohmm corrected sample file as input.

# To invoke:

# Rscript crossover_detection.R \
#   ${chromosome}_scaffold_boundary_file.txt \
#   ${chromosome}.maternal.shapeit.duohmm-corrected.haps \
#   ${chromosome}.maternal.shapeit.duohmm.generr \
#   ${chromosome}.maternal.shapeit.duohmm-corrected.sample \
#   mum \
#   ${chromosome} \
#   <output_prefix>;

# Further details and description are outlined in the supplementary methods of the paper.


library(Matrix)
rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}
inc <- function(x,step)
{
 eval.parent(substitute(x <- x + step))
}

opposite_haplotype<-function(hap){

  if(hap==hap1)
    return(hap2)
  else return(hap1)

}

#Change the values between positions p1 and p2 for offspring x to the filter name
apply_filter<-function(x,p1,p2,filter_name){

  for(k in p1:p2){
    eval.parent(substitute(x[k]<-filter_name))
  }

}

#Count the SNPs that we lost for all the filters
count_filter_snp<-function(x){

  nbSNP<-0
  for(i in 2:nrow(x)){
    for(j in 1:ncol(x)){

      if(substring(x[i,j], 1, 1)=='F')
        inc(nbSNP,1)

    }
  }
  return(nbSNP)

}

#Returns the end of a clean haplotype block(=no filtered positions and no scaff boundires)
find_clean_hap_end_pos<-function(x,pos){

  hap<-x[pos]
  while( x[pos]==hap ){
    if(pos==length(x))
      return(pos)
    inc(pos,1)
  }
  return(pos-1)
}

#Returns the end of a haplotype block(=with filtered position & scaff boundries)
find_hap_end_pos<-function(x,pos){

  hap<-x[pos]
  op_hap<-opposite_haplotype(hap)
  while(x[pos]!=op_hap & pos<length(x))
    inc(pos,1)
  while(x[pos]!=hap)
    inc(pos,-1)

  return(pos)
}

#Filters small clean haplotype blocks with low nb of SNPs
filter_cleanhap_blocks<-function(x,snp_limit,snp_filter_name){

  #For each offspring
  for(i in 2:nrow(x)){

    j<-2
    while(j<=ncol(x)){

      if(x[i,j]==hap1 | x[i,j]==hap2){#we are in a haplotype block

        #find the end of the clean hapotype block
        end_block<-find_hap_end_pos(x[i,],j)

        if(end_block-j+1<=snp_limit){#the block contains to few SNPs
          apply_filter(x[i,],j,end_block,snp_filter_name);
        }

        j<-end_block+1 #jump the block

      }else{#we are in a filtered position or a scaffold bound
        j<-j+1
      }

    }


  }
  return(x)

}




#Filters all the small haplotype blocks based on a min bp length
filter_gene_switch<-function(x,bp_limit,bp_filter_name){

  #For each offspring
  for(i in 2:nrow(x)){

    j<-2
    while(j<=ncol(x)){

      if(x[i,j]==hap1 | x[i,j]==hap2){#we are in a haplotype block

        #find the end of the clean hapotype block
        end_block<-find_hap_end_pos(x[i,],j)

        if( ( (as.numeric(x[1,end_block])-as.numeric(x[1,j]))<bp_limit ) & #the block is too short#
            (x[i,j-1]!=scaff_beg & x[i,end_block]!=scaff_end) ){ #and we are not at scaff boundries

          apply_filter(x[i,],j,end_block,bp_filter_name);

        }

        j<-end_block+1 #jump the block

      }else{#we are in a filtered position or a scaffold bound
        j<-j+1
      }

    }


  }
  return(x)

}

#Returns the end of a haplotype block(=with filtered position and scaff boundries)
hap_snps_details<-function(x,pos){

  hap<-x[pos]
  op_hap<-opposite_haplotype(hap)
  phased_snp_count<-0
  while(x[pos]!=op_hap & pos<length(x)){
    if(x[pos]==hap)
      inc(phased_snp_count,1)
    inc(pos,1)
  }
  while(x[pos]!=hap)
    inc(pos,-1)
  return(c(phased_snp_count,pos))
}

#Filter all the haplotype blocks that have too few phased SNPs
filter_sparse_blocks<-function(x,phased_snp_limit,unphased_switch_filter_name){

  #For each offspring
  for(i in 2:nrow(x)){

    j<-2
    while(j<=ncol(x)){

      if(x[i,j]==hap1 | x[i,j]==hap2){#we are in a haplotype block

        #find the number of phased SNPs and the end of the hapotype block
        snps_details<-hap_snps_details(x[i,],j)

        if( (snps_details[1]/(snps_details[2]-j+1))<phased_snp_limit){

          apply_filter(x[i,],j,snps_details[2],unphased_switch_filter_name)

        }

        j<-snps_details[2]+1 #jump the block

      }else{#we are in a filtered position or a scaffold bound
        j<-j+1
      }

    }




  }
  return(x)

}



#MAIN START

#values for the matrix
hap1<-"A"
hap2<-"B"

scaff_beg<-"b"
scaff_end<-"e"

generr_filter<-"F0"
small_switch_filter_name_snp<-"F1"
small_switch_filter_name_bp<-"F2"
unphased_switch_filter_name<-"F3"

#parameters
error_prob_threhold<-0.9
small_snp_limit<-50
small_bp_limit<-50000
phased_snp_limit<-0.6

#input files: - scaffold boundries for the chromosome (scaffName \t begPos \t endPos)
#             - haps file
#             - generr file
#             - sample file
#             - parent: dad or mum


argv <- commandArgs(TRUE)

scaffBound<-as.matrix(read.table(argv[1],sep="\t", header=F))
haps<-as.matrix(read.table(argv[2],sep=" ", header=F))
generrSNPs<-as.matrix(read.table(argv[3],sep="\t",header=T))
sample<-as.matrix(read.table(argv[4],sep=" ",header=T))
parent<-argv[5]
chromosome<-argv[6]

outputName<-argv[7]
bedOutputFile<-paste(outputName,".bed",sep="")
statOutputFile<-paste(outputName,".txt",sep="")
#pdfOutputFile<-paste(outputName,".pdf",sep="")
#plot name ".pdf"
#pdf(pdfOutputFile)

offsping_to_print<-10
#selecting the informative columns

if(parent=="dad"){
  info_cond<-"haps[,6] != haps[,7] & haps[,8] == haps[,9]"
  off_col_start<-10
  p_col_start<-6
}else if(parent=="mum"){
  info_cond<-"haps[,6] == haps[,7] & haps[,8] != haps[,9]"
  off_col_start<-11
  p_col_start<-8
}



SNP_positions<-subset(haps, eval(parse(text=info_cond)), select=c(3))

offsprings<-subset(haps,eval(parse(text=info_cond)),select=c(3,seq(off_col_start,ncol(haps), by=2)))
colnames(offsprings)<-c("pos",as.character(sample[4:length(sample[,2]),2]))

parent<-subset(haps, eval(parse(text=info_cond)), select=c(p_col_start:(p_col_start+1)))


#offsprings
#paste(parent[,1],collapse=' ');
#paste(parent[,2],collapse=' ');
#paste(offsprings[,offsping_to_print],collapse=' ');

#Transforming offsprings columns into haplotype
for (i in 2:length(offsprings[1,])) {
  offsprings[which(offsprings[,i]==parent[,1]),i]<-hap1;
  offsprings[which(offsprings[,i]==parent[,2]),i]<-hap2;
}
#paste(offsprings[,offsping_to_print],collapse=' ');
#Filter0: Filtering the SNPs with genotype error
for (i in 2:length(offsprings[1,])) {

    generr_pos<-generrSNPs[which((generrSNPs[,1]==as.character(colnames(offsprings)[i]))&(generrSNPs[,6]>error_prob_threhold)),5]
    offsprings[offsprings[,1] %in% generr_pos,i]<-generr_filter
}

#offsprings
#paste(offsprings[,offsping_to_print],collapse=' ');
#Adding scaffold boundires in each offspring column
scaff_b<-rep(scaff_beg,length(offsprings[1,])-1)
scaff_b_mat<-rep.row(scaff_b,length(scaffBound[,2]))
scaff_e<-rep(scaff_end,length(offsprings[1,])-1)
scaff_e_mat<-rep.row(scaff_e,length(scaffBound[,3]))

s_b<-cbind(scaffBound[,2],scaff_b_mat)
s_e<-cbind(scaffBound[,3],scaff_e_mat)

#s_b
colnames(s_b)<-c("pos",as.character(sample[4:length(sample[,2]),2]))
colnames(s_e)<-c("pos",as.character(sample[4:length(sample[,2]),2]))

options(warn=1)

total<-as.matrix(rbind(s_e,s_b,offsprings))
total_order<-total[order(as.numeric(total[,1])), ]


total_order_t<-t(total_order)
#nrow(total_order_t)
#ncol(total_order_t)


#Filter1: filtering clean haplotype blocks with a low number of SNPs = phasing/sequencing/mapping error
total_order_t<-filter_cleanhap_blocks(total_order_t,small_snp_limit,small_switch_filter_name_snp)

#Filter2: filtering short haplotype blocks = gene switch
total_order_t<-filter_gene_switch(total_order_t,small_bp_limit,small_switch_filter_name_bp)

#Filter3: filtering haplotype blocks with too few phased SNPs = sparse block
total_order_t<-filter_sparse_blocks(total_order_t,phased_snp_limit,unphased_switch_filter_name)


#print bed output : chr \t ind \t star \t stop

for(i in 2:nrow(total_order_t)){

  nb_of_XO<-0
  prevHap<-"0"
  for(j in 2:ncol(total_order_t)){
      if(total_order_t[i,j]==hap1 | total_order_t[i,j]==hap2){#we are in a haplotype block
          if(prevHap=="0")
            prevHap<-total_order_t[i,j]
          else if(prevHap!=total_order_t[i,j]){
            prevHap<-total_order_t[i,j]
            inc(nb_of_XO,1)
            write(paste(c(chromosome, colnames(offsprings)[i], total_order_t[1,beg], total_order_t[1,j]), collapse = "\t"), file=bedOutputFile,append=TRUE)
            #write(paste(c(chromosome, colnames(offsprings)[i], beg, j, total_order_t[1,beg] ,total_order_t[1,j]), collapse = "\t"), file=statOutputFile,append=TRUE)
            matrix_line=total_order_t[i,beg:j]
            write(paste(c(chromosome, colnames(offsprings)[i], beg, j, total_order_t[1,beg] ,total_order_t[1,j], matrix_line), collapse = "\t"), file=statOutputFile,append=TRUE)

          #  cat(sprintf("XO index: %i %i genome position: %s %s  \n",beg,j, total_order_t[1,beg] , total_order_t[1,j]))
          }
          beg<-j
      }
  }

#  cat(sprintf("offspring %i %s : %i \n", i, colnames(offsprings)[i], nb_of_XO))

}
#close(bedOutputFile)
#close(statOutputFile)


#Plotting

#matrix_color <- matrix(, nrow = nrow(total_order_t), ncol = ncol(total_order_t))

#for (i in 1:length(total_order_t[1,])) {
#  matrix_color[which(total_order_t[,i]==hap1),i]<-"red";
#  matrix_color[which(total_order_t[,i]==hap2),i]<-"blue";

#  matrix_color[which(total_order_t[,i]==generr_filter),i]<-"green";

#  matrix_color[which(total_order_t[,i]==small_switch_filter_name_snp),i]<-"yellow";

#  matrix_color[which(total_order_t[,i]==small_switch_filter_name_bp),i]<-"orange";

#  matrix_color[which(total_order_t[,i]==unphased_switch_filter_name),i]<-"purple1";

#  matrix_color[which(total_order_t[,i]==scaff_beg),i]<-"black";

#  matrix_color[which(total_order_t[,i]==scaff_end),i]<-"black";
#}

#matrix_color
#par(mar=c(4,5,3,1)+.1)
#chrSize<-length(total_order_t[1,])
#plot(NULL,xlab="Variants position (Mb)", xlim=c(0,chrSize),ylim=c(0,145),ylab="Offspings", yaxt="n",xaxt="n", cex.lab=2,cex.axis=2)
#for (i in 2:length(total_order_t[,1])) {
#  ypos<-rep(i,length(total_order_t[1,]))
#  xpos<-c(1:length(total_order_t[1,]))
  #print(xpos)
  #print(total_order_t[i,])
  #print(matrix_color[i,])
#  points(xpos,ypos,col=matrix_color[i,],pch=15, cex=0.5);
#}

#dev.off()
