
# filter to take UMIs with certain number of reads and with certain proportion assigned to top BC
# find top BC for each cellBC-UMI combination
# filter to keep cells with certain proportion of UMIs having the same barcode
# output one VDJ barcode per cell along with confidence metrics
# Anne-Marie Lyne, Jun 2022

rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
wd = as.character(args[1])
n.const_v = as.numeric(args[2])
n.const_j = as.numeric(args[3])
samp.name = as.character(args[4])
min.read = as.numeric(args[5])
min.prop = as.numeric(as.character(args[6]))
min.prop.umi = as.numeric(as.character(args[7]))
min.umi = as.numeric(args[8])
min.prop.1umi = as.numeric(as.character(args[9]))

setwd(wd)

library(stringdist)
library(seqinr)
library(Biostrings)

#make relevant directory name
min.prop.str = gsub(".","",min.prop,fixed=T)
if(min.umi==1){
  min.prop.1umi.str = gsub(".","",min.prop.1umi,fixed=T)
  out.dir = paste0("output_",min.read,"_",min.prop.str,"_",min.umi,"_",min.prop.1umi.str)
}else{
  out.dir = paste0("output_",min.read,"_",min.prop.str,"_",min.umi)
}
dir.create(out.dir)

#names of input and output files
infile = paste0("agrep_10xbc_and_vbc_",samp.name)
outfile1 = paste0(out.dir,"/top_bc_per_UMI_VDJ_barcodes_",samp.name,"_",n.const_v,"_",n.const_j)
outfile2 = paste0(out.dir,"/top_bc_per_cell_VDJ_barcodes_",samp.name,"_",n.const_v,"_",n.const_j)


if(file.exists(paste0(infile,"_both.txt"))){
  nr = length(count.fields(paste0(infile,"_both.txt"), sep="\t"))
  if(nr>0){
    #read in files with CB, UB VDJ barcode etc
    maxfields = max(count.fields(paste0(infile,"_both.txt"), sep="\t"))
    data1 = read.table(paste0(infile,"_both.txt"), nrows = nr, comment.char="", 
                       colClasses = "character", fill=T, header=F, sep="\t")
    
    #get all 10x cell barcodes
    cb = unique(data1[,1])
    n.cb = length(cb)
    
    #count total number of unique sequences
    n.seq.ini = length(unique(data1[,3]))
    
    df1 = data.frame(Cell.bc = NULL, UMI = NULL, Con.seq = NULL,
                        Prop.top.bc = NULL, N.reads = NULL, stringsAsFactors = F)
    
    ##############################################################
    ### Filter UMIs based on no. reads and proportion for top BC
    ### Output top BC for each retained UMI
    ##############################################################
    
    #iterate through each BC UMI combination
    # count number of reads per UMI and proportion of reads for top barcode
    for(i in 1:n.cb){
      #take rows with same cell barcode
      rel.rows = data1[data1[,1]==cb[i],]
      #check how many UMIs this cell has and how many reads each has associated
      rel.umi = unique(rel.rows[,2])
      N.reads.umi = sapply(rel.umi, function(x) {sum(rel.rows[,2]==x)})
      #find top bc for each UMI
      top.bc = sapply(rel.umi, function(x) {
        rel.rows2 = rel.rows[rel.rows[,2]==x,]
        return(names(sort(table(rel.rows2[,3]),decreasing=T))[1])
      })
      #find proportion of reads for top barcode for each UMI
      max.prop = sapply(rel.umi, function(x) {
        rel.rows2 = rel.rows[rel.rows[,2]==x,]
        Nreads.max = sort(table(rel.rows2[,3]),decreasing=T)[1]
        return(Nreads.max/nrow(rel.rows2))
      })
      
      df.add = data.frame(Cell.bc = cb[i], UMI = rel.umi, Con.seq = top.bc,
                          Prop.top.bc = max.prop, N.reads = N.reads.umi, stringsAsFactors = F)
      df1 = rbind(df1, df.add)
    }
    
    #filter to take only UMIs with minimum min.read reads and min.pro proportion for top barcode
    # then filter cells to keep only those with min.umi UMIs
    logi = sapply(seq(nrow(df1)), function(x) {
      df1[x,5]>=min.read && df1[x,4]>=min.prop
    })
    df2 = df1[logi,]
    #count how many UMIs each cell has
    cb.filt = unique(df2$Cell.bc)
    Numi = sapply(cb.filt, function(x) {
      rel.rows2 = df2[df2$Cell.bc==x,]
      return(nrow(rel.rows2))
    })
    cb2 = cb.filt[Numi>=min.umi]
    df3 = df2[which(df2$Cell.bc %in% cb2),]
    
    #order by cell barcode
    df3 = df3[order(df3[,1]),]
    write.csv(df3, file = paste0(outfile1,"_both.csv"), quote = F, row.names = F)
    
    
    
    ##############################################################
    ### Filter cells based on consensus across UMIs
    ##############################################################
    
    #require at least min.prop.umi of barcodes across UMIs to agree
    all.cells = unique(df3[,1])
    
    #count number of cells and sequences
    n.cell.filt1 = length(all.cells)
    n.seq.filt1 = length(unique(df3[,3]))
    
    cells.to.keep = all.cells[sapply(all.cells, function(x) {
      rel.rows4 = df3[df3$Cell.bc==x,]
      #find top VDJ bc
      top.bc = names(sort(table(rel.rows4$Con.seq), decreasing = T))[1]
      return((sum(rel.rows4$Con.seq==top.bc)/length(rel.rows4$Con.seq))>=min.prop.umi)
    })]
    df4 = df3[which(df3$Cell.bc %in% cells.to.keep),]
    
    #get final VDJ barcode
    Con.seq = sapply(cells.to.keep, function(x) {
      rel.rows5 = df4[df4$Cell.bc==x,]
      return(names(sort(table(rel.rows5$Con.seq), decreasing = T))[1])
    })
    
    #compute some metrics to assess quality
    #number of UMIs for each cell before any filtering
    Numi.tot = sapply(cells.to.keep, function(x) {
      length(unique(data1[data1[,1]==x,2]))
    })
    
    #number of UMIs remaining after filtering
    Numi.filt = sapply(cells.to.keep, function(x) {
      nrow(df4[df4$Cell.bc==x,])
    })
    
    #proportion of unfiltered UMIs supporting final call
    UMI.prop = sapply(seq(length(Con.seq)), function(x) {
      #get all top bcs for this cell pre-filtering
      rel.rows6 = df1[df1[,1]==cells.to.keep[x],]
      sum(rel.rows6[,3]==Con.seq[x])/Numi.tot[x]
    })
    
    #proportion of filtered UMIs supporting final call
    UMI.prop.filt = sapply(seq(length(Con.seq)), function(x) {
      nrow(df4[df4$Cell.bc==cells.to.keep[x] & df4$Con.seq==Con.seq[x],])/Numi.filt[x]
    })
    
    #proportion of all unfiltered reads supporting final call
    read.prop = sapply(seq(length(cells.to.keep)), function(x) {
      nrow(data1[data1[,1]==cells.to.keep[x] & data1[,3]==Con.seq[x],])/
        nrow(data1[data1[,1]==cells.to.keep[x],])
    })
    
    df5 = data.frame(Cell.bc = cells.to.keep, Con.seq = Con.seq, Numi.tot = Numi.tot,
                     Numi.filt = Numi.filt, UMI.prop = UMI.prop, UMI.prop.filt = UMI.prop.filt, 
                     read.prop = read.prop, stringsAsFactors = F)
    
    ##############################################################
    ### If we are keeping cells with only 1 UMI, ensure the top BC  
    ### has a large enough proportion of the reads
    ##############################################################
    
    if(min.umi==1){
      #for cells with only one UMI, ensure top.bc is overall top BC in reads
      #and has at least min.prop.1umi proportion of total reads
      logi2 = sapply(seq(nrow(df5)), function(x) {
        if(df5$Numi.filt[x]>1){
          return(TRUE)
        }else{
          return(df5$read.prop[x]>=min.prop.1umi &&
                   names(sort(table(data1[data1[,1]==df5$Cell.bc[x],3]),decreasing=T))[1]==df5$Con.seq[x])
        }
      })
      df5 = df5[logi2,]
    }
    
    write.csv(df5, file = paste0(outfile2,"_both.csv"), quote = F, row.names = F)
  }
}

print(paste0("At UMI consensus step, ",length(cells.to.keep),
             " cells are retained out of ",length(all.cells)))
if(exists("logi2")){
  print(paste0("At 1UMI filtering stage, ",sum(logi2)," cells are retained out of ",
               length(logi2)))
}




