

#' Extracting transcript name and spanning reads for circRNAs
#'
#' This function makes a text file containing circle ids, transcript name and the spanning redas for each circRNAs
#' @param circ_out The output from the circRNA prediction tool CIRCexplorer
#' @param chimeric_out The junction file obtained from the chimeric alignment produced by STAR
#' @param output The name of the output folder
#'
#' @return The circle ids, transcript name and the spanning reads will be written in the file 'output' folder
#'
#' @importFrom utils write.table
#' @importFrom utils read.table
#' @export
#'
#' @examples
#'
#' #Loading an example output from the circRNA prediction tool CIRCexplorer
#' circ_out<-data(CIRCexplorer_output)
#' circ_out<-CIRCexplorer_output
#'
#' #Loading an example junction file obtained from the chimeric alignment
#' #produced by STAR
#' chimeric_out<-data(Chimeric.out.junction)
#' chimeric_out<-Chimeric.out.junction
#'
#' #A temporary directory is created as an output directory
#' output<-tempdir()
#'
#' #Extracting circle ids, transcript name and spanning reads
#' #of circRNAs. The output will be written in 'output' directory
#' extract_circle_ids(circ_out, chimeric_out, output)
#'

extract_circle_ids<-function(circ_out, chimeric_out, output){
  A<-paste(circ_out$V1, circ_out$V2, sep=":")
  circle_id_1<-paste(A, circ_out$V3, sep="|")
  indx1<-which(chimeric_out$V3=="+")
  indx1_chimericout<-chimeric_out[indx1,]
  A1<-paste(indx1_chimericout$V1, indx1_chimericout$V5, sep=":")
  junction_id1<-paste(A1, indx1_chimericout$V2-1, sep="|")
  junction_p<-cbind(junction_id1, indx1_chimericout)
  indx2<-which(chimeric_out$V3=="-")
  indx2_chimericout<-chimeric_out[indx2,]
  A2<-paste(indx2_chimericout$V1, indx2_chimericout$V2, sep=":")
  junction_id2<-paste(A2, indx2_chimericout$V5-1, sep="|")
  junction_n<-cbind(junction_id2, indx2_chimericout)
  c1<-1:14
  c2<-"V"
  colname<-paste0(c2, c1)
  colnames(junction_p)<-colname
  colnames(junction_n)<-colname
  junction<-rbind(junction_p,junction_n)
  for( i in 1:length(circle_id_1)){
    ind1<-which(junction$V1==circle_id_1[i])
    if(length(ind1)>0){
      reads<-paste(junction[ind1,11],collapse=",")
      circle_ids<-cbind(circle_id_1[i], as.character(circ_out$V16[i]), reads)
      write.table(circle_ids, file.path(output, "circle_ids.txt"), row.names=F, col.names=F, sep="\t", quote=F, append=T)
    }
  }
}


#' Extracting individual alignment file for each circRNA from chimeric alignment bam file
#'
#' This function generates alignment file for each circRNA from the chimeric alignment produced by STAR
#' @param circle_id A data frame containing circle ids, transcript name and spanning reads of circRNAs (obtained from function \code{\link[circRNAFull]{extract_circle_ids}})
#' @param bamfile A chimeric alignment bam file read as a data frame
#' @param outfolder The name of the output directory
#'
#' @return Individual alignment files for each circRNAs
#'
#' @importFrom Rsamtools scanBam
#' @export
#'
#' @examples
#' #loading an example circle_id file
#' data(circle_ID)
#' circle_id<-circle_ID
#'
#' \donttest{
#' \dontrun{
#' #Please upload your chimeric alignment bam file.Suppose you have the chimeric
#' #alignemnt bam file 'Chimeric.out.bam' in you working directory.
#' #Then run:
#' bam <- scanBam("Chimeric.out.bam")
#' bamfile <- as.data.frame(bam)
#' }
#' }
#'
#' #Creating an output directory
#' outfolder<-tempdir()
#'
#' \donttest{
#' \dontrun{
#' #Extracting individual alignment file for each circRNAs. The individual alignment file
#' #will be generated in the 'outfolder' directory.
#' extract_reads_from_bam<-function(circle_id, bamfile, outfolder)
#' }
#' }


extract_reads_from_bam<-function(circle_id, bamfile, outfolder){
  ##################################################################
  #############Start of function read_bam_as_dataframe()############
  ##################################################################
  read_bam_as_dataframe<-function(bamfile){
    bam <- scanBam(bamfile)
    bam_df <- as.data.frame(bam)
    return(bam_df)
  }
  ##################################################################
  #############End of function read_bam_as_dataframe()##############
  ##################################################################
  for(i in 1:nrow(circle_id)){
    circle_read<-unlist(strsplit(as.character(circle_id[i,3]), ","))
    v=NULL
    for (j in 1:length(circle_read)){
      v1=which(bamfile$qname==circle_read[j])
      v=c(v,v1)
    }
    extracted_reads<-data.frame(bamfile[v,])
    index<-which(extracted_reads$mapq<3)
    if (length(index)>0)extracted_reads_mapq<-extracted_reads[-index,] else extracted_reads_mapq<-extracted_reads
    circle_id_t<-chartr(":|", "__", circle_id[i,1])
    write.table(extracted_reads_mapq, file.path(outfolder, paste(circle_id_t, "_", length(circle_read),"reads.txt", sep="")), sep = "\t", row.names = FALSE, quote= FALSE)
  }
}





#' Extracting exon for the circRNAs
#'
#' This function extracts exons for the circRNAs from the CIGAR value of the chimeric alignment
#' @param bed_file A data frame containing the exon coordinated obtained from reference genome annotation
#' @param folder_name Location of the individual alignments (obtained from function \code{\link[circRNAFull]{extract_reads_from_bam}})
#' @param circle_ids A data frame containing the circle id, transcript name and spanning reads (obtained from function \code{\link[circRNAFull]{extract_circle_ids}})
#' @param output_folder The output folder name
#'
#' @return Exons for circRNAs
#' @export
#'
#' @examples
#' #Loading an example bed_file
#' data(exon_coordinate)
#' bed_file<-exon_coordinate
#'
#' #Example of folder_name containing the individual alignments
#' folder_name<-system.file('extdata', package = 'circRNAFull')
#'
#' #Loading an example circle_ids file
#' data(circle_ID)
#' circle_ids<-circle_ID
#'
#' #creating a temporary output_folder
#' output_folder<-tempdir()
#'
#' #Extracting exon for all circRNAs
#' Extract_exon(bed_file, folder_name, circle_ids, output_folder)

Extract_exon<-function(bed_file, folder_name, circle_ids, output_folder){

  #############Start of function cigar_as_vector()#########
  cigar_as_vector<-function(cigar){
    A=NULL
    for(i in 1:nchar(cigar)){
      A[i]<-substring(cigar, i,i)
    }
    return(A)
  }
  #############End of function cigar_as_vector()###########

  #############Start of function cigar_as_letter()#########
  cigar_as_letter<-function(cigar_vector){
    a<-sub("^([[:alpha:]]*).*", "\\1", cigar_vector)
    a<-a[-which(a=="")]
    return(a)
  }
  #############End of function cigar_as_letter()###########

  #############Start of function end_from_cigar()##########
  end_from_cigar<-function(s, cigar){
    cigar<-as.character(cigar)
    cigar_vector<-cigar_as_vector(cigar)
    cigar_letter<-cigar_as_letter(cigar_vector)
    start=NULL
    end=NULL
    for (i in 1:length(cigar_letter)){
      p<-which(cigar_vector==cigar_letter[i])[1]
      M<-as.numeric(paste(cigar_vector[1:p-1], collapse=""))
      if(cigar_letter[i]=="M" | cigar_letter[i]=="N") {
        start[i]=s-1; end [i]=start[i]+M
      } else {
        start[i]=s-1; end [i]=s-1
      }
      cigar_vector<-cigar_vector[-c(1:p)]
      s=end[i]+1
    }
    end=end[length(end)]
    return(end)
  }
  #############End of function end_from_cigar()############

  #############Start of function number_unique_read()######
  number_unique_read<-function(reads){
    dup<-duplicated(reads)
    uni_read<-reads[which(dup==F)]
    number_uni_read<-length(uni_read)
    return(number_uni_read)
  }
  #############End of function number_unique_read()########

  #################Start of function exon_count_bam()##############
  exon_count_bam<-function(bam_df, bed_df, circle_ids, file_list){
    transcript=NULL
    for (i in 1:dim(bam_df)[1]){
      start=bam_df$pos[i]
      cigar=as.character(bam_df$cigar[i])
      end=end_from_cigar(start, cigar)
      #b<-matrix(unlist(strsplit(as.character(file_list), "_")),ncol=4,byrow=T)
      b<-unlist(strsplit(as.character(file_list), "_"))
      t_id<-paste(b[1], b[2], sep=":")
      t_id_f<-paste(t_id, b[3], sep="|")
      transcript_id<-circle_ids$V2[which(circle_ids$V1==t_id_f)]
      index<-which(bed_df[,1]==as.character(bam_df$rname[i]) & bed_df[,4]==as.character(transcript_id) & bed_df[,2]<=end & bed_df[,3]>=start)
      if (length(index)>0){
        trans=bed_df[index,]
        a<-matrix(unlist(strsplit(as.character(trans[,5]), "_")),ncol=4,byrow=T)
        exon_id<-a[,4]
        transcript_exon<-paste(transcript_id, exon_id, sep="_")
      }else{
        trans=data.frame(0,0, 0, 0, 0, 0, 0)
        colnames(trans)<-c("V1", "V2", "V3", "V4", "V5", "V6", "V7")
        transcript_id=0
        exon_id=0
        transcript_exon=0
      }
      trans_1<-data.frame(bam_df[i,1], bam_df[i,3], bam_df[i,4], bam_df[i,7], bam_df[i,5]-1,end,trans, transcript_id, exon_id, transcript_exon)
      h<- sapply(trans_1, is.factor)
      trans_1[h] <- lapply(trans_1[h], as.character)
      transcript=rbind(transcript, trans_1)
    }
    dup<-duplicated(transcript$transcript_exon)
    transcript_id_uni<-transcript$transcript_exon[which(dup==F)]
    transcript_id=NULL
    exon_id=NULL
    start=NULL
    end=NULL
    reads=NULL
    strand=NULL
    n_uni_reads=NULL
    for(i in 1:length(transcript_id_uni)){
      a<-which(transcript$transcript_exon==transcript_id_uni[i])
      transcript_id[i]=as.character(transcript$transcript_id[a[1]])
      exon_id[i]=as.character(transcript$exon_id[a[1]])
      start[i]=transcript$V2[a[1]]
      end[i]=transcript$V3[a[1]]
      reads[i]=paste(transcript[a,1],collapse=",")
      strand[i]=as.character(transcript$V6[a[1]])
      n_uni_reads[i]<-number_unique_read(transcript[a,1])
    }
    exon_count<-data.frame(transcript_id, exon_id, start, end, strand, reads, n_uni_reads)
    return(exon_count)
  }
  #################End of function exon_count_bam()################

  file_list <- list.files(folder_name, pattern="*.txt")
  exon_count=NULL
  for (i in 1:length(file_list)){
    bam<-read.table(file.path(folder_name, file_list[i]), header=T, fill=T)[,1:8]
    transcript_all<-exon_count_bam(bam, bed_file, circle_ids, file_list[i])
    transcript_all<-transcript_all[which(transcript_all$transcript_id!=0),]
    if(nrow(transcript_all)>0){
      select_transcript_data<-transcript_all
      circle_split<-unlist(strsplit(file_list[i], "_"))[1:3]
      circle_id<-paste(circle_split, collapse="_")
      transcript_id<-transcript_all$transcript_id
      exon_id<-select_transcript_data$exon_id
      chrom<-circle_split[1]
      start<-select_transcript_data$start
      end<-select_transcript_data$end
      strand<-select_transcript_data$strand
      unique_reads<-select_transcript_data$n_uni_reads
      data_frame<-data.frame(circle_id, transcript_id, exon_id, chrom, start, end, strand, unique_reads)
      data_frame_sort<-data_frame[order(as.numeric(as.character(data_frame$start))),]
      write.table(data_frame_sort, file.path(output_folder, "exon_count.txt"),sep = "\t", col.names=F, row.names = FALSE, quote= FALSE, append=T)
    }
  }
}



#' Detecting skip exon
#'
#' This function detects skip exon
#' @param circle_reads_folder Location of the individual alignments (obtained from function \code{\link[circRNAFull]{extract_reads_from_bam}})
#' @param exon_count_file A data frame containing exons for circRNAs (obtained from function \code{\link[circRNAFull]{Extract_exon}})
#' @param output_folder The name of the output directory
#'
#' @return A file containing the skip exons written in 'output_folder'
#' @export
#'
#' @examples
#' #Example of circle_reads_folder containing the individual alignments
#' circle_reads_folder<-system.file('extdata', package = 'circRNAFull')
#'
#' #Loading an example exon_count_file containing exons for circRNAs
#' data(exon_count)
#' exon_count_file<-exon_count
#'
#' #Creating a temporary output directory
#' output_folder<-tempdir()
#'
#' #Detecting skip exon
#' skip_exon_detection(circle_reads_folder, exon_count_file, output_folder)

skip_exon_detection<-function(circle_reads_folder, exon_count_file, output_folder){

  #################Start of function cigar_as_vector()###########
  cigar_as_vector<-function(cigar){
    A=NULL
    for(i in 1:nchar(cigar)){
      A[i]<-substring(cigar, i,i)
    }
    return(A)
  }
  #################End of function cigar_as_vector()###########


  #################Start of function cigar_as_letter()###########
  cigar_as_letter<-function(cigar_vector){
    a<-sub("^([[:alpha:]]*).*", "\\1", cigar_vector)
    a<-a[-which(a=="")]
    return(a)
  }
  #################End of function cigar_as_letter()###########


  ########Start of function extract_intron_1()##########
  extract_intron_1<-function(bam){
    intron_data=NULL
    for(read in 1:dim(bam)[1]){
      cigar<-as.character(bam$cigar[read])
      cigar_vector<-cigar_as_vector(cigar)
      index_N<-which(cigar_vector=="N")
      if (length(index_N)>0){
        cigar_letter<-cigar_as_letter(cigar_vector)
        start=NULL
        end=NULL
        s=bam$pos[read]
        for (i in 1:length(cigar_letter)){
          p<-which(cigar_vector==cigar_letter[i])[1]
          M<-as.numeric(paste(cigar_vector[1:p-1], collapse=""))
          if(cigar_letter[i]=="M" | cigar_letter[i]=="N") {
            start[i]=s; end [i]=start[i]+M-1
            s=end[i]+1
          } else {
            start[i]=s; end [i]=s
          }
          cigar_vector<-cigar_vector[-c(1:p)]
        }
        all<-data.frame(start, end)
        u<-which(cigar_letter=="N")
        intron_start<-start[u];intron_end<-end[u]
      } else{
        intron_start<-0;intron_end<-0
      }
      intron_s_e<-data.frame(bam[read,1], bam[read,3], bam[read,7],intron_start, intron_end)
      intron_data<-rbind(intron_data, intron_s_e)
    }
    colnames(intron_data)<-c("qname", "chrom", "mapq", "intron_starts", "intron_ends")
    return(intron_data)
  }
  ######End of function extract_intron_1()#############

  ########Start of function unique_read()##############
  unique_read<-function(reads){
    dup<-duplicated(reads)
    uni_reads<-reads[which(dup==F)]
    return(uni_reads)
  }
  ########End of function unique_read()################

  exon_count<-exon_count_file
  file_list <- list.files(circle_reads_folder, pattern="*.txt")
  detect_skip_exon=data.frame(circle_id=character(), transcript=character(), skip_exon_start=numeric(), skip_exon_end=numeric(), intron_start=numeric(), intron_end=numeric(), unique_reads=character())
  write.table(detect_skip_exon, file.path(output_folder,"detect_skip_exon.txt"),sep="\t", row.names=F, quote=F)
  for (i in 1:length(file_list)){
    bam<-read.table(file.path(circle_reads_folder, file_list[i]), header=T, fill=T)
    intron_data<-extract_intron_1(bam)
    circle_split<-unlist(strsplit(file_list[i], "_"))[1:3]
    circle_id<-paste(circle_split, collapse="_")
    index_exon_count<-which(exon_count$V1==circle_id)
    exon_count_circle_id<-exon_count[index_exon_count,]
    ind_0<-which(intron_data$intron_end!=0)
    if(length(ind_0)>0){
      intron<-intron_data[ind_0,]
      skip_exon=data.frame(circ_id=character(0), trans_id=character(0), skip_exon_start=numeric(0), skip_exon_end=numeric(0), intron_start=numeric(0), intron_end=numeric(0), read=character(0))
      for(j in 1:nrow(intron)){
        index<-which(exon_count_circle_id$V5<=intron$intron_end[j]-1 & exon_count_circle_id$V6>=intron$intron_start[j])
        if (length(index)>0){
          skip_exon_start=NULL
          skip_exon_end=NULL
          for(t in 1:length(index)){
            skip_exon_start[t]<-max(exon_count_circle_id$V5[index[t]],intron$intron_start[j])
            skip_exon_end[t]<-min(exon_count_circle_id$V6[index[t]],intron$intron_end[j])
          }
          circ_id<-circle_id
          trans_id<-exon_count_circle_id$V2[1]
          intron_start<-intron$intron_start[j]
          intron_end<-intron$intron_end[j]
          read<-intron$qname[j]
          skip_exon_temp<-data.frame(circ_id, trans_id, skip_exon_start, skip_exon_end, intron_start, intron_end, read)
        }else{
          skip_exon_temp<-data.frame(circ_id=character(0), trans_id=character(0), skip_exon_start=numeric(0), skip_exon_end=numeric(0), intron_start=numeric(0), intron_end=numeric(0), read=character(0))
        }
        skip_exon<-rbind(skip_exon,skip_exon_temp)
      }
      if(nrow(skip_exon)>0){
        skip_exons<-paste(skip_exon[,3],skip_exon[,4], sep="_")
        dup<-duplicated(skip_exons)
        skip_exon_uni<-skip_exons[which(dup==F)]
        if (length(skip_exon_uni)==1) {
          skip_exon_1=data.frame(skip_exon[1,c(1:6)],unique_reads=paste(unique_read(skip_exon[,7]),collapse=","))
        }else{
          skip_exon_start=NULL
          skip_exon_end=NULL
          intron_start=NULL
          intron_end=NULL
          reads=NULL
          for( k in 1:length(skip_exon_uni)){
            a<-which(skip_exons==skip_exon_uni[k])
            skip_exon_start[k]=skip_exon[a[1],3]
            skip_exon_end[k]=skip_exon[a[1],4]
            intron_start[k]=skip_exon[a[1],5]
            intron_end[k]=skip_exon[a[1],6]
            reads[k]=paste(unique_read(skip_exon[a,7]), collapse=",")
          }
          skip_exon_1=data.frame(skip_exon[1,1],skip_exon[1,2],skip_exon_start, skip_exon_end,intron_start, intron_end, unique_reads=reads)
        }
      }else skip_exon_1=data.frame(circle_id=character(0), transcript=character(0), skip_exon_start=numeric(0), skip_exon_end=numeric(0), intron_start=numeric(0), intron_end=numeric(0), unique_reads=character(0))
    }else skip_exon_1=data.frame(circle_id=character(0), transcript=character(0), skip_exon_start=numeric(0), skip_exon_end=numeric(0), intron_start=numeric(0), intron_end=numeric(0), unique_reads=character(0))
    skip_exon_1_sort<-skip_exon_1[order(skip_exon_1[,3], skip_exon_1[,4]),]
    write.table(skip_exon_1_sort, file.path(output_folder,"detect_skip_exon.txt"),sep="\t", row.names=F, col.names=F, quote=F, append=T)
  }
}


#' Extracting exon after deleting skip exon
#'
#' This function extract the remaining exons after deleting skip exons
#' @param exon_count_file A data frame containing exons for circRNAs (obtained from function \code{\link[circRNAFull]{Extract_exon}})
#' @param skip_exon_file A data frame containing skip exons for circRNAs (obtained from function \code{\link[circRNAFull]{skip_exon_detection}})
#' @param output_folder The name of the output directory
#'
#' @return A text file containing the exons after deleting skip exon
#' @export
#'
#' @examples
#' #loading an example exon count file
#' data(exon_count)
#' exon_count_file<-exon_count
#'
#' #Loading an example skip exon file
#' data(skip_exon)
#' skip_exon_file<-skip_exon
#'
#' #Creating a temporary output directory
#' output_folder<-tempdir()
#'
#' #Extracting exon after deleting skip exon
#' extract_exon_after_delete_skip_exon(exon_count_file, skip_exon_file, output_folder)

extract_exon_after_delete_skip_exon<-function(exon_count_file, skip_exon_file, output_folder){
  exon_count<-exon_count_file
  skip_exon<-skip_exon_file
  for (i in 1:nrow(skip_exon)){
    ind<-which(exon_count$V1==as.character(skip_exon$circle_id[i]))
    circ_exon<-exon_count[ind,]
    index<-which(circ_exon$V5<=skip_exon$skip_exon_end[i] & circ_exon$V6>=skip_exon$skip_exon_start[i])
    if (length(index)>0){
      circ_exon_intersect_skip_exon<-circ_exon[index,]
      s=circ_exon_intersect_skip_exon$V5
      e=circ_exon_intersect_skip_exon$V6
      if(s==skip_exon$skip_exon_start[i] & e==skip_exon$skip_exon_end[i]){
        circ_exon_new=circ_exon[-index,]
      }else if (s>=skip_exon$skip_exon_start[i]){
        s_new<-skip_exon$skip_exon_end[i]+1
        circ_exon_new<-circ_exon
        circ_exon_new$V5[index]<-s_new
      }else if (e<=skip_exon$skip_exon_end[i]){
        e_new<-skip_exon$skip_exon_start[i]-1
        circ_exon_new<-circ_exon
        circ_exon_new$V6[index]<-e_new
      }else{
        s_new_1<-s
        e_new_1<-skip_exon$skip_exon_start[i]-1
        s_new_2<-skip_exon$skip_exon_end[i]+1
        e_new_2<-e
        s_e_new<-data.frame(V1=circ_exon[index,1], V2=circ_exon[index,2], V3=circ_exon[index,3], V4=circ_exon[index,4], V5=c(s_new_1,s_new_2), V6=c(e_new_1,e_new_2), V7=circ_exon[index,7], V8=circ_exon[index,8])
        circ_exon_new_1<-circ_exon[-index,]
        circ_exon_new<-rbind(circ_exon_new_1,s_e_new)
        circ_exon_new<-circ_exon_new[order(circ_exon_new$V3),]
      }
      exon_count<-exon_count[-ind,]
      exon_count<-rbind(exon_count,circ_exon_new)
    }
  }
  write.table(exon_count, file.path(output_folder,"extract_exon_after_skipping_exon.txt"),sep="\t", row.names=F, quote=F)
}





#' Reconstruction of full length circRNA sequences
#'
#' This function can reconstruct the full length circRNA
#' sequences from the output of the circular RNA predictions tools
#'
#' @param seq_name A vector containing the name of chromosomes of the reference genome
#' @param sequence A vector containing the sequences of the reference genome
#' @param exon_count_file A data frame containing exons after deleting skip exons
#' @param output The output folder name
#'
#' @return The reconstructed circRNA sequences will be written in folder 'output'
#'
#' @importFrom seqinr write.fasta
#' @importFrom stringi stri_reverse
#' @importFrom Biostrings readDNAStringSet
#'
#' @export
#'
#' @examples
#' \donttest{
#' \dontrun{
#' #Please download the reference genome. suppose you have downloaded
#' #the reference genome 'hg38.fa' in you working directory.
#' #Then run:
#' ref_genome<-"hg38.fa"
#' fastaFile <- readDNAStringSet(ref_genome)
#' seq_name = sub('\\ .*', '', names(fastaFile))
#' sequence = paste(fastaFile)
#' }
#' }
#'
#' #Loading an exon count file after deleting skip exon
#' exon_count<-data(exon_count_after_skipping_exon)
#' exon_count<-exon_count_after_skipping_exon
#'
#' #Creating an output directory
#' output<-tempdir()
#'
#' \donttest{
#' \dontrun{
#' #Reconstructing the circRNA sequences. The circRNA sequences will be generated
#' #in the 'output' directory.
#' extract_sequence(seq_name, sequence, exon_count, output)
#' }
#' }

extract_sequence<-function(seq_name, sequence, exon_count_file, output){
  dup<-duplicated(exon_count_file$V1)
  circle_ids_uni<-exon_count_file$V1[which(dup==F)]
  for (i in 1:length(circle_ids_uni)){
    circle_split<-unlist(strsplit(as.character(circle_ids_uni[i]), "_"))[1:3]
    chr<-circle_split[1]
    start<-as.numeric(circle_split[2])
    end<-as.numeric(circle_split[3])
    v<-which(chr==as.character(seq_name))
    B=substr(sequence[v], start, end)
    exon_start_end<-exon_count_file[which(exon_count_file$V1==circle_ids_uni[i]),c(5:7)]
    if(exon_start_end$V7[1]=="-")B<-chartr("acgtACGT", "tgcaTGCA", B)
    n_exon<-nrow(exon_start_end)
    e_sizes<-exon_start_end$V6-exon_start_end$V5+1
    e_offsets<-exon_start_end$V5-exon_start_end$V5[1]
    C=NULL
    for (k in 1:n_exon){
      C[k]<-substr(B, e_offsets[k]+1, e_offsets[k]+e_sizes[k])
    }
    D<-paste(C,collapse="")
    if(exon_start_end$V7[1]=="-") D<-stri_reverse(D)
    write.fasta(D, circle_ids_uni[i], file.path(output,"circRNAseq.fasta"), open = "a", nbchar = 100, as.string = FALSE)
  }
}
