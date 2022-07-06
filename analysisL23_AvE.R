library(eulerr)

allen_rna <- readRDS("/home/luoqinhuan/public-atac/adult-brain/allen_brain.rds")
allen_rna <- FindVariableFeatures(
  object = allen_rna,
  nfeatures = 5000
)

L23.markers <- FindMarkers(allen_rna, only.pos = TRUE,ident.1 = "L2/3 IT",min.pct = 0.2 ,group.by = "subclass")

peak <- rownames(`peak_cluster_filtered_L2/3 IT`)
motif_result = Find_feature_motif(peak,brain)
motif_result = motif_result[motif_result$pvalue<0.01,]

fragment <- Fragment_regularization(peak)
colnames(fragment) =  c("seq_name","start","end")
index_L1 <- Find_index_class(mm10_LINE_filtered,fragment)

peak_L1 <- peak[index_L1[[1]]]
peak_noL1 <- peak[which(!peak%in%peak_L1)]

## accessibility of these peak

accessibility_L1 = numeric()
for(i in 1:length(peak_L1)){
  accessibility_L1[i] = peak_adult_brain_quantile[peak_L1[i],2]
}
accessibility_noL1 = numeric()
for(i in 1:length(peak_noL1)){
  accessibility_noL1[i] = peak_adult_brain_quantile[peak_noL1[i],2]
}

accessibility = data.frame(access = c(accessibility_L1,accessibility_noL1), 
                           group = c(rep("LINE",length(accessibility_L1)),rep("no-LINE",length(accessibility_noL1))))
  
p1_3_1 =ggplot(accessibility,aes(x=factor(group),y=access,color=group))+
    geom_violin(aes(fill=group))+
    geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
    geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
    scale_fill_brewer(palette = "Pastel1")+
    #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
    labs(x="", y="Accessibility of Peak",title="LINE vs non-LINE")+theme_bw()+
    stat_compare_means(label.x = 1.35, label.y = 0.65,method = "t.test")+
    stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.62,method = "t.test")+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

gene_LINE = Find_Genes_with_Peak(peak_L1,brain_co,0.48)$id
gene_noLINE = Find_Genes_with_Peak(peak_noL1,brain_co,0.48)$id



vd <- Venn_two(
  lst.A <- gene_LINE,
  lst.B <- gene_noLINE
)

plot(vd,
     labels = list(
       labels = c('LINE', 'non-LINE'),
       col = "gray20", font = 2
     ), 
     edges = list(col="gray60", lex=1),
     fills = list(fill = c("#297CA0", "#BBBBBB"), alpha = 0.6),
     quantities = list(cex=.8, col='gray20')
)

## 比对A,B,AB三组基因list中，哪一组和cell marker中关系最为密切，以及overlapped gene 中谁的co-accessibility 更大
gene_L23 <- rownames(L23.markers)

genes_LINE<-mapIds(x=org.Mm.eg.db,
                        keys = gene_LINE,
                        keytype = "ENSEMBL",
                        column = "SYMBOL")
genes_LINE<-na.omit(genes_LINE)


genes_noLINE<-mapIds(x=org.Mm.eg.db,
                    keys = gene_noLINE,
                    keytype = "ENSEMBL",
                    column = "SYMBOL")
genes_noLINE<-na.omit(genes_noLINE)

gene_common_marker = gene_L23[which(gene_L23%in%genes_LINE)]
gene_common_marker = gene_common_marker[gene_common_marker%in%genes_noLINE]
gene_common_marker_EN =names(genes_noLINE[genes_noLINE%in%gene_common_marker])
index_common_marker = which(gene_common%in%gene_common_marker_EN)
  

vd <- Venn_three(
  lst.A <- genes_LINE,
  lst.B <- genes_noLINE,
  lst.C <- gene_L23
)

enrich.go.BP<-enrichGO(gene = gene_LINE,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       readable = T)
dotplot(enrich.go.BP)


p3_3=plot(vd,
     labels = list(
       labels = c('LINE', 'non-LINE', 'Marker'),
       col = "gray20", font = 2
     ), 
     edges = list(col="gray60", lex=1),
     fills = list(fill = c("#297CA0", "#BBBBBB", "#E9EA77"), alpha = 0.6),
     quantities = list(cex=.8, col='gray20')
)

links = brain_co@assays$peaks@links
names(links) = 1:length(links)
links_start = GRanges(seqnames = links@seqnames,
                             ranges = IRanges(start(links),end= start(links),names = 1:length(links)),
                             strand = links@strand)
links_end = GRanges(seqnames = links@seqnames,
                           ranges = IRanges(end(links),end= end(links),names = 1:length(links)),
                           strand = links@strand)
links1_index = which(start(links)%in%start(links_start[findOverlaps(links_start,genes_granges,type = "within")@from]))
links2_index = which(end(links)%in%end(links_end[findOverlaps(links_end,genes_granges,type = "within")@from]))
links_genes = unique(links[unique(c(links1_index,links2_index))])


gene_common = gene_LINE[which(gene_LINE%in%gene_noLINE)]

score_LINE <- numeric()
score_noLINE <- numeric()

for(i in 1:length(gene_common)){
  score_LINE[i] = Coaccessibility_Evaluation(gene_common[i],links_genes,peak_L1,genes_granges)
  score_noLINE[i] = Coaccessibility_Evaluation(gene_common[i],links_genes,peak_noL1,genes_granges)
  print(i)
}

gene_common_data = data.frame(score=c(score_LINE,score_noLINE),group=as.factor(c(rep("LINE",length(score_LINE)),rep("non-LINE",length(score_LINE)))),
                              set = c(1:length(score_LINE),1:length(score_noLINE)))
p1_3_2= ggplot(data=gene_common_data,aes(x=group,y=score,color=group))+
      geom_violin(aes(fill=group))+
      geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
      geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
      scale_fill_brewer(palette = "Pastel1")+
      #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
      labs(x="LINE vs non-LINE", y="Coaccessibility_Score")+theme_bw()+
      stat_compare_means(label.x = 1.35, label.y = 0.7,method = "t.test",paired = TRUE)+
    stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.65,method = "t.test",paired = TRUE)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())



score_marker = data.frame(score=c(score_LINE[index_common_marker],score_noLINE[index_common_marker]),group=as.factor(c(rep("LINE",length(index_common_marker)),rep("non-LINE",length(index_common_marker)))),
                           set = c(1:length(index_common_marker),1:length(index_common_marker)))

p1_3_4=ggplot(data=score_marker,aes(x=group,y=score,color=group))+
  geom_violin(aes(fill=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_fill_brewer(palette = "Pastel1")+
  #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
  labs(x="LINE vs non-LINE marker gene", y="Coaccessibility Score")+theme_bw()+
  stat_compare_means(label.x = 1.35, label.y = 0.62,method = "t.test",paired = TRUE)+
  stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.59,method = "t.test",paired = TRUE)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
## Functional Analysis
gene_LINE_unique = gene_LINE[which(!gene_LINE%in%gene_noLINE)]
gene_LINE_unique_sym <-mapIds(x=org.Mm.eg.db,
                        keys = gene_LINE_unique,
                        keytype = "ENSEMBL",
                        column = "SYMBOL")
gene_LINE_unique_sym <-na.omit(gene_LINE_unique_sym)

marker_LINE_unique = gene_LINE_unique_sym[which(gene_LINE_unique_sym%in%gene_L23)]

enrich.go.BP<-enrichGO(gene = gene_LINE,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       readable = T)



fragment <- Fragment_regularization(peak_L1)
colnames(fragment) =  c("seq_name","start","end")
index<-Find_index_class(mm10_LINE_filtered,fragment)
fragment_new = Fragment_Partition(fragment,index)
fragment_new[is.na(fragment_new)]=0 

## Repeat the fragment_new two times to produce +/-chain
fragment.expanded <- fragment_new[rep(row.names(fragment_new), each=2), 1:4]
colnames(fragment.expanded) = c("seq","start","end","status")
fragment.expanded$strand = rep(c("+","-"),dim(fragment.expanded)[1]/2)
fragment.expanded.LINE = fragment.expanded[fragment.expanded$status==1,]
fragment.expanded.background = fragment.expanded[fragment.expanded$status==0,]
fragment.expanded.background = fragment.expanded.background[,-4]
fragment.expanded.LINE = fragment.expanded.LINE[,-4]

rownames(fragment.expanded.LINE) = c(1:dim(fragment.expanded.LINE)[1])
rownames(fragment.expanded.background) = c((dim(fragment.expanded.LINE)[1]+1):(dim(fragment.expanded.background)[1]+dim(fragment.expanded.LINE)[1]))

## tab to seperate file
write.table(fragment.expanded.LINE,"/home/luoqinhuan/LS_project/L2_neurons/peak_LINE.txt",col.names = FALSE, quote = FALSE,sep = "\t")
write.table(fragment.expanded.background,"/home/luoqinhuan/LS_project/L2_neurons/peak_Background.txt",col.names = FALSE, quote = FALSE,sep = "\t")


neurons = c("L2")
file_path1 = "/home/luoqinhuan/LS_project/"
file_path2 = "L2_neurons/"
pattern_interest = "POU"
correspondance_table = data.frame(matrix(0,2,2))
correspondance_table[1,]=c("OCT2","POU2F2")
correspondance_table[2,]=c("OCT6","POU3F1")
correspondance_table[3,]=c("BRN1","POU3F3")
correspondance_table[4,]=c("BRN2","POU3F2")

path = paste(file_path1,file_path2,sep = "",collapse = "")
result1 = summaryHomer(path)[,c(1,2,4,5,6)]
colnames(result1)[2] = c("TF")
result2 = summaryHomerKnown(path)
result_homer = rbind(result1,result2)
index_keep = numeric()
for (i in 1:dim(result_homer)[1]) {
  if(result_homer$TF[i]%in%correspondance_table$X1){
    result_homer$TF[i] =correspondance_table[which(correspondance_table$X1==result_homer$TF[i]),2]
  }
  if(result_homer$TF[i]%in%motif_L2){
    index_keep[i] = 1
  }
  else index_keep[i] = 0
}
result_homer = result_homer[index_keep==1,]

TF_name = result_homer[grep(pattern = pattern_interest,result_homer[,2]),2]
motif = rownames(result_homer[grep(pattern = pattern_interest,result_homer[,2]),])
homer_motif = motif[grep(pattern = "motif",motif)]
known_motif = motif[grep(pattern = "known",motif)]
motif = c(homer_motif,known_motif)
if(length(homer_motif)>0) cmd1 = homer_motif_analysis(homer_motif,path) else cmd1=character()
if(length(known_motif)>0) cmd2 = known_motif_analysis(known_motif,path) else cmd2 = known_motif_analysis(known_motif,path)
cmd = c(cmd1,cmd2)
write.table(cmd,file = paste(path,"command_homer.sh",sep="",collapse = ""),sep = "/n",quote = FALSE,row.names = FALSE,col.names = FALSE)

links = brain_co@assays$peaks@links
###After execute the command, then analyze the result of command
gene_list = list()
fragment_list = list()
for(i in 1:length(motif)){
  print(i)
  file_path = paste(path,motif[i],"_peak.txt",sep="",collapse = "")
  result <- read.delim(file_path)
  id <- result$PositionID
  peak_motif <- fragment.expanded.LINE[id,]
  motif_peak <- character()
  for(j in 1:dim(peak_motif)[1]){
    motif_peak[j] = paste0(peak_motif[j,1:3],collapse = "-",sep="")
  }
  motif_peak <- unique(motif_peak)
  
  ## get fragment contain these motif_peak
  peak_motif <- Fragment_regularization(motif_peak)
  colnames(peak_motif) = c("seq_name","start","end")
  index = Find_index(fragment,peak_motif)
  fragment_motif <- index[[3]]
  fragment_motif <- fragment_motif[!duplicated(fragment_motif),]
  colnames(fragment_motif) = c("seq_name","start","end")
  threshold_coaccessibility = 0.48
  links_subset <- links[which(links$score>threshold_coaccessibility)]
  subset_links <- as.data.frame(links_subset)
  subset_links <- subset_links[,1:3]
  
  ## Filter out links that involved with peak cared about
  test = Find(subset_links,fragment_motif)
  peak_contain_motif_link = test[[2]]
  links_subset_peak = links_subset[test[[1]],]
  
  #### Find gene involved with this potential motif
  links_peak_motif <- as.data.frame(links_subset_peak)
  links_peak_motif = links_peak_motif[,1:3]
  colnames(links_peak_motif) = c("seq_name","start","end")
  test = Find(links_peak_motif,genes_range)
  genes = unique(all_genes[which(gene_location%in%test[[2]]),]$gene_name)
  gene_list[[i]] = genes
  names(gene_list)[i] = TF_name[i]
  test_total = Find(links_peak_motif,fragment_total)
  test_motif = Find(links_peak_motif,fragment_motif)
  peak_total_link = test_total[[2]]
  peak_motif_link = test_motif[[2]]
  peak_rm_motif  = peak_total_link[!peak_total_link%in%peak_motif_link]
  fragment_motif_link = Peak_Expansion(peak_rm_motif)
  fragment_list[[i]] = fragment_motif_link
  names(fragment_list)[i] = TF_name[i]
}


for(i in 1:length(gene_list)){
  if(names(gene_list[i])%in%correspondance_table$X1){
    names(gene_list)[i] = correspondance_table[which(correspondance_table$X1==names(gene_list[i])),2]
  }
}
name_gene = unique(names(gene_list))
gene_list_merge = list()
for(i in 1:length(name_gene)){
  gene_list_merge[i]=list(as.character(unlist((gene_list[which(names(gene_list)==name_gene[i])]))))
  names(gene_list_merge)[i] = name_gene[i]
}
variable_name = paste("gene_list","_",pattern_interest,sep="",collapse = "")
assign(variable_name,gene_list_merge)
saveRDS(get(variable_name),paste(path,variable_name,".rds",sep="",collapse = ""))
for(i in 1:length(gene_list_merge)){
  if(length(gene_list_merge[[i]])!=0){
    genes = gene_list_merge[[i]]
    genes.entrez_id<-mapIds(x=org.Mm.eg.db,
                            keys = genes,
                            keytype = "SYMBOL",
                            column = "ENTREZID")
    genes.entrez_id<-na.omit(genes.entrez_id)
    
    erich.go.BP<-enrichGO(gene = genes.entrez_id,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          readable = T)
    assign(paste("plot",i,sep="",collapse = ""),dotplot(erich.go.BP,title = paste(name_gene[i],"BP",sep = "_",collapse = "")))
    if(dim(get(paste("plot",i,sep="",collapse = ""))[[1]])[1]>0){
      plot_path = paste(path,paste(pattern_interest,"_family_","plot",i,".png",sep = "",collapse = ""),sep = "",collapse = "")
      png(file=plot_path, width = 1000,height = 800)
      print(dotplot(erich.go.BP,title = paste(name_gene[i],"BP",sep = "_",collapse = "")))
      dev.off()
    }
  }
}

gene = character()
for(i in 1:length(gene_list)){
  gene = c(gene,gene_list[[i]])
}
gene = unique(gene)
paste(gene,sep = "",collapse = " ")
## Protein Protein interaction Analysis: Finding long-distance motif regulatory relationships
# For each transcription factor interested in
## Firstly, prepare homer prediction command
name_fragment = unique(names(fragment_list))
fragment_list_merge = list()
for(i in 1:length(name_gene)){
  int = data.frame(matrix(0,0,4))
  fragment_list_sub = fragment_list[which(names(fragment_list)==name_fragment[i])]
  for(j in 1:length(fragment_list_sub)){
    int = rbind(int,fragment_list_sub[[j]])
  }
  fragment_list_merge[[i]]= int
  names(fragment_list_merge)[i] = name_gene[i]
}

fragment_list = fragment_list_merge
pattern_interest = "POU"
cmd_fragment = character()
num = 1
for(i in 1:length(fragment_list)){
  fragment_file = paste(path,"fragment",i,"_",name_fragment[i],"_",pattern_interest,".txt",sep = "",collapse = "")
  write.table(fragment_list[[i]],fragment_file,col.names = FALSE, quote = FALSE,sep = "\t")
  output_path = paste(path,"fragment",i,"_",name_fragment[i],"_",pattern_interest,sep="",collapse = "")
  system(paste("mkdir ", output_path,sep = "",collapse = ""))
  cmd_fragment[i] = paste("findMotifsGenome.pl ",fragment_file," mm10 ",output_path," -bg /home/luoqinhuan/LS_project/L2_neurons/peak_Background_genebody.txt",sep = "",collapse = "")
  
}
write.table(cmd_fragment,file = paste(path,"command_homer","_",pattern_interest,".sh",sep="",collapse = ""),sep = "/n",quote = FALSE,row.names = FALSE,col.names = FALSE)

## After exec
result_interaction = list()
for(i in 1:length(fragment_list)){
  output_path = paste(path,"fragment",i,"_",name_fragment[i],"_",pattern_interest,sep="",collapse = "")
  related_tf = Homer_Result_Disposal(output_path)[,2]
  for(j in 1:length(related_tf)){
    if(related_tf[j]%in%correspondance_table$X1){
      related_tf[j] = correspondance_table[which(correspondance_table$X1==related_tf[j]),2]
    }
  }
  tf = name_fragment[i]
  tf = unlist(lapply(strsplit(tf, split = "[^[:alnum:]]"),function(x) x[[1]][1]))
  tf = capitalize(tolower(tf))
  if(tf=="Brn1"){
    tf = "Pou3f3"
  }
  interaction_protein = Find_protein(tf)
  label = rep(0,length(interaction_protein[,1]))
  for(j in 1:length(label)){
    if(length(grep(pattern=interaction_protein[j,1],related_tf))){
      label[j] = 1
    }
  }
  result_interaction[[i]]= interaction_protein[which(label==1),1:2]
  names(result_interaction)[i]=tf
}

variable_name = paste("result_interaction","_",pattern_interest,sep="",collapse = "")
assign(variable_name,result_interaction)
saveRDS(get(variable_name),paste(path,variable_name,".rds",sep="",collapse = ""))



## Create Background fragment in gene body
peak_total_adult_gr = GRanges(seqnames = fragment_total$seq_name,
                  ranges = IRanges(start = as.numeric(fragment_total$start),end=as.numeric(fragment_total$end),names = c(1:dim(fragment_total)[1])),
                  strand = rep("*",dim(fragment_total)[1]))

fragment_in_gene <- peak_total_adult_gr[findOverlaps(peak_total_adult_gr,genes_granges,type = "within")@from]
fragment_in_gene_adult <- data.frame(seq_name=fragment_in_gene@seqnames,start = start(fragment_in_gene),end = end(fragment_in_gene))

fragment_interest <- rbind(fragment_list[[1]],fragment_list[[2]],fragment_list[[3]],fragment_list[[4]])

fragment_in_gene_adult%in%fragment_interest[,1:3]
fragment_interest = fragment_interest[,1:3]
fragment_interest$start = as.numeric(fragment_interest$start)
fragment_interest$end = as.numeric(fragment_interest$end)
fragment_backgroud = setdiff(fragment_in_gene_adult,fragment_interest)


fragment_backgroud.expanded <- fragment_backgroud[rep(row.names(fragment_backgroud), each=2), 1:3]
fragment_backgroud.expanded$strand = rep(c("+","-"),dim(fragment_backgroud.expanded)[1]/2)
rownames(fragment_backgroud.expanded) = c(1:dim(fragment_backgroud.expanded)[1])
write.table(fragment.expanded.background,"/home/luoqinhuan/LS_project/L2_neurons/peak_Background_genebody.txt",col.names = FALSE, quote = FALSE,sep = "\t")




pattern_interest = "POU"
cmd_fragment = character()
num = 1
for(i in 1:length(fragment_list)){
  fragment_file = paste(path,"fragment",i,"_",motif[i],"_",pattern_interest,".txt",sep = "",collapse = "")
  write.table(fragment_list[[i]],fragment_file,col.names = FALSE, quote = FALSE,sep = "\t")
  output_path = paste(path,"fragment",i,"_",motif[i],"_",pattern_interest,sep="",collapse = "")
  system(paste("mkdir ", output_path,sep = "",collapse = ""))
  cmd_fragment[i] = paste("findMotifsGenome.pl ",fragment_file," mm10 ",output_path," -bg /home/luoqinhuan/LS_project/L2_neurons/peak_Background_genebody.txt",sep = "",collapse = "")
  
}
write.table(cmd_fragment,file = paste(path,"command_homer","_",pattern_interest,".sh",sep="",collapse = ""),sep = "/n",quote = FALSE,row.names = FALSE,col.names = FALSE)
result_interaction = list()
for(i in 1:length(fragment_list)){
  output_path = paste(path,"fragment",i,"_",motif[i],"_",pattern_interest,sep="",collapse = "")
  related_tf = Homer_Result_Disposal(output_path)[,2]
  for(j in 1:length(related_tf)){
    if(related_tf[j]%in%correspondance_table$X1){
      related_tf[j] = correspondance_table[which(correspondance_table$X1==related_tf[j]),2]
    }
  }
  tf = TF_name[i]
  tf = unlist(lapply(strsplit(tf, split = "[^[:alnum:]]"),function(x) x[[1]][1]))
  tf = capitalize(tolower(tf))
  if(tf=="Brn1"){
    tf = "Pou3f3"
  }
  interaction_protein = Find_protein(tf)
  label = rep(0,length(interaction_protein[,1]))
  for(j in 1:length(label)){
    if(length(grep(pattern=interaction_protein[j,1],related_tf))){
      label[j] = 1
    }
  }
  result_interaction[[i]]= interaction_protein[which(label==1),1:2]
  names(result_interaction)[i]=tf
}

variable_name = paste("result_interaction","_",pattern_interest,"_not_merge",sep="",collapse = "")
assign(variable_name,result_interaction)
saveRDS(get(variable_name),paste(path,variable_name,".rds",sep="",collapse = ""))

#### Comparative Analysis 

##Overlap peak 
Excitatory.markers <- FindMarkers(allen_rna, only.pos = TRUE,ident.1 = "Glutamatergic",min.pct = 0.2 ,group.by = "class")
E18_Excitatory_peak = FindMarkers(object = E18_brain,ident.1 = "A:Excitatory neurons" ,min.pct = 0.1,logfc.threshold = 0.1,
                                  test.use = 'LR',latent.vars = 'peak_region_fragments',only.pos = TRUE,group.by = "group")
Adult_Excitatory_peak = FindMarkers(object = brain,ident.1 = "A:Excitatory neurons" ,min.pct = 0.1,logfc.threshold = 0.2,
                                    test.use = 'LR',latent.vars = 'peak_region_fragments',only.pos = TRUE,group.by = "group")

peak_Excitatory_adult <- rownames(Adult_Excitatory_peak)
peak_Excitatory_adult_LINE <- Peak_with_LINE(peak_Excitatory_adult,mm10_LINE_filtered)

peak_Excitatory_E18 <- rownames(E18_Excitatory_peak)
peak_Excitatory_E18_LINE <- Peak_with_LINE(peak_Excitatory_E18,mm10_LINE_filtered)


peak_E18_Ex_L_gr <- Peak_to_GRanges(peak_Excitatory_E18_LINE)
peak_Adult_Ex_L_gr <- Peak_to_GRanges(peak_Excitatory_adult_LINE)
peak_correspondence_table = data.frame(findOverlaps(peak_E18_Ex_L_gr,peak_Adult_Ex_L_gr,type = "any")@from,findOverlaps(peak_E18_Ex_L_gr,peak_Adult_Ex_L_gr,type = "any")@to)
colnames(peak_correspondence_table) = c("E18","Adult")
peak_adult_brain_quantile <-FindTopFeatures(object = brain_co[['peaks']][])
peak_E18_brain_quantile <-FindTopFeatures(object = E18_link[['peaks']][])



accessibility_table = peak_correspondence_table
for(i in 1:dim(peak_correspondence_table)[1]){
  peak1 = peak_Excitatory_E18_LINE[peak_correspondence_table[i,1]]
  peak2 = peak_Excitatory_adult_LINE[peak_correspondence_table[i,2]]
  accessibility_table[i,1] = peak_E18_brain_quantile[peak1,2]
  accessibility_table[i,2] = peak_adult_brain_quantile[peak2,2]
  print(i)
}

accessibility_table_reshape = data.frame(accessibility_score=c(accessibility_table$E18,accessibility_table$Adult),
                                         group = c(rep("E18",length(accessibility_table$E18)),rep("Adult",length(accessibility_table$E18))))

score_table = peak_correspondence_table
for(i in 1:dim(peak_correspondence_table)[1]){
  peak1 = peak_Excitatory_E18_LINE[peak_correspondence_table[i,1]]
  peak2 = peak_Excitatory_adult_LINE[peak_correspondence_table[i,2]]
  score_table[i,1] = mean(Coaccessibility_Calculation_score(peak1,E18_link))
  score_table[i,2] = mean(Coaccessibility_Calculation_score(peak2,brain_co))
  print(i)
}


p4_1= ggplot(data=accessibility_table_reshape,aes(x=group,y=accessibility_score,color=group))+
    geom_violin(aes(fill=group))+
    geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
    geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
    scale_fill_brewer(palette = "Pastel1")+
    #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
    labs(x="Adult vs E18 overlapped LINE peak", y="Accessibility Score")+theme_bw()+
    stat_compare_means(label.x = 1.35, label.y = 0.7,method = "t.test",paired = TRUE)+
    stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.68,method = "t.test",paired = TRUE)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p4_2=ggplot(data=accessibility_table_reshape,aes(x=group,y=coaccessibility_score,color=group))+
  geom_violin(aes(fill=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_fill_brewer(palette = "Pastel1")+
  #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
  labs(x="Adult vs E18 overlapped LINE peak", y="Coaccessibility Score")+theme_bw()+
  stat_compare_means(label.x = 1.35, label.y = 0.44,method = "t.test",paired = TRUE)+
  stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.42,method = "t.test",paired = TRUE)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())


accessibility_table_reshape$set = rep(1:dim(peak_correspondence_table)[1],2)

p4_3=ggplot(data=accessibility_table_reshape)+
  geom_point(size=2.5,aes(x=accessibility_score,y=coaccessibility_score,color=group,shape=group))+
  scale_color_brewer(palette = "Pastel1")+
  geom_line(aes(x=accessibility_score,y=coaccessibility_score,group=set),linetype=2,size=0.1)

accessibility_table_reshape$coaccessibility_score <- c(score_table$E18,score_table$Adult)

ggplot(data=accessibility_table_reshape)+
  geom_point(size=2.5,aes(x=accessibility_score,y=coaccessibility_score,color=group))+
  scale_color_brewer(palette = "Pastel1")+
  geom_line(aes(x=accessibility_score,y=coaccessibility_score,group=set),linetype=2,size=0.1)

accessibility_table_reshape_new = accessibility_table_reshape[245:488,]
accessibility_table_reshape_new$accessibility_score = accessibility_table_reshape_new$accessibility_score -
  accessibility_table_reshape[1:244,]$accessibility_score
accessibility_table_reshape_new$coaccessibility_score = accessibility_table_reshape_new$coaccessibility_score -
  accessibility_table_reshape[1:244,]$coaccessibility_score

accessibility_table_reshape_new$group = "static"
accessibility_table_reshape_new[which(accessibility_table_reshape_new$accessibility_score>0.01),]$group = "Up-regulated"
accessibility_table_reshape_new[which(accessibility_table_reshape_new$accessibility_score<0.01*-1),]$group = "Down-regulated"
accessibility_table_reshape_new[which(accessibility_table_reshape_new$coaccessibility_score>0.01&accessibility_table_reshape_new$accessibility_score>0.01),]$group = "Both-Positive"
accessibility_table_reshape_new[which(accessibility_table_reshape_new$accessibility_score<0.01*-1&accessibility_table_reshape_new$coaccessibility_score<0.01*-1),]$group = "Both-Negative"


p=ggplot(data=accessibility_table_reshape_new)+
  geom_point(size=2.5,aes(x=accessibility_score,y=coaccessibility_score,color=group))+
  scale_color_brewer(palette = "Accent")+
  geom_hline(yintercept = c(0.01),color="red")+
  geom_vline(xintercept = c(0.01),color="red")+
  theme_bw()+
  theme(legend.position = 'None',panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2=ggMarginal(p,type="densigram",xparams = list(fill ="#d6eaff",color="lightgrey"),yparams = list(fill ="#a8d3ff",color= "skyblue"))
ggsave("~/p2_5.jpg",p2,width = 4,height = 4,dpi=300)



accessibility_table$difference = accessibility_table$Adult - accessibility_table$E18
upregulated_peak = peak_Excitatory_adult_LINE[peak_correspondence_table[which(accessibility_table$difference>0.01),2]]
gene1 = Find_Genes_with_Peak(upregulated_peak,brain_co,0.48)$id

enrich.go.BP<-enrichGO(gene = gene1,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       readable = T)

dotplot(enrich.go.BP,title="The Pathway linked with Up-Regulated Peak")
p4_4 = dotplot(enrich.go.BP,title="The Pathway linked with Up-Regulated Peak")




accessibility_score_total = c(peak_E18_brain_quantile[peak_Excitatory_E18_LINE,2],peak_adult_brain_quantile[peak_Excitatory_adult_LINE,2])
group = c(rep("E18",length(peak_Excitatory_E18_LINE)),rep("Adult",length(peak_Excitatory_adult_LINE)))
data = data.frame(accessibility_score=accessibility_score_total,group)

p4_5=ggplot(data=data,aes(x=group,y=accessibility_score,color=group))+
  geom_violin(aes(fill=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_fill_brewer(palette = "Pastel1")+
  #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
  labs(x="Adult vs E18 LINE peak", y="Accessibility Score")+theme_bw()+
  stat_compare_means(label.x = 1.35, label.y = 0.7,method = "t.test")+
  stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.68,method = "t.test")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())


score1 = as.numeric()
score2 = as.numeric()
for(i in 1:length(peak_Excitatory_E18_LINE)){
  score1[i] = mean(Coaccessibility_Calculation_score(peak_Excitatory_E18_LINE[i],E18_link))
  print(i)
}

for(i in 1:length(peak_Excitatory_adult_LINE)){
  score2[i] = mean(Coaccessibility_Calculation_score(peak_Excitatory_adult_LINE[i],brain_co))
  print(i)
}

data$coaccessibility_score=c(score1,score2)
p4_6= ggplot(data=data,aes(x=group,y=coaccessibility_score,color=group))+
  geom_violin(aes(fill=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_fill_brewer(palette = "Pastel1")+
  #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
  labs(x="Adult vs E18 LINE CREs", y="Coaccessibility Score")+theme_bw()+
  stat_compare_means(label.x = 1.35, label.y = 0.45,method = "t.test")+
  stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.42,method = "t.test")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p2_1 = (p4_5+p4_6)/(p4_1+p4_2)
ggsave("~/p2_1.jpg",p2_1,width = 12,height = 10,dpi=300)

ggplot(data=data)+
  geom_point(size=2.5,aes(x=accessibility_score,y=coaccessibility_score,color=group,shape=group))+
  scale_color_brewer(palette = "Pastel1")+
  labs(title = "E18")+theme_bw()




png(filename ="~/LS_project/Figure3_3.jpg",width = 1800,height = 2400,res = 300)
p3_3
dev.off()
