Adult_Excitatory_peak = FindMarkers(object = brain,ident.1 = "A:Excitatory neurons" ,min.pct = 0.1,logfc.threshold = 0.2,
                                    test.use = 'LR',latent.vars = 'peak_region_fragments',only.pos = TRUE,group.by = "group")

Adult_Inhibitory_peak = FindMarkers(object = brain,ident.1 = "B:Inhibitory Interneurons" ,min.pct = 0.05,logfc.threshold = 0.2,
                                    test.use = 'LR',latent.vars = 'peak_region_fragments',only.pos = TRUE,group.by = "group")

E18_Excitatory_peak = FindMarkers(object = E18_brain,ident.1 = "A:Excitatory neurons" ,min.pct = 0.1,logfc.threshold = 0.1,
                                  test.use = 'LR',latent.vars = 'peak_region_fragments',only.pos = TRUE,group.by = "group")
E18_Inhibitory_peak = FindMarkers(object = E18_brain,ident.1 = "B:Inhibitory Interneurons" ,min.pct = 0.05,logfc.threshold = 0.1,
                                  test.use = 'LR',latent.vars = 'peak_region_fragments',only.pos = TRUE,group.by = "group")

Excitatory.markers <- rownames(FindMarkers(allen_rna, only.pos = TRUE,ident.1 = "Glutamatergic",min.pct = 0.2 ,group.by = "class"))
Inhibitory.markers <- rownames(FindMarkers(allen_rna, only.pos = TRUE,ident.1 = "GABAergic",min.pct = 0.1 ,group.by = "class"))

peak_Excitatory_adult <- rownames(Adult_Excitatory_peak)
peak_Excitatory_adult_LINE <- Peak_with_LINE(peak_Excitatory_adult,mm10_LINE_filtered)

peak_Inhibitory_adult <- rownames(Adult_Inhibitory_peak)
peak_Inhibitory_adult_LINE <- Peak_with_LINE(peak_Inhibitory_adult,mm10_LINE_filtered)


peak_Excitatory_E18 <- rownames(E18_Excitatory_peak)
peak_Excitatory_E18_LINE <- Peak_with_LINE(peak_Excitatory_E18,mm10_LINE_filtered)

peak_Inhibitory_E18 <- rownames(E18_Inhibitory_peak)
peak_Inhibitory_E18_LINE <- Peak_with_LINE(peak_Inhibitory_E18,mm10_LINE_filtered)

score_inhibitory = numeric()
for(i in 1:length(peak_Inhibitory_adult_LINE)){
  score_inhibitory[i] = mean(Coaccessibility_Calculation_score(peak_Inhibitory_adult_LINE[i],brain_co))
  print(i)
}

score_inhibitory_E18 = numeric()
for(i in 349:length(peak_Inhibitory_E18_LINE)){
  score_inhibitory_E18[i] = mean(Coaccessibility_Calculation_score(peak_Inhibitory_E18_LINE[i],E18_link))
  print(i)
}
score_inhibitory_E18 = na.omit(score_inhibitory_E18)

accessibility_score_inhibitory = peak_adult_brain_quantile[peak_Inhibitory_adult_LINE,2]
accessibility_score_excitatory = peak_adult_brain_quantile[peak_Excitatory_adult_LINE,2]
data_InvEx <- data.frame(accessibility_score = c(accessibility_score_inhibitory,accessibility_score_excitatory),
                         coaccessibility_score = c(score_inhibitory,score2),
                         group =c(rep("Inhibitory",length(peak_Inhibitory_adult_LINE)),
                                  rep("Excitatory",length(peak_Excitatory_adult_LINE))))
accessibility_score_inhibitory_E18 = peak_E18_brain_quantile[peak_Inhibitory_E18_LINE,2]
accessibility_score_inhibitory_E18 =accessibility_score_inhibitory_E18[-95]
data_InvEx$Time = "Adult"
data_InvEx2 <- data.frame(accessibility_score = c(accessibility_score_inhibitory_E18,peak_E18_brain_quantile[peak_Excitatory_E18_LINE,2]),
                         coaccessibility_score = c(score_inhibitory_E18,score1),
                         group =c(rep("Inhibitory",(length(peak_Inhibitory_E18_LINE)-1)),
                                  rep("Excitatory",length(peak_Excitatory_E18_LINE))))
data_InvEx2$Time = "E18"

p1=ggplot(data=data_InvEx,aes(x=group,y=accessibility_score,color=group))+
  geom_violin(aes(fill=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_fill_brewer(palette = "Pastel1")+
  #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
  labs(x="Excitatory vs Inhibitory LINE peak (Adult)", y="Accessibility Score")+theme_bw()+
  stat_compare_means(label.x = 1.35, label.y = 0.2,method = "t.test")+
  stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.12,method = "t.test")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p4=ggplot(data=data_InvEx2,aes(x=group,y=accessibility_score,color=group))+
  geom_violin(aes(fill=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_fill_brewer(palette = "Pastel1")+
  #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
  labs(x="Excitatory vs Inhibitory LINE peak (E18)", y="Accessibility Score")+theme_bw()+
  stat_compare_means(label.x = 1.35, label.y = 0.45,method = "t.test")+
  stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.4,method = "t.test")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p2=ggplot(data=data_InvEx,aes(x=group,y=coaccessibility_score,color=group))+
  geom_violin(aes(fill=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_fill_brewer(palette = "Pastel1")+
  #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
  labs(x="Excitatory vs Inhibitory LINE peak (Adult)", y="Coaccessibility Score")+theme_bw()+
  stat_compare_means(label.x = 1.35, label.y = 0.45,method = "t.test")+
  stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.4,method = "t.test")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p5=ggplot(data=data_InvEx2,aes(x=group,y=coaccessibility_score,color=group))+
  geom_violin(aes(fill=group))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_fill_brewer(palette = "Pastel1")+
  #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
  labs(x="Excitatory vs Inhibitory LINE peak (E18)", y="Coaccessibility Score")+theme_bw()+
  stat_compare_means(label.x = 1.35, label.y = 0.42,method = "t.test")+
  stat_compare_means(label = "p.signif",label.x = 1.5, label.y = 0.38,method = "t.test")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

data = rbind(data_InvEx,data_InvEx2)
data_Ex = data[data$group=="Excitatory",]
data_In = data[data$group=="Inhibitory",]
p3=ggplot(data=data_In)+
  geom_point(size=2,aes(x=accessibility_score,y=coaccessibility_score,color=Time,shape=Time),alpha=0.7)+
  scale_color_brewer(palette = "Pastel1")+
  labs(title = "Inhibitory Interneurons")+
  coord_cartesian(xlim = c(0.15,1))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p6 = ggplot(data=data_Ex)+
  geom_point(size=2,aes(x=accessibility_score,y=coaccessibility_score,color=Time,shape=Time),alpha=0.7)+
  scale_color_brewer(palette = "Pastel1")+
  labs(title = "Excitatory Interneurons")+
  #coord_cartesian(xlim = c(0.7,1))+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p3_total = (p4+p1)/(p5+p2)/(p6+p3)
ggsave("~/p3.jpg",p3_total,width = 12,height = 14,dpi = 300)



data_total = rbind(data_InvEx,data_InvEx2)
ggplot(data=data_total,aes(x=group,y=coaccessibility_score,color=Time))+
  geom_violin(aes(fill=Time))+
  geom_boxplot(alpha=0.6, outlier.size=0, size=0.9, width=0.8)+
  #geom_jitter(aes(fill=Time),position=position_jitter(0.05), size=1, alpha=0.7)+
  scale_fill_brewer(palette = "Pastel1")+
  #geom_line(aes(group=set) ,size=0.8,colour="#9C9C9C",linetype=2)+
  labs(x="Excitatory vs Inhibitory LINE peak", y="Coaccessibility Score")+theme_bw()


gene_inhibitory = Find_Genes_with_Peak(peak_Inhibitory_adult_LINE,brain_co,0.48)$id
gene_excitatory = Find_Genes_with_Peak(peak_Excitatory_adult_LINE,brain_co,0.48)$id
genes_inhibitory<-mapIds(x=org.Mm.eg.db,
                   keys = gene_inhibitory,
                   keytype = "ENSEMBL",
                   column = "SYMBOL")
genes_inhibitory<-na.omit(genes_inhibitory)


genes_excitatory<-mapIds(x=org.Mm.eg.db,
                         keys = gene_excitatory,
                         keytype = "ENSEMBL",
                         column = "SYMBOL")
genes_excitatory<-na.omit(genes_excitatory)

vd <- Venn_four(
  lst.A <- genes_excitatory,
  lst.B <- genes_inhibitory,
  lst.C <- Excitatory.markers,
  lst.D <- Inhibitory.markers
)

plot(vd,
     labels = list(
       labels = c('Excitatory_LINE', 'Inhibitory_LINE', 'Excitatory_Marker','Inhibitory_Marker'),
       col = "gray20", font = 2
     ), 
     edges = list(col="gray60", lex=1),
     fills = list(fill = c("#297CA0", "#BBBBBB", "#E9EA77","#B8B8FF"), alpha = 0.6),
     quantities = list(cex=.8, col='gray20')
)



