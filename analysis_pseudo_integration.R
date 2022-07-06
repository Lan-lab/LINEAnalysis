Excitatory_E18 <-subset(E18_brain,subset=group=="A:Excitatory neurons")

Excitatory_E18_pseudo = Excitatory_E18
DefaultAssay(Excitatory_E18_pseudo) = "RNA"

transfer.anchors <- FindTransferAnchors(
  reference = E18_seurat,
  query = Excitatory_E18_pseudo,
  reduction = 'cca',
  dims = 1:30
)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = E18_seurat$celltype,
  weight.reduction = Excitatory_E18_pseudo[['lsi']],
  dims = 2:30
)
id = grep("neuron",Excitatory_E18_pseudo$predicted.id)
Excitatory_E18_pseudo <- AddMetaData(object = Excitatory_E18_pseudo, metadata = predicted.labels)
Excitatory_E18_pseudo = Excitatory_E18_pseudo[,colnames(Excitatory_E18_pseudo)[id]]

Excitatory_E18_pseudo$cell_allen = (Excitatory_E18$predicted.id)[id]









