#args<-commandArgs(TRUE)
#args[1]=guide_list
#args[2]=bs_list
#args[3]=cds_list
setwd('~/Documents/deskgen_projects/guide_picker_plot')
args=c('smarca4_guidelist_cds_pos.txt', 'muc_5_blackswans.txt', 'Smarca4_protein_exon_scores.txt')


#bs_guides <- read.table(file='cnio_blackswans.txt', header=T, sep="\t")
#cdss <- read.table(file='cnio_genelist_exon_scores.txt', header=F, sep="\t")
#all_guides <- read.table(file=args[1], header=T, sep="\t")

guides_all <- read.table(file=args[1], header=T, sep="\t")
#guides_bs <- read.table(file=args[2], header=T, sep="\t")
cdss <- read.table(file=args[3], header=F, sep="\t")
#gene_list <- levels(guides_all[,2])
#pos = ncol(guides_all)

cds_pos_all = guides_all
#guides_bs = data.frame()
cdss_list = cdss
#gene_list <- levels(guides_all[,1])

draw_gene_graphs <- function(cds_pos_bs, cds_pos_all, cdss_list) {
  for (gene in levels(cds_pos_bs$gene_id)) {
    gene_name = as.character(unique(cds_pos_bs[which(cds_pos_bs$gene_id==gene),2]))
    gene_guides = cds_pos_bs[which(cds_pos_bs$gene_id==gene),]
    all_gene_guides = cds_pos_all[which(cds_pos_all$gene_id==gene),]
    gene_cdss = cdss_list[which(cdss_list[,1] == gene),]
    gene_cdss[,8] = gene_cdss[,4]-gene_cdss[,3]+1
    if (gene_cdss[1,5] == 1) {
      bars = gene_cdss[,8]
      scores = gene_cdss[,7]
    } else {
      bars = rev(gene_cdss[,8])
      scores = rev(gene_cdss[,7])
    }
    # visualise per gene guides
    pdf(paste('gene_pdfs/', gene_name,'_guides.pdf', sep=''), width=15, height=5)
    par(mfrow=c(1,1))
    par(mar=c(6,6,4,4))
    plot(all_gene_guides$cds_pos, all_gene_guides$doench, cex=1,
         lwd=2, col='azure4', pch=18, main=paste(gene_name,' guides'),
         xlab="5'-> 3' along all CDSs", ylab='Doench score', xaxt='n',
         cex.axis=1, cex.lab=1, xlim=c(0,1), ylim=c(0,1))
    par(new=T)
    opar <- par(lwd = 1)
    barplot(scores, bars, pch=17, space=0,
            lwd=1, axes=F, col='black', border='white', cex=1,
            #xlab="Doench Score", ylab='Number of Guides (In Lethal 100 Coding Regions)',
            cex.axis=1, cex.lab=1, ylim=c(0,10), xlim=c(0,sum(gene_cdss[,8])))
    axis(1,lwd.ticks=1,lwd=1.5, cex.axis=1, xpd=F)
    #axis(4,lwd.ticks=1,lwd=2,cex.axis=1, xpd=F)
    mtext(" ^  # transcripts", side=4, line=0.5, adj=0, cex=1.2)
    par(new=T)
    plot(gene_guides$cds_pos, gene_guides$doench, cex=2,
         lwd=2, col='black', bg='cyan3', pch=23, yaxt='n', ylab='', xlab='',
         xaxt='n', cex.axis=1, cex.lab=1, xlim=c(0,1), ylim=c(0,1))
    dev.off()
  }
}

draw_one_gene_graph_all_guides <- function(gene, cds_pos_all, cdss_list) {
  cdss_list[,8] = cdss_list[,4]-cdss_list[,3]+1
  if (cdss_list[1,5] == 1) {
    bars = cdss_list[,8]
    scores = cdss_list[,7]
  } else {
    bars = rev(cdss_list[,8])
    scores = rev(cdss_list[,7])
  }
  pdf(paste('gene_pdfs/', gene,'_guides.pdf', sep=''), width=10, height=5)
  par(mfrow=c(1,1))
  par(mar=c(6,6,4,4))
  plot(cds_pos_all$cds_pos, cds_pos_all$doench, cex=1,
       lwd=2, col='azure4', pch=18, main=paste(gene,' guides'),
       xlab="5'-> 3' along all CDSs", ylab='Doench score', xaxt='n',
       cex.axis=1, cex.lab=1, xlim=c(0,1), ylim=c(0,1))
  par(new=T)
  opar <- par(lwd = 1)
  barplot(scores, bars, pch=17, space=0,
          lwd=1, axes=F, col='black', border='white', cex=1,
          #xlab="Doench Score", ylab='Number of Guides (In Lethal 100 Coding Regions)',
          cex.axis=1, cex.lab=1, ylim=c(0,20), xlim=c(0,sum(cdss_list[,8])))
  axis(1,lwd.ticks=1,lwd=1.5, cex.axis=1, xpd=F)
  #axis(4,lwd.ticks=1,lwd=2,cex.axis=1, xpd=F)
  mtext(" ^  # transcripts", side=4, line=0.5, adj=0, cex=1.2)
  dev.off()
}

draw_gene_graphs(guides_bs, guides_all, cdss)

cdss_list[,7] = 1
draw_one_gene_graph_all_guides('Smarca4', cds_pos_all, cdss_list)

# On-target vs off-target scatter plot
# aggregate off-target score (coding off-target X10 + non-coding)
cds_pos_all[,16] = (cds_pos_all$coding_offt*10 + cds_pos_all$genomic_offt)/2
cds_pos_all[,17] = cds_pos_all$coding_offt*10 + cds_pos_all$genomic_offt
for (i in 1:nrow(cds_pos_all)) {
  if (cds_pos_all[i,16] > 1) {
    cds_pos_all[i,16] <- 1.5
  } 
}
names(cds_pos_all)[16] <- 'aggregate_offt'
gene = 'Smarca4'

pdf(paste('gene_pdfs/', gene,'_gpp.pdf', sep=''), width=7, height=7)
  par(mfrow=c(1,1))
  par(mar=c(5,5,3,3))
  plot(cds_pos_all$aggregate_offt, cds_pos_all$doench,
       xlim = c(0,1), ylim = c(0,1),
       xlab = 'Off-Target Aggregated Score', ylab = 'Doench On-Target Score',
       pch = 18, col = 'azure4', cex = 2)
dev.off()

quantile(cds_pos_all[,17], probs=c(0, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
