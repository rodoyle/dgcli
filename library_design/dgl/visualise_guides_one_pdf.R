args<-commandArgs(TRUE)
#args[1]=guide_list
#args[2]=bs_list
#args[3]=cds_list
#args[4]=gene_list
#setwd('~/Documents/deskgen_projects/mouse_library')
#args=c('all_guides_cds_pos.txt', 'cnio_blackswans_36k_150423.txt', 'cnio_genelist_exon_scores.txt', 'cnio_genelist.txt')

#bs_guides <- read.table(file='cnio_blackswans.txt', header=T, sep="\t")
#cdss <- read.table(file='cnio_genelist_exon_scores.txt', header=F, sep="\t")
#all_guides <- read.table(file=args[1], header=T, sep="\t")

guides_all <- read.table(file=args[1], header=T, sep="\t")
guides_bs <- read.table(file=args[2], header=T, sep="\t")
cdss <- read.table(file=args[3], header=F, sep="\t")
gene_list <- read.table(file=args[4], header=F, sep="\t", stringsAsFactors=FALSE)

draw_gene_graphs <- function(cds_pos_bs, cds_pos_all, cdss_list, gene_list) {
  graph_length = round(length(levels(cds_pos_bs[,2]))/2)
  pdf('gene_pdfs/all_guides.pdf', width=16, height=5*graph_length)
  par(mfrow=c(graph_length,2))
  for (gene in levels(cds_pos_bs[,2])) {
    gene_name = gene_list[which(gene_list[,2]==gene),1]
    gene_guides = cds_pos_bs[which(cds_pos_bs[,2]==gene),]
    all_gene_guides = cds_pos_all[which(cds_pos_all[,2]==gene),]
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
    par(mar=c(5,5,4,4))
    plot(all_gene_guides[,12], all_gene_guides[,7], cex=1,
         lwd=2, col='azure4', pch=18,
         xlab="5'-> 3' along all CDSs", ylab='Doench score', xaxt='n',
         cex.axis=1, cex.lab=1, xlim=c(0,1), ylim=c(0,1))
    par(new=T)
    opar <- par(lwd = 2)
    barplot(scores, bars, pch=17, space=0,
            lwd=2, axes=F, col='black', border='white', cex=1,
            #xlab="Doench Score", ylab='Number of Guides (In Lethal 100 Coding Regions)',
            cex.axis=1, cex.lab=1, ylim=c(0,10), xlim=c(0,sum(gene_cdss[,8])))
    axis(1,lwd.ticks=1,lwd=2, cex.axis=1, xpd=F)
    #axis(4,lwd.ticks=1,lwd=2,cex.axis=1, xpd=F)
    mtext(" ^  # transcripts", side=4, line=0.5, adj=0, cex=1.2)
    par(new=T)
    plot(gene_guides[,12], gene_guides[,7], cex=2,
         lwd=2, col='black', bg='cyan3', pch=23, main=paste(gene_name,' guides'),
         xlab="5'-> 3' along all CDSs", ylab='Doench score', xaxt='n',
         cex.axis=1, cex.lab=1, xlim=c(0,1), ylim=c(0,1))
  }
  dev.off()
}

draw_gene_graphs(guides_bs, guides_all, cdss, gene_list)

