args<-commandArgs(TRUE)
#args[1]=guide_list
#args[2]=cdss
#args[3]=outfile
# setwd('~/Documents/deskgen_projects/muc_library')
# args=c('cnio_all_guides.txt', 'cnio_genelist_exon_scores.txt', 'all_guides_cds_pos.txt')

get_cds_positions <- function(guides_file, cdss_file) {
  # get positions within concatenated cds
  guides <- read.table(file=guides_file, header=T, sep="\t")
  cdss <- read.table(file=cdss_file, header=F, sep="\t")
  cds_pos = data.frame()
  pos=ncol(guides)+1
  for (gene in levels(guides[,3])) {
    gene_guides = guides[which(guides[,3] == gene),]
    gene_guides[,pos] = -1
    gene_cdss = cdss[which(cdss[,1] == gene),]
    for (guide in 1:nrow(gene_guides)) {
      if (gene_guides[guide,6] == 1) {
        cut_site = gene_guides[guide,5] + 21
      } else {
        cut_site = gene_guides[guide,5] + 9
      }
      gene_cdss[,8] = gene_cdss[,4]-gene_cdss[,3]+1
      for (cds in 1:nrow(gene_cdss)) {
        if (cut_site >= gene_cdss[cds,3] & cut_site <= gene_cdss[cds,4]) {
          gene_guides[guide,pos] = (cut_site-gene_cdss[cds,3] + sum(gene_cdss[1:cds-1,8]))/sum(gene_cdss[,8])
        }
        if (cut_site > gene_cdss[nrow(gene_cdss),4]) {
          gene_guides[guide,pos] = 1
        }
        if (cut_site < gene_cdss[1,3]) {
          gene_guides[guide,pos] = 0
        }
        #if (gene_guides[guide,11] == -1 & gene_guides[guide,12] >= 0) {
        #  gene_guides[guide,12] = 1-gene_guides[guide,12]
        #}
      } }
    for (guide in 1:nrow(gene_guides)) {
      if (gene_guides[guide,pos-1] == -1 & gene_guides[guide,pos] >= 0) {
        gene_guides[guide,pos] = 1-gene_guides[guide,pos] }
    }
    gene_guides[is.na(gene_guides)] <- 0
    names(gene_guides)[pos] = 'cds_pos'
    if (nrow(cds_pos) == 0) {
      cds_pos = gene_guides
    } else {
      cds_pos = rbind(gene_guides, cds_pos)
    }
  }
  return(cds_pos)
}

cds_pos_all <- get_cds_positions(args[1], args[2])
write.table(cds_pos_all, file=args[3], sep="\t", row.names=F, col.names=T, quote=F)

