args<-commandArgs(TRUE)
#args[1]=guide_list
#args[2]=outfile
#args[3]=no of best guides
#setwd('~/Documents/deskgen_projects/mouse_library')
#args=c('cnio_all_guides_NGG_CDSpos_filtered.txt', 'cnio_blackswans_150507.txt', 5)

# CRITERIA FOR SELECTING THE ‘BEST’ GUIDES:
#
# G = non-coding off-target hits << LEAST IMPORTANT
# D = on-target Doench score             V
# C = coding off-target hits     <<  MOST IMPORTANT
# 
# find n guides in consensus exon
# if n > 100:
#   top 10% G
#   top 10% D
#   top 10% C
#  then
#   top 25% G
#   top 10% D
#   top 10% C
#  then
#   top 25% G
#   top 25% D
#   top 10% C
# else:           1>>>>>>>>>>>>>>4
#     top 25% G | 25 | 50 | 50 | 50
#     top 25% D | 25 | 25 | 50 | 50
#     top 25% C | 25 | 25 | 25 | 50
#
#   then if 2nd consensus exon has score >66% of 1st:
#         do same as previous
#
# else back to consensus exon:
#        top G% | 75 | 75 | 75
#        top D% | 50 | 75 | 75
#        top C% | 50 | 50 | 75
#
# if any are still <5, do by hand

all_guides <- read.table(file=args[1], header=T, sep="\t")
# order guides based on CDS orientation
guides = data.frame()
#genes = levels(all_guides[,2])
genes = naughty_genes
for (gene in genes) {
  gene_guides = all_guides[which(all_guides[,2]==gene),]
  if (gene_guides[1,11] == 1) {
    ordered_guides = gene_guides[order(gene_guides[,4], decreasing=FALSE),]
  } else {
    ordered_guides = gene_guides[order(gene_guides[,4], decreasing=TRUE),]
  }
  if (nrow(guides) == 0) {
    guides = ordered_guides
  } else {
    guides = rbind(guides, ordered_guides)
  }
}

guide_parameters <- function(good_g, exon_g, g, d, c, n) {
  if (nrow(exon_g) == 0) {
    return(list(exon_g, good_g))
  } else {
    index=vector()
    rest=vector()
    for (i in 1:nrow(exon_g)) {
      if (exon_g[i,9] <= g_quant[g] & exon_g[i,9] >= g_quant[g-1] &
            exon_g[i,7] >= d_quant[d] & exon_g[i,7] <= d_quant[d-1] &
            exon_g[i,8] <= c_quant[c] & exon_g[i,8] >= c_quant[c-1] &
            (length(index)+nrow(good_g)) < n) {
        index = c(index,i) } else {rest=c(rest,i)} }
    if (nrow(good_g) == 0) {
      good_g = exon_g[index,]
    } else {
      good_g = rbind(good_g, exon_g[index,])
    }
    exon_g = exon_g[rest,]
    guide_list = list(exon_g,good_g)
    return(guide_list)
  }
}

# get n BLACK SWANS for each gene
black_swans = data.frame()
ten=2
twenty_five=3
fifty=4
seventy_five=5
n=as.numeric(args[3])
for (gene in genes) {
  gene_guides = guides[which(guides[,2]==gene),]
  exon_guides = gene_guides[which(gene_guides[,10]==max(gene_guides[,10])),]
  good_guides = data.frame()
  # calculate quantiles of guides
  g_quant = quantile(gene_guides[,9], probs=c(0, 0.1, 0.25, 0.5, 0.75))
  d_quant = quantile(gene_guides[,7], probs=c(1, 0.9, 0.75, 0.5, 0.25))
  c_quant = quantile(gene_guides[,8], probs=c(0, 0.1, 0.25, 0.5, 0.75))
  ### try progressively top 10-50% of all scores in consensus exon
  guide_list = guide_parameters(good_guides, exon_guides, ten, ten, ten, n)
  exon_guides = data.frame(guide_list[1])
  good_guides = data.frame(guide_list[2])
  if (nrow(data.frame(guide_list[2])) < n) {
    guide_list = guide_parameters(good_guides, exon_guides, twenty_five, ten, ten, n)
    exon_guides = data.frame(guide_list[1])
    good_guides = data.frame(guide_list[2])
  }
  if (nrow(data.frame(guide_list[2])) < n) {
    guide_list = guide_parameters(good_guides, exon_guides, twenty_five, twenty_five, ten, n)
    exon_guides = data.frame(guide_list[1])
    good_guides = data.frame(guide_list[2])
  }
  if (nrow(data.frame(guide_list[2])) < n) {
    guide_list = guide_parameters(good_guides, exon_guides, twenty_five, twenty_five, twenty_five, n)
    exon_guides = data.frame(guide_list[1])
    good_guides = data.frame(guide_list[2])
  }
  if (nrow(data.frame(guide_list[2])) < n) {
    guide_list = guide_parameters(good_guides, exon_guides, fifty, twenty_five, twenty_five, n)
    exon_guides = data.frame(guide_list[1])
    good_guides = data.frame(guide_list[2])
  }
  if (nrow(data.frame(guide_list[2])) < n) {
    guide_list = guide_parameters(good_guides, exon_guides, fifty, fifty, twenty_five, n)
    exon_guides = data.frame(guide_list[1])
    good_guides = data.frame(guide_list[2])
  }
  if (nrow(data.frame(guide_list[2])) < n) {
    guide_list = guide_parameters(good_guides, exon_guides, fifty, fifty, fifty, n)
    exon_guides = data.frame(guide_list[1])
    good_guides = data.frame(guide_list[2])
  }
  ### try exons within 50% of best exon score - look for top 10-50%, else try top 75% of consensus exon
  if (nrow(good_guides) < n) {
  gene_guides = guides[which(guides[,2]==gene),]
  exon_guides = gene_guides[which(gene_guides[,10] < max(gene_guides[,10]) &
                                  gene_guides[,10] >= 0.50),]
  if (nrow(exon_guides) > 0) {
    # calculate quantiles of guides
    #g_quant = quantile(exon_guides[,9], probs=c(0.1, 0.25, 0.5, 0.75))
    #d_quant = quantile(exon_guides[,7], probs=c(0.9, 0.75, 0.5, 0.25))
    #c_quant = quantile(exon_guides[,8], probs=c(0.1, 0.25, 0.5, 0.75))
    guide_list = guide_parameters(good_guides, exon_guides, ten, ten, ten, n)
    exon_guides = data.frame(guide_list[1])
    good_guides = data.frame(guide_list[2])
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, twenty_five, ten, ten, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, twenty_five, twenty_five, ten, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, twenty_five, twenty_five, twenty_five, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, fifty, twenty_five, twenty_five, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, fifty, fifty, twenty_five, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, fifty, fifty, fifty, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
  } else {
    # try top 75% of all scores in consensus exon
    guide_list = guide_parameters(good_guides, exon_guides, seventy_five, fifty, fifty, n)
    exon_guides = data.frame(guide_list[1])
    good_guides = data.frame(guide_list[2])
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, seventy_five, seventy_five, fifty, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, seventy_five, seventy_five, seventy_five, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
  }
  }
  ### check for genes where all guides appear in 3' 1/3 of gene - design two guides in 5' 1/3
  if (min(good_guides[,12]) >= 0.66) {
    good_ordered = good_guides[order(good_guides[,12]),]
    good_guides = good_ordered[1:3,]
    third_guides = gene_guides[which(gene_guides[,12] <= 0.333),]
    exon_guides = third_guides[which(third_guides[,10]==max(third_guides[,10])),]
    # calculate quantiles of guides
    #g_quant = quantile(exon_guides[,9], probs=c(0, 0.1, 0.25, 0.5, 0.75))
    #d_quant = quantile(exon_guides[,7], probs=c(1, 0.9, 0.75, 0.5, 0.25))
    #c_quant = quantile(exon_guides[,8], probs=c(0, 0.1, 0.25, 0.5, 0.75))
    # try progressively top 10-50% of all scores in consensus exon
    guide_list = guide_parameters(good_guides, exon_guides, ten, ten, ten, n)
    exon_guides = data.frame(guide_list[1])
    good_guides = data.frame(guide_list[2])
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, twenty_five, ten, ten, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, twenty_five, twenty_five, ten, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, twenty_five, twenty_five, twenty_five, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, fifty, twenty_five, twenty_five, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, fifty, fifty, twenty_five, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, fifty, fifty, fifty, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, seventy_five, fifty, fifty, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, seventy_five, seventy_five, fifty, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
    if (nrow(data.frame(guide_list[2])) < n) {
      guide_list = guide_parameters(good_guides, exon_guides, seventy_five, seventy_five, seventy_five, n)
      exon_guides = data.frame(guide_list[1])
      good_guides = data.frame(guide_list[2])
    }
  }
  if (nrow(black_swans) == 0) {
    black_swans = good_guides
  } else {
    black_swans = rbind(black_swans, good_guides)
  }
  
}
write.table(black_swans, file=args[2], sep="\t", row.names=F, col.names=T, quote=F)
