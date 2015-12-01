source('~/Desktop/MacArthur/de-novo-splice/analysis2/barplot_functions.R')

# orange is for exonic decoys
# periwinkle is for intronic decoys

dist_cols = c("IntronicDist", "ProximalDist", "DistalDist", "ExonicDist")

do_barplot <- function(cols, sre, colors, dens=T) {
  intron = subset(use, intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols)
  exon = subset(use, !intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols)
  m_bar1 = get_plot_data2(exon, c("Distal", "Proximal", "Intronic"), sre, pois=dens)
  m_bar2 = get_plot_data2(intron, c("Exonic", "Proximal", "Distal"), sre, pois=dens)
  plot_bars(m_bar1, m_bar2, sprintf('%s comparison', sre), xax=c("up", "middle", "down"), colors=colors)
}

weighted = F
use_density = F
# --- ESS ---
if (weighted) {
  cols = c(dist_cols, "IntronicESSWeighted", "ProximalESSWeighted", "DistalESSWeighted", "ExonicESSWeighted")
  sre = "ESSWeighted"
} else {
  cols = c(dist_cols, "IntronicESS", "ProximalESS", "DistalESS", "ExonicESS")
  sre = "ESS"
}
do_barplot(cols, sre, c('#FA9911', '#A9A1FF'), dens=use_density)

# --- ESE  ---
if (weighted) {
  cols = c(dist_cols, "IntronicESEWeighted", "ProximalESEWeighted", "DistalESEWeighted", "ExonicESEWeighted")
  sre = "ESEWeighted"
} else {
  cols = c(dist_cols, "ExonicESE", "IntronicESE", "ProximalESE", "DistalESE")
  sre = "ESE"
}
do_barplot(cols, sre, c('#FA9911', '#A9A1FF'), dens=use_density)

# --- G-triplets ---
cols = c(dist_cols, "ExonicGtrip", "ProximalGtrip", "DistalGtrip", "IntronicGtrip")
sre = "Gtrip"  
do_barplot(cols, sre, c('#FA9911', '#A9A1FF'), dens=use_density)
