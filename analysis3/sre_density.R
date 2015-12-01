# source('~/Desktop/MacArthur/de-novo-splice/analysis2/weight_hexamers.R')
source('~/Desktop/MacArthur/de-novo-splice/analysis2/plot_functions.R')

weighted = T
# aberrant = cryptic_feats
# aberrant = novo_feats
aberrant = rbind(novo_feats, cryptic_feats)

# --- ESS for exonic alternative ---
if (weighted) {
  cols = c("IntronicDist", "ProximalDist", "DistalDist", "IntronicESSWeighted", "ProximalESSWeighted", "DistalESSWeighted")
  sre = "ESSWeighted"
} else {
  cols = c("IntronicDist", "ProximalDist", "DistalDist", "IntronicESS", "ProximalESS", "DistalESS")
  sre = "ESS"
}
subset(aberrant, !intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols) %>% colSums -> novo
subset(use, !intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols) %>% colSums -> decoy
sections = c("Distal", "Proximal", "Intronic")
m_bar = get_plot_data(novo, decoy, sections, sre)
plot_bars(m_bar[,1:2], m_bar[,3:4], sprintf('%s for exonic alternative', sre), xax=sections)

# --- ESS for intronic alternative ---
if (weighted) {
  cols = c("ProximalDist", "DistalDist", "ExonicDist", "ExonicESSWeighted", "ProximalESSWeighted", "DistalESSWeighted")
  sre = "ESSWeighted"  
} else {
  cols = c("ProximalDist", "DistalDist", "ExonicDist", "ExonicESS", "ProximalESS", "DistalESS")
  sre = "ESS"
}
subset(aberrant, intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols) %>% colSums -> novo
subset(use, intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols) %>% colSums -> decoy
sections = c("Exonic", "Proximal", "Distal")
m_bar = get_plot_data(novo, decoy, sections, sre)
plot_bars(m_bar[,1:2], m_bar[,3:4], sprintf('%s for intronic alternative', sre), xax=sections)

# --- ESE for exonic alternative ---
if (weighted) {
  cols = c("ProximalDist", "DistalDist", "IntronicDist", "IntronicESEWeighted", "ProximalESEWeighted", "DistalESEWeighted")
  sre = "ESEWeighted"
} else {
  cols = c("ProximalDist", "DistalDist", "IntronicDist", "IntronicESE", "ProximalESE", "DistalESE")
  sre = "ESE"
}
subset(aberrant, !intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols) %>% colSums -> novo
subset(use, !intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols) %>% colSums -> decoy
sections = c("Distal", "Proximal", "Intronic")
m_bar = get_plot_data(novo, decoy, sections, sre)
plot_bars(m_bar[,1:2], m_bar[,3:4], sprintf('%s for exonic alternative', sre), xax=sections)

# --- ESE for intronic alternative ---
if (weighted) {
  cols = c("ProximalDist", "DistalDist", "ExonicDist", "ExonicESEWeighted", "ProximalESEWeighted", "DistalESEWeighted")
  sre = "ESEWeighted"
} else {
  cols = c("ProximalDist", "DistalDist", "ExonicDist", "ExonicESE", "ProximalESE", "DistalESE")
  sre = "ESE"  
}
subset(aberrant, intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols) %>% colSums -> novo
subset(use, intronic & AltMES - RefMES >= -5 & AltMES - RefMES <= 5, select=cols) %>% colSums -> decoy
sections = c("Exonic", "Proximal", "Distal")
m_bar = get_plot_data(novo, decoy, sections, sre)
plot_bars(m_bar[,1:2], m_bar[,3:4], sprintf('%s for intronic alternative', sre), xax=sections)

