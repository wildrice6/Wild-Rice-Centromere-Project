ltr$V3 <- factor(ltr$V3, levels=c('XX_Lhex_hap4', 'XX_Lhex_hap3', 'XX_Lhex_hap2', 'XX_Lhex_hap1', 'GG_Omey_hap2', 'GG_Omey_hap1', 'FF_Obra_hap2', 'FF_Obra_hap1', 'EE_Oaus_hap2', 'EE_Oaus_hap1', 'CC_Ooff_hap2', 'CC_Ooff_hap1', 'BB_Opun_hap2', 'BB_Opun_hap1', 'AA_Olon_hap2', 'AA_Olon_hap1', 'AA_Oglu_hap2', 'AA_Oglu_hap1', 'AA_Ogla_hap2', 'AA_Ogla_hap1', 'AA_Oniv_hap2', 'AA_Oniv_hap1', 'AA_Oruf_hap2', 'AA_Oruf_hap1', 'AA_Osat_ind', 'AA_Osat_jap'))

genome_colors <- c(
    "AA" = "#76D273",
    "BB" = "#F2AE2C",
    "CC" = "#684E94",
    "EE" = "#4E84C3",
    "FF" = "#D55F6F",
    "GG" = "#9D5427",
    "XX" = "#595959"
)

ggplot(ltr, aes(x=V1, y=V3, color=V4))+
    geom_point(aes(size=V2), alpha=0.8)+
    scale_color_manual(values = genome_colors)+
    theme_minimal() +
    labs(x='LTR type', y='Materials', color='Genome type', size='Proportion')+
    theme(axis.text.x=element_text(angle=30, hjust=1))