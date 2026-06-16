# Adjustments required: 1. The satellite_gap region needs to be adjusted based on satellite modifications.
# 2. Reorder cultivated rice accessions by swapping HZ and NIP.

gapte$V3 <- factor(gapte$V3, levels=c('XX_Lhex_hap4', 'XX_Lhex_hap3', 'XX_Lhex_hap2', 'XX_Lhex_hap1', 'GG_Omey_hap2', 'GG_Omey_hap1', 'FF_Obra_hap2', 'FF_Obra_hap1', 'EE_Oaus_hap2', 'EE_Oaus_hap1', 'CC_Ooff_hap2', 'CC_Ooff_hap1', 'BB_Opun_hap2', 'BB_Opun_hap1', 'AA_Olon_hap2', 'AA_Olon_hap1', 'AA_Oglu_hap2', 'AA_Oglu_hap1', 'AA_Ogla_hap2', 'AA_Ogla_hap1', 'AA_Oniv_hap2', 'AA_Oniv_hap1', 'AA_Oruf_hap2', 'AA_Oruf_hap1', 'AA_Osat_ind', 'AA_Osat_jap'))

gapte$V1 <- factor(gapte$V1, levels = c('Helitron', 'TIR','SINE','LINE','LTR'))

ggplot(gapte, aes(x=V3,y=V2, fill=V1))+
    geom_bar(stat = "identity", position = "stack")+
    coord_flip()+
    scale_fill_manual(values = c('#9BC985','#F7D58B','#B595BF','#982C2C','#159FD7'))+
    labs(y='Percentage', x='material',fill='type')+
    theme_light()