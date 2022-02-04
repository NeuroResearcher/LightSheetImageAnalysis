pacman::p_load(rstatix)

#Import the BrainsDataFunction 
source(file = 'BrainsDataFunction.R')

#Import the RegionsTblFunction
source(file = 'RegionsTblFunction.R')

# Using the BrainsDataFunction from line 4, a dataframe is imported into the 
# environment containing all data from a given directory as indicated by the  
# .path argument in the function

#Create a list of dataframes of all sample data for cFos stains
CFosGFP<-BrainsData(.path = 'CFosGFP')

#Create a list of dataframes of all sample data for NeuN stains
NeuN<-BrainsData(.path = 'NeuN')

#Create a list of dataframes of all sample data for co-expression of NeuN and cFosGFP
Coexpression<-BrainsData(.path = 'Coexpression')

#Create a list of dataframes that have whole-brain data as output from the 
#BrainsDataFunction, and annotate each experimental treatmet group based on the 
#BrainId column, which corresponds to the individual csv file from which each sample's data
#was recorded
StainResults<-lst(CFosGFP$WholeBrain, NeuN$WholeBrain, Coexpression$WholeBrain) %>%
  map(function (d){
    d %>%
      #some brains were statistical outliers and were removed analysis
      filter(!str_detect(BrainId,'FosGFP_FVB[1579]'))%>%
      mutate(BrainId = str_remove(BrainId, '_CFosGFP$|_NeuN$|_Coexpression$')) %>%
      mutate(Treatment = as_factor(case_when(
               str_detect(BrainId, 'FVB[15]|0723_3|June23_[345]') ~ 'HT',
               str_detect(BrainId, 'FVB7|(June|07)23_[12]') ~ 'RT',
               str_detect(BrainId, 'FVB9|0723_[5678]') ~ 'LPS + HT'))) %>%
      
      #any sample data that was in the directory from which original dataframes 
      #were imported, but is not part of this study, is now removed from the 
      #dataframe
      filter(! is.na(Treatment))
  })


#Each experimental group may have a differential expression of cFosGFP or NeuN
#It is important to look at the order of brain regions enriched in either marker
#RegionsTbl function is used to look at each experimental group for each stain
RegionsTbls<-StainResults %>%
  map(function (d){
    TF<- d %>% filter(str_detect(Treatment,'^HT$')) %>% RegionsTbl()
    
    IF<- d %>% filter(str_detect(Treatment, 'LPS')) %>% RegionsTbl()
    
    RT<- d %>% filter(str_detect(Treatment, 'RT')) %>% RegionsTbl()
     
    return(lst(TF, IF, RT))
  })

#Look at all regions with significant variance across experimental groups for
#cell density of markers of interest
AnovaOutput<- StainResults %>%
  map(function (d){
    d %>%
      group_by(Name) %>%
      anova_test(formula = DensityCellsMm3 ~ Treatment) %>%
      as_tibble() %>%
      filter(p<0.05)
  })

#Run a TukeyHSD postHoc for the same linear model as ANOVA
TukeyHsdOutput<- StainResults %>%
  map(function (d){
    d %>%
      group_by(Name) %>%
      tukey_hsd(formula = DensityCellsMm3 ~ Treatment) %>%
      filter(! str_detect(p.adj.signif, 'ns'))
  })


openxlsx::write.xlsx(x = AnovaOutput, 'Stats/FeverConditionAnovas.xlsx', overwrite = T)
openxlsx::write.xlsx(x = TukeyHsdOutput, 'Stats/FeverConditionTukeyHSD.xlsx', overwrite = T)
openxlsx::write.xlsx(x = RegionsTbls, 'Stats/FeverConditionRegionAbundance.xlsx', overwrite = T)

saveRDS(TukeyHsdOutput, 'WrittenData/TukeyHSDOutputs.RDS')
saveRDS(CFosGFP, 'WrittenData/CFosGFP.RDS')
saveRDS(NeuN, 'writtenData/NeuN.RDS')
saveRDS(Coexpression, 'writtenData/Coexpression.RDS')

#Write an RDS file of the munged cFosGFP, NeuN, and co-labeled dataframes so that
#the they can be used in a script for visualization of specific results
saveRDS(StainResults, 'WrittenData/AllStainsWholeBrainResults.RDS')
