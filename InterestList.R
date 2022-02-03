#Load in necessary packages for script function
pacman::p_load(tidyverse, ggpubr, rstatix)

#Load data results from light-sheet imaging with  BrainId, experimental treatment,
#and anitgen density in the dataframe
StainResults<-readRDS('WrittenData/AllStainsWholeBrainResults.RDS')


#Generate a list of brain regions of interest by name
#Each of the following regions will be taken as a subset of all global regions from the loaded dataframe

InterestList<-c('cerebrum related', 'Cerebral cortex', 'Isocortex', 
                'Frontal pole, cerebral cortex','Somatomotor areas', 'Primary motor area', 
                'Secondary motor area', 'Somatosensory areas','Primary somatosensory area', 
                'Primary somatosensory area, barrel field, layer 2/3',
                'Primary somatosensory area, barrel field, layer 4', 'Auditory areas', 
                'Primary auditory area', 'Visual areas', 'Laterointermediate area',
                'Anterior cingulate area', 'Prelimbic area', 'Infralimbic area', 
                'Orbital area', 'Agranular insular area', 'Retrosplenial area', 
                'Temporal association areas', 'Postrhinal area', 'Perirhinal area', 
                'Entorhinal area')

#The imported list of dataframes from line 5 needs to be subsetted
#Use map as a pipe-friendly way to apply the same function to all dataframes in a list
#Function is a custom function of input 'd' being each dataframe from the imported list
#Custom function and map returns subsetted list of dataframes containing only data for
  #brain regions in the list of areas of interest

InterestListStainsData<-StainResults %>%
  map(function(d){
    d %>%
      filter(Name %in% InterestList) %>%
      dplyr::group_by(Name)
  })



#Run ANOVA on each subsetted dataframe in the list to see if there is significant variance of Cell Density across experimental groups
#Because each dataframe is a grouped datafame, an ANOVA is run for each 'Name' group or brain region group in the dataframe

AnovaResults<- InterestListStainsData %>%
  map(rstatix::anova_test, formula = DensityCellsMm3 ~ Treatment)

#Use the same linear model from the ANOVA to run a TukeyHSD posthoc to analyze pairwise differences between experimental groups
TukeyResults<- InterestListStainsData %>%
  map(rstatix::tukey_hsd, formula = DensityCellsMm3 ~ Treatment) %>%
  map(mutate, y.position = 25000)



openxlsx::write.xlsx(x = AnovaResults, file='Stats/GlobalListAnova.xlsx', createWorkbook = T)
openxlsx::write.xlsx(x = TukeyResults, file='Stats/GlobalListTukeyHsd.xlsx', createWorkbook = T)


#Visualize group trends/differences in cFosGFP staining for all areas of interest
ggbarplot(data = InterestListStainsData$`CFosGFP$WholeBrain`, 
          x =  'Treatment',
          y = 'DensityCellsMm3',
          facet.by = 'Name',
          add = c('mean_se', 'jitter'),
          add.params = list(shape = 1, size = .9),
          color = 'Treatment',
          ylab = 'Cells / mm^3',
          width = .65,
          ylim = c(0,63000),
          title = 'FosGFP') + 
  stat_pvalue_manual(data = TukeyResults$`CFosGFP$WholeBrain`,
                     hide.ns = T,
                     step.increase = 0.075,
                     bracket.size = 0.01,
                     tip.length = 0)

#Visualize group trends/differences in NeuN staining for all areas of interest
ggbarplot(data = MichelleListStainsData$`NeuN$WholeBrain`, 
          x =  'Treatment',
          y = 'DensityCellsMm3',
          facet.by = 'Name',
          add = c('mean_se', 'jitter'),
          add.params = list(shape = 1, size = .9),
          color = 'Treatment',
          ylab = 'Cells / mm^3',
          width = .65,
          ylim = c(0,75000),
          title = 'NeuN') + 
  stat_pvalue_manual(data = TukeyResults$`NeuN$WholeBrain`,
                     hide.ns = T,
                     step.increase = 0.075,
                     bracket.size = 0.01,
                     tip.length = 0)

#Visualize group trends/differences in cFosGFP+NeuN+ co-labeled staining for all areas of interest
ggbarplot(data = MichelleListStainsData$`Coexpression$WholeBrain`, 
          x =  'Treatment',
          y = 'DensityCellsMm3',
          facet.by = 'Name',
          add = c('mean_se', 'jitter'),
          add.params = list(shape = 1, size = .9),
          color = 'Treatment',
          ylab = 'Cells / mm^3',
          width = .65,
          ylim = c(0,35000),
          title = 'Coexpression') + 
  stat_pvalue_manual(data = TukeyResults$`Coexpression$WholeBrain`,
                     hide.ns = T,
                     step.increase = 0.075,
                     bracket.size = 0.01,
                     tip.length = 0)
