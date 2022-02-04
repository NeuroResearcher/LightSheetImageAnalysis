pacman::p_load(tidyverse)

BrainsData<- function(.path){
  #Import all csv data from a given directory in the computer
  #Helps to organize your files in the directory so that only the files of interest are in
  #that given directory; however, you could always filter out data associated with a file
  #that you do not intend to analyze.
  ExpressionData <-rio::import_list(
      dir(
        path.expand(.path),
        pattern = '.csv',
        full.names = T),
      rbind = T,
      setclass = 'tbl_df') %>%
      select(c(1,2,6:9)) %>%
      janitor::clean_names('upper_camel') %>%
    mutate(File = str_remove(File,'^.*/')) %>%
    mutate(File = str_remove(File,'\\.csv$')) %>%
    
    # 'direct tectospinal pathway' was only measured on the right
    # hemisphere, so it's being removed from the dataset because
    # we can't analyze both sides of it
    
    filter(!+str_detect(Name, 'background|direct tectospinal pathway')) %>%
    mutate(Name = as.factor(Name),
           Id = as.factor(Id),
           File = str_remove(File, '_GFP_RFP')) %>%
    mutate(File = as.factor(str_remove(File, '.csv$'))) %>%
    # order the dataset by id so that the first half of the dataframe
    # is a mirror image of the second half when hemispheres are removed
    arrange(Id) %>%
    rename(BrainId = File,
           RegionId=Id)
  
  #Make a dataset for the left hemisphere values
  LeftBrain<-slice_head(ExpressionData, n=.5*nrow(ExpressionData))
  
  #Make a dataset for the right hemisphere values
  RightBrain<-slice_tail(ExpressionData, n=.5*nrow(ExpressionData))
  
  # Convert the original dataframe so that it equals the combined
  # counts, volumes, and densities from the 
  WholeBrain<-LeftBrain %>%
    mutate(Count = Count + RightBrain$Count,
           Name = as.character(Name),
           RegionId = as.numeric(RegionId),
           VolumeMm3 = (VolumeMm3 + RightBrain$VolumeMm3)) %>%
    mutate(DensityCellsMm3 = Count/VolumeMm3) %>%
    select(2:6) %>%
    mutate(Name = str_remove(Name, 'left ')) %>%
    filter(str_detect(Name, '^root$|^Basic cell groups and regions$|^Cerebrum$') ==F) %>%
    mutate(Name = as_factor(Name))
  
  
  return(lst(WholeBrain,LeftBrain,RightBrain))
}
