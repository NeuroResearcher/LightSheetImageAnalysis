pacman::p_load(tableone, tidyverse)


# Make a summary dataset to show descending order of D1 D2 receptor
# coexpression

# First step is to make a table and convert it to tibble
RegionsTbl<- function(data){
  RegionsTbl<- as_tibble(print(CreateTableOne(data = data,
                                              vars = 'DensityCellsMm3',
                                              strata = 'Name',
                                              test = F))) %>%
    # Make wide tibble long                 
    pivot_longer(cols = everything()) %>%
  
    # means and SD's are in the bottom half of the dataset, so subset that
    slice_tail(n= 0.5 * nrow(.)) %>%
  
    # add column of n values back in (easier to add this because constant value)
    mutate(n = n_distinct(data$BrainId)) %>%
  
    #make mean and SD into different columns
    separate(col = value, into = c('mean', 'SD'), sep = '\\(') %>%
  
    # remove parenthesis from SD column so that we can convert to numeric
    mutate(SD = as.numeric(str_remove(SD, '\\)')),
         mean = as.numeric(mean)) %>%
  
    # formally change name to 'Region' b/c it's a brain region
    rename(Region = name) %>%
  
    # arrange by descending mean value b/c the goal is to see where 
    # coexpression is highest in the brain
    arrange(desc(mean))%>%
  
    # remove all 0 values b/c we don't care about regions that don't express
    filter(mean > 0) %>%
  
    # add a computed column of std error of the mean
    mutate(SEM = SD/sqrt(n))
  
    return(RegionsTbl)
}