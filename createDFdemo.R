# function to aggregate and load the data for the MSYNC project (c) IS Plank
# INPUT: 
#   * df.indi: data frame with one row per subject, must contain the columns listed in ls.vars and grp.var
#   * ls.vars: list of column names in df.indi to be included in overview; referenced columns must contain doubles
#   * grp.var: column name for the grouping variable; must contain only two different values based on which comparisons are computed
#
# OUTPUT: 
#   * df.demo: overview containing one row per variable in ls.vars and the columns measuremend, first group, second group and the log.bf

createDFdemo = function(df.indi, ls.vars, grp.var) {
  
  # load library
  library(tidyverse)
  library(BayesFactor)
  
  # add group column which is the same as the grp.var column
  df.indi$group = df.indi[[grp.var]]
  
  # get the group options
  ls.grp = unique(df.indi$group)
  
  # convert the measures to long which we include in the participant table
  df = df.indi %>%
    select(subID, group, all_of(c(ls.vars, grp.var))) %>%
    pivot_longer(cols = where(is.numeric)) %>%
    mutate_if(is.character, as.factor)
  
  # initialise the data frame for posthoc tests
  df.demo = data.frame()
  
  # now we loop through our measurements to create our demographics table
  for (m in unique(df$name)) {
    # select the relevant part of df.sub
    df.rel = df %>% filter(name == m)
    # check which of the group's data is not normally distributed
    df.sht = df.rel %>% 
      group_by(group) %>%
      shapiro_test(value) %>%
      filter(p < 0.05)
    # if more than zero is not normally distributed...
    if (nrow(df.sht) > 0) {
      # rank transform the data
      df.rel = df.rel %>% ungroup() %>% mutate(value = rank(value))
    }
    # compute the ANOVA
    aov = anovaBF(value ~ group, data = df.rel)
    # get back the original, untransformed values 
    df.rel = df %>% filter(name == m)
    # put all the information into the demographics table
    for (g in ls.grp) {
      df.demo = rbind(df.demo, 
                      data.frame(
                        measurement = m,
                        name        = g,
                        value       = sprintf("%.2f ±%.2f (%.0f to %.0f)", 
                                                mean(df.rel[df.rel$group == g,]$value, na.rm = T), 
                                                sd(df.rel[df.rel$group   == g,]$value, na.rm = T), 
                                                min(df.rel[df.rel$group  == g,]$value, na.rm = T), 
                                                max(df.rel[df.rel$group  == g,]$value, na.rm = T)
                        ),
                        bf.log      = aov@bayesFactor[["bf"]])
                      )
    }
  }
  
  df.demo = df.demo %>%
    pivot_wider() %>% relocate(bf.log, .after = last_col())
  
  return(df.demo)
  
}