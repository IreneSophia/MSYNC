# function to aggregate and load the data for the MSYNC project (c) IS Plank
# Takes as input the path where all the data is stored. No output. 

createMSYNCdata = function(dt.path) {
  
  # load library
  library(tidyverse)
  
  # get demo info for subjects
  df.indi = merge(
    read_csv(file.path(dt.path, "df_centraXX_final.csv")) %>%
      rename("PID" = "internalStudyMemberID") %>% select(-dyad),
    read_csv(file.path(dt.path, "MSYNC_subID_info.csv")) %>%
      select(dyad, PID, order, role, position, subID)) %>%
    rename(
      'enjoyment_MG'  = 'enjoyment',
      'responsive_MG' = 'responsive',
      'continue_MG'   = 'continue'
    ) %>% select(-PID) %>%
    mutate(
      # manually add the CFT iq estimate for one participant >
      # collected in different study in our lab
      CFT_iq = if_else(is.na(CFT_iq), 98, CFT_iq),
      # recode enjoyment > reverse
      enjoyment_MG   = 6 - enjoyment_MG,
      # recode the vision information
      vision = case_match(vision,
                          "Brille" ~ "glasses", 
                          "Kontaktlinsen" ~ "contacts",
                          .default = "none"
                          ),
      gender = tolower(gender)
    ) %>%
    select(subID, dyad, role, order, age, gender, gender_identity, edu, vision, handedness,
           RAADS_final, IRI_final, ADC_final, ends_with("_MG"), ends_with("_conv"), 
           CFT_iq, anx, avoid)
  
  # combine social experience
  df.exp = df.indi %>%
    select(subID | ends_with("_MG") | ends_with("_conv")) %>%
    mutate(
      ## adjust all to the same scale [0 to 1]
      # Mood Cohesien Trust goes from 0 to 10
      close.sMG      = close_MG/10,
      similar.sMG    = similar_MG/10,
      connect.sMG    = connect_MG/10,
      trust.sMG      = trust_MG/10,
      close.sCT      = close_conv/10,
      similar.sCT    = similar_conv/10,
      connect.sCT    = connect_conv/10,
      trust.sCT      = trust_conv/10,
      # IOS 1 to 7
      IOS.sMG        = (IOS_MG-1)/6,
      IOS.sCT        = (IOS_conv-1)/6,
      # rapport 0 to 6
      smooth.sMG     = smooth_MG/6,
      comfort.sMG    = comfort_MG/6,
      smooth.sCT     = smooth_conv/6,
      comfort.sCT    = comfort_conv/6,
      # post game questionnaire 1 to 5
      enjoyment.sMG  = (enjoyment_MG-1)/4,
      responsive.sMG = (responsive_MG-1)/4,
      continue.sMG   = (continue_MG-1)/4
    ) %>%
    mutate(
      soc.exp.MG     = rowMeans(select(., ends_with(".sMG"))),
      soc.exp.CT     = rowMeans(select(., ends_with(".sCT")))
    )
  
  # merge together
  df.indi = merge(df.indi, df.exp) %>% relocate(subID, dyad)
  
  # read in info on passing pseudosync test for OF
  df.key.pseudo = read_csv(file.path(dt.path, "preprocessedOF", 
                                     "df_psync_key_OF.csv")) %>%
    filter(relevant == "credible")
  ls.keys = unique(df.key.pseudo$key)
  
  # add information on the dyads
  df.dyad = df.indi %>%
    group_by(dyad) %>%
    mutate(subID = row_number(), 
           RAADS.mean = mean(RAADS_final), IRI.mean = mean(IRI_final), 
           RAADS.diff = abs(diff(RAADS_final)), IRI.diff = abs(diff(IRI_final))) %>%
    select(dyad, order, subID, gender, age, starts_with("RAADS."), starts_with("IRI.")) %>%
    pivot_wider(names_from = subID, values_from = c(gender, age)) %>%
    mutate(
      dyad.gender = case_when(
        gender_1 == "female" & gender_2 == "female" ~ "female",
        gender_1 == "male"   & gender_2 == "male" ~ "male",
        T ~ "mixed"
      ),
      dyad.age = abs(age_1 - age_2)
    ) %>%
    # merge with which dyad's synchrony exceeded pseudo
    merge(., 
          read_csv(file.path(dt.path, "preprocessedMEA", "df_psync_sub_MEA.csv")) %>%
            select(dyad, relevant))
  
  ## EYE-TRACKING DATA
  df.indi.mg.et = read_csv(file.path(dt.path, "MSYNC_ET_MG_indi.csv")) %>%
    pivot_wider(names_from = ROI, values_from = c(dur.ms, dwell)) %>%
    mutate(
      dyad = substr(subID, 1, 8)
    ) %>% merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(subID, dyad, order) %>%
    mutate_if(is.character, as.factor)
  df.indi.ct.et = read_csv(file.path(dt.path, "MSYNC_ET_CT_indi.csv")) %>%
    pivot_wider(names_from = ROI, values_from = c(dur.ms, dwell)) %>%
    mutate(
      dyad = substr(subID, 1, 8)
    ) %>% merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(subID, dyad, order) %>%
    mutate_if(is.character, as.factor)
  df.dyad.mg.et = read_csv(file.path(dt.path, "MSYNC_ET_MG_dyad.csv")) %>%
    select(-comment) %>%
    merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(dyad, order) %>%
    mutate_if(is.character, as.factor)
  df.dyad.ct.et = read_csv(file.path(dt.path, "MSYNC_ET_CT_dyad.csv")) %>%
    select(-comment) %>%
    merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(dyad, order) %>%
    mutate_if(is.character, as.factor)
  
  ## MEA DATA
  df.indi.ct.mea = read_csv(file.path(dt.path, "MSYNC_mea_CT.csv")) %>%
    filter(position != "B") %>%
    rename("MEA.total.mov" = "total.mv", "MEA.head.mov" = "head.mv",
           "MEA.body.mov"  = "body.mv") %>%
    pivot_wider(names_from = measure, values_from = MEA.sync, names_prefix = "MEA.") %>%
    mutate(
      subID = paste0(dyad, "_", position),
      position = NULL, phase = NULL
    ) %>% merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(subID, dyad, order) %>%
    mutate_if(is.character, as.factor)
  df.dyad.ct.mea = read_csv(file.path(dt.path, "MSYNC_mea_CT.csv")) %>%
    filter(position == "B") %>%
    rename("MEA.total.mov" = "total.mv", "MEA.head.mov" = "head.mv",
           "MEA.body.mov"  = "body.mv") %>%
    pivot_wider(names_from = measure, values_from = MEA.sync, names_prefix = "MEA.") %>%
    mutate(
      position = NULL, phase = NULL
    ) %>%
    merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(dyad, order) %>%
    mutate_if(is.character, as.factor)
  df.indi.mg.mea = read_csv(file.path(dt.path, "MSYNC_mea_MG.csv")) %>%
    filter(position != "B") %>%
    rename("MEA.total.mov" = "total.mv") %>%
    pivot_wider(names_from = measure, values_from = MEA.sync, names_prefix = "MEA.") %>%
    mutate(
      subID = paste0(dyad, "_", position),
      phase = sprintf("p%i", phase),
      position = NULL
    ) %>% merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(subID, dyad, order) %>%
    mutate_if(is.character, as.factor)
  df.dyad.mg.mea = read_csv(file.path(dt.path, "MSYNC_mea_MG.csv")) %>%
    filter(position == "B") %>%
    rename("MEA.total.mov" = "total.mv") %>%
    pivot_wider(names_from = measure, values_from = MEA.sync, names_prefix = "MEA.") %>%
    mutate(
      position = NULL, 
      phase = sprintf("p%i", phase)
    ) %>%
    merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(dyad, order) %>%
    mutate_if(is.character, as.factor)
  
  ## OPENPOSE [!MISSING: only some keys?]
  df.dyad.mg.op = read_csv(file.path(dt.path, "MSYNC_OP_MG.csv")) %>%
    filter(position == "B") %>%
    rename("OP.total.mov" = "QNTmov") %>% 
    pivot_wider(names_from = measure, values_from = IPSmov, names_prefix = "OP.") %>%
    mutate(
      phase = sprintf('p%i', phase),
      position = NULL
    ) %>%
    merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(dyad, order) %>%
    mutate_if(is.character, as.factor)
  df.indi.mg.op = read_csv(file.path(dt.path, "MSYNC_OP_MG.csv")) %>%
    filter(position != "B") %>%
    rename("OP.total.mov" = "QNTmov") %>% 
    pivot_wider(names_from = measure, values_from = IPSmov, names_prefix = "OP.") %>%
    mutate(
      phase = sprintf('p%i', phase),
      subID = paste0(dyad, "_", position),
      position = NULL
    ) %>% merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(subID, dyad, order) %>%
    mutate_if(is.character, as.factor)
  
  ## OPENFACE > only keys which exceeded pseudosynchrony
  df.dyad.ct.of = readRDS(file.path(dt.path, "MSYNC_AU_sync_CT.rds")) %>%
    ungroup() %>% filter(position == "B" & key %in% ls.keys) %>% 
    select(-position, -phase) %>%
    pivot_wider(names_from = measure, values_from = OF.sync, names_prefix = "OF.") %>%
    merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(dyad, order) %>%
    mutate_if(is.character, as.factor)
  df.indi.ct.of = readRDS(file.path(dt.path, "MSYNC_AU_sync_CT.rds")) %>%
    ungroup() %>% filter(position != "B" & key %in% ls.keys) %>% 
    mutate(
      subID = paste0(dyad, "_", position), 
      position = NULL, phase = NULL
    ) %>%
    pivot_wider(names_from = measure, values_from = OF.sync, names_prefix = "OF.") %>%
    merge(., 
          readRDS(file.path(dt.path, "MSYNC_AU_intensity_CT.rds")) %>%
            mutate(
              subID = gsub("_F", "_", ID)
            ) %>% ungroup() %>% select(-ID, -speaker)) %>%
    rename("OF.exp" = "exp") %>% merge(., df.dyad %>% select(dyad, order)) %>%
    relocate(subID, dyad, order) %>%
    mutate_if(is.character, as.factor)
  
  ls.vars = ls()
  
  # save it all
  save(list = ls.vars[grepl("^df.dyad.*|^df.indi.*", ls.vars)],
       file = "MSYNC_data.RData")
  
}