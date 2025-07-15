# (C) Irene Sophia Plank (10planki@gmail.com)

# This script takes all the scripts produced in this pipeline and merges them to
# two csv files, one containing all values of the individual speakers and one
# containing all the values of the dyad. Additionally, it allows for the 
# exclusion of speakers' turn-based information if the categorisation into 
# silences and sounding instances by the uhm-o-meter was suboptimal. The script
# also computes the silence-to-turn ratio of the dyad based on the phonation
# time and the duration of the silence. It also summarises the turn-taking gaps
# by taking the median of all turn-taking gaps for the speaker and the dyad. 

library(tidyverse)

dt.path = c("/Users/vilya/Documents/MSYNC/data/preprocessedAudio", 
            "/Users/vilya/Documents/MSYNC/data")

# Load data ---------------------------------------------------------------

# synchronisation data - information per dyad
df.dyad = read_csv(file.path(dt.path[1], "OUT_sync-dyad.csv"))

df.dyad[c("dyad1", "dyad2", "task")] = str_split_fixed(df.dyad$name, "_", 3)
df.dyad$dyad = paste(df.dyad$dyad1, df.dyad$dyad2, sep = "_")
df.dyad = df.dyad %>% select(-c(name, dyad1, dyad2))

# Pitch intensity and prosodic information --------------------------------

# pitch-intensity and prosodic - information per individual
df.pint = read_csv(file.path(dt.path[1], "OUT_pitch_intensity.csv"))
df.pros = list.files(path = dt.path[1], pattern = "*_prosodic_extraction.csv", full.names = T) %>%
  map_df(~read_csv(., show_col_types = F))
colnames(df.pros) = c("soundname", "nsyll", "npause", "dur", "pho", "spr",
                      "art", "asd", "prepro")

# merge pitch-intensity and prosodic information
df.indi = merge(df.pint, df.pros, all.x = T)
rm(df.pint, df.pros)

# split filenames
df.indi = df.indi %>%
  separate(col = soundname, into = c("dyad1", "dyad2", "task", "prep", "ch", "side"), remove = F) %>%
  mutate(
    dyad = paste0(dyad1, "_", dyad2)
  ) %>% select(-dyad1, -dyad2, -ch, -prep)

# calculate speech rate for dyad
df.spr = df.indi %>% 
  group_by(dyad, task) %>% 
  summarise(
    nsyll = sum(nsyll, na.rm = T), 
    dur   = mean(dur)
  ) %>%
  distinct() %>%
  mutate(
    spr = nsyll / dur
  ) %>%
  select(dyad, task, spr)

# merge speech rate with dyad data frame
df.dyad = merge(df.dyad, df.spr)
rm(df.spr)

# Individual synchronisation ----------------------------------------------

# load data
df.sync = read_csv(file.path(dt.path[1], "OUT_sync-indi.csv")) %>%
  separate(col = name, into = c("dyad1", "dyad2", "task"), remove = F) %>%
  mutate(
    dyad = paste0(dyad1, "_", dyad2),
    side = speaker
  ) %>% select(-dyad1, -dyad2)

# sync data frames
df.indi = merge(df.indi,df.sync, all.x = T)
rm(df.sync)

# Silence-to-turn ratio ---------------------------------------------------

str.dyad = df.indi %>%
  group_by(dyad, task) %>%
  summarise(
    dur = mean(dur),
    pho = sum(pho)
  ) %>% 
  mutate(
    sil  = dur - pho) %>%  # silence time
  mutate(
    str  = sil/pho         # pause to turn ratio
  )

df.dyad = merge(df.dyad, str.dyad)
rm(str.dyad)

# Turn-taking gaps --------------------------------------------------------

# load data
df.turn = read_csv(file.path(dt.path[1], "OUT_turns.csv")) %>%
  mutate(
    name = paste0(dyad, "_", task)
  )

# summarise turn-taking gaps and number of turns
ttg.indi = df.turn %>%
  group_by(dyad, task, speaker) %>%
  summarise(
    ttg = median(ttg, na.rm = T)
  )
df.indi = merge(df.indi, ttg.indi, all.x = T)

ttg.dyad = df.turn %>%
  group_by(dyad, task) %>%
  summarise(
    ttg      = median(ttg, na.rm = T),
    no_turns = n()
  )
df.dyad = merge(df.dyad, ttg.dyad, all.x = T)

rm(ttg.indi, ttg.dyad, df.turn)


# Save --------------------------------------------------------------------

# save it all
write_csv(df.dyad, file.path(dt.path[1], 'OUT_dyad.csv'))
write_csv(df.indi, file.path(dt.path[1], 'OUT_indi.csv'))

