# MSYNC

This repository contains all scripts used to preprocess and analyse the data from the MSYNC project. This project was conducted at the LMU Klinikum in Munich, Germany. Project members are Yasemin Abra, Christine M. Falter-Wagner and Irene Sophia Plank. 

In this project, participants completed two tasks: 

* Mirror game (MG): participants were asked to mirror each others behaviour, without talking
* Conversation task (CT): participants had a ten-minute conversation about foods and drinks they both dislike

In both tasks, we collected video recordings and mobile eye-tracking data. During the conversation task, we additionally recorded speech. After each task, we asked participants to report their subjective experience using questionnaires. The order of tasks was manipulated such that half of the dyads first had the conversation and the other half first completed the mirror game. 

We share anonymised, preprocessed data in `MSYNC_data.RData` as well as all scripts that were used in the preprocessing and analysis of the data. Most scripts are located in the `helpers` folder with the exception of RMarkdown scripts combining code with documentation. These scripts should all run through if the whole repository was downloaded: 

* `MSYNC_checkmodel.Rmd` contains code to perform simulation-based calibration of the models used to analyse the data
* `MSYNC_pseudosync.Rmd` contains code to create pseudosynchrony values and compare these values to real synchrony values
* `MSYNC_supps.Rmd` contains the main analysis code and results. 

`MSYNC_data.RData` contains the following dataframes which can be loaded into R: 

* `df.demo`: group comparisons between the two orders (MG first or CT first)
* `df.dyad`: information on the dyad composition, columns:
	* `dyad`: dyad identifier
	* `order`: order of tasks (MG-CT or CT-MG)
	* `RAADS.mean` / `IRI.mean`: average RAADS / IRI value of the two dyad patners
	* `RAADS.diff` / `IRI.diff`: absolute difference in RAADS / IRI values of the two dyad partners
	* `gender_1` / `gender_2` / `age_1` / `age_2`: gender and age of the two interaction partners, order is random
	* `dyad.gender`: dyad gender composition (mixed, female, male)
	* `dyad.age`: absolute difference in age of the two dyad partners
	* `relevant`: whether or not the pseudosynchrony exceeded MEA synchrony in the CT (credible or not credible)
* `df.dyad.ct.et` / `df.dyad.mg.et`: eye-tracking data on the level of the dyad, columns: 
	* `dyad`: dyad identifier
	* `order`: order of tasks (MG-CT or CT-MG)
	* `shared.face`: shared face attention (ms)
	* `FI`: failed eye contact initiations
	* `SI`: successful eye contact initiations
	* `TI`: total eye contact initiations
	* `SI.ratio`: `SI` divided by `TI`
	* `SI.FI`: `SI` divided by `FI`
	* `phase`: MG phase (only for MG)
* `df.dyad.ct.mea` / `df.dyad.mg.mea`: synchrony data based on MEA on the level of the dyad, columns: 
	* `dyad`: dyad identifier
	* `order`: order of tasks (MG-CT or CT-MG)
	* `phase`: MG phase (only for MG)
	* `MEA.total.mov`: percent of frames on which someone moved
	* `MEA.head.mov`: percent of frames on which someone moved their head (only for CT)
	* `MEA.body.mov`: percent of frames on which someone moved their body (only for CT)
	* `MEA.mean`: grandaverage synchrony based on MEA values
	* `MEA.peak`: average of peak synchrony based on MEA values
	* `MEA.plag`: average lag of the peak (s) - how many seconds away from perfect synchrony was the peak?
	* `MEA.zero`: average synchrony at lag 0, perfect synchrony
* `df.dyad.ct.of`: synchrony data based on OpenFace on the level of the dyad, columns: 
	* `dyad`: dyad identifier
	* `order`: order of tasks (MG-CT or CT-MG)
	* `key`: action unit based on the FACS (AU06_r or AU12_r)
	* `OF.mean`: grandaverage synchrony 
	* `OF.peak`: average of peak synchrony
	* `OF.plag`: average lag of the peak (s) - how many seconds away from perfect synchrony was the peak?
	* `OF.zero`: average synchrony at lag 0, perfect synchrony
* `df.dyad.mg.op`: synchrony data based on OpenPose on the level of the dyad, columns: 
	* `dyad`: dyad identifier
	* `order`: order of tasks (MG-CT or CT-MG)
	* `phase`: MG phase (only for MG)
	* `OP.total.mov`: total movement of all pose
	* `axis`: axis of the pose coordinate (x or y)
	* `key`: which pose coordinate, e.g., hip, elbowL (left elbow) etc.
	* `OP.mean`: grandaverage synchrony  
	* `OP.peak`: average of peak synchrony
	* `OP.plag`: average lag of the peak (s) - how many seconds away from perfect synchrony was the peak?
	* `OP.zero`: average synchrony at lag 0, perfect synchrony

The naming scheme for the individual data is the same, but while all data in dyad dataframes contains interpersonal synchrony, data in the individual dataframes contains interpersonal adaptation; thus, an inidivdual feature and not a dyadic feature. 
The individual data frames contain the following additional columns: 

* all dataframes: `subID` as a subject identifier
* `df.indi.ct.et` / `df.indi.mg.et`
	* `dwell.face`: dwell time to the face of the interaction partner, regardless of where they look
* `df.indi`: subject information, including the following additional columns:
	* `enjoyment_MG`, `responsive_MG`, `continue_MG`: post-game questionnaire (1 to 5)
	* `IOS_MG` / `IOS_CT`: inclusion of the self in the other rating (1 to 7)
	* `mood_MG` / `close_MG` / `similar_MG` / `connect_MG` / `trust_MG` / `mood_conv` / `close_conv` / `similar_conv` / `connect_conv` / `trust_conv`: Mood Cohesion Trust questionnaire for MG and CT (conv; 1 to 10)
	* `likeable_MG` / `friendly_MG` / `attentive_MG` / `smooth_MG` / `comfort_MG` / `likeable_conv` / `friendly_conv` / `attentive_conv` / `smooth_conv` / `comfort_conv`: rapport questionnaire (0 to 6) 
	* `role`: role in the MG (either `Leader1` = leader in the first phase or `Leader2` = leader in the second phase)
	* `age`: age in years
	* `gender`: gender description
	* `gender_identity`: whether gender matches AGAB (`Ja` = cis, `Nein` = trans)
	* `edu`: education level (1 = none to 5 = university degree)
	* `vision`: whether no vision impairment or how corrected to normal (none, glasses or contacts)
	* `handedness`: self-reported handedness (`link` = left-handed, `rechts` = right-handed, none chose ambidextrous)
	* `RAADS_final` / `IRI_final` / `ADC_final` / `ECR` : questionnaire scores
	* `CFT_iq`: intelligence estimate based on CFT
	* all columns ending in `.sMG` and `.sCT` are rescaled experience ratings (0 to 1) from above which are then averaged to social experience scores (`soc.exp.MG` and `soc.exp.CT`)

