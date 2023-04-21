###################################################################
#                                                                 #
# Documentation and vignettes for the 'EHR' and 'medExtractR'     #
# packages can be found at:                                       #
#                                                                 #
# https://cran.r-project.org/web/packages/EHR/index.html          #
#                                                                 #
# https://cran.r-project.org/web/packages/medExtractR/index.html  #
#                                                                 #
###################################################################

# After installing 'EHR' and 'medExtractR' packages, load into environment.
#install.packages("EHR")
#install.packages("medExtractR")
library(EHR)
library(medExtractR)

##### Extract-Med Module #####

# Load in example notes included in 'EHR' package
tac_notes <- list(
  system.file("examples", "tacpid1_2008-06-26_note1_1.txt", package = "EHR"),
  system.file("examples", "tacpid1_2008-06-26_note2_1.txt", package = "EHR"),
  system.file("examples", "tacpid1_2008-12-16_note3_1.txt", package = "EHR")
)

# Execute Extract-Med Module
tac_mxr <- extractMed(note_fn = tac_notes, #File name or vector/list of notes
                      drugnames = c("tacrolimus", "prograf", #Names of drugs we want to extract
                                    "tac", "tacro", 
                                    "fk", "fk506"), 
                      drgunit = "mg", #Unit of drug
                      windowlength = 60, #Length of search window after drug name
                      max_edit_dist = 2, #Define edit distance to capture misspellings in drug names (e.g tacrolimius)
                      strength_sep = NULL, #Special character(s) to separate doses administered multiple times a day
                      lastdose = TRUE) #Input for whether last dose time entity should be extracted

# Show first six rows of tac_mxr dataframe
head(tac_mxr)

# Output of Extract-Med module is a dataframe containing four columns:
# (1) filename: name of clinical note file corresponding to results
# (2) entity: label of entity for extracted expression (e.g. DrugName, Strength, DoseAmt)
# (3) expr: expression extracted from clinical note (e.g. Prograf, 1 mg, 3)
# (4) pos: position of the extracted expression within the note, formatted as start:stop

# Export to csv
write.csv(tac_mxr, file="tac_mxr.csv", row.names=FALSE)

##### Tuning medExtractR #####

# Read in the annotations - might be specific to annotation method/software
ann_filenames <- list(system.file("mxr_tune", "tune_note1.ann", package = "EHR"),
                      system.file("mxr_tune", "tune_note2.ann", package = "EHR"),
                      system.file("mxr_tune", "tune_note3.ann", package = "EHR"))

# Format annotations 
tune_ann <- do.call(rbind, lapply(ann_filenames, function(fn){
  annotations <- read.delim(fn, 
                            header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                            col.names = c("id", "entity", "annotation"))
  
  # Label with file name
  annotations$filename <- sub(".ann", ".txt", sub(".+/", "", fn), fixed=TRUE)
  
  # Separate entity information into entity label and start:stop position
  # Format is "entity start stop"
  ent_info <- strsplit(as.character(annotations$entity), split="\\s")
  annotations$entity <- unlist(lapply(ent_info, '[[', 1))
  annotations$pos <- paste(lapply(ent_info, '[[', 2), 
                           lapply(ent_info, '[[', 3), sep=":")
  
  annotations <- annotations[,c("filename", "entity", "annotation", "pos")]
  
  return(annotations)
}))

head(tune_ann)

# NOTE: THIS EXAMPLE ONLY USES 3 NOTES TO TUNE PARAMETERS. RECOMMEND USING >=10 NOTES THAT CONTAIN ENTITIES THAT WILL BE EXTRACTED WHEN APPLIED TO ALL NOTES. THESE TRAINING NOTES SHOULD ALSO BE MANUALLY ANNOTATED BY REVIEWERS TO IDENTIFY RELEVANT ENTITIES TO PRODUCE A GOLD STANDARD OF ANNOTATIONS (SEE BRAT RAPID ANNOTATION TOOL (BRAT) SOFTWARE -- "Stenetorp P, Pyysalo S, Topic G, Ohta T, Ananiadou S, Tsujii JI. BRAT: a web-based tool for NLP-assisted text annotation. InProceedings of the Demonstrations at the 13th Conference of the European Chapter of the Association for Computational Linguistics 2012 Apr 23 (pp. 102-107). Association for Computational Linguistics.")

# Define ranges of possible values for windowlength and max_edit_dist
wind_len <- seq(30, 120, 30)
max_edit <- seq(0, 2, 1)
tune_pick <- expand.grid("window_length" = wind_len, 
                         "max_edit_distance" = max_edit)

# Run the Extract-Med module on the tuning notes
note_filenames <- list(system.file("mxr_tune", "tune_note1.txt", package = "EHR"),
                       system.file("mxr_tune", "tune_note2.txt", package = "EHR"),
                       system.file("mxr_tune", "tune_note3.txt", package = "EHR"))

# List to store output for each parameter combination
mxr_tune <- vector(mode="list", length=nrow(tune_pick))

for(i in 1:nrow(tune_pick)){
  mxr_tune[[i]] <- extractMed(note_filenames,
                              drugnames = c("tacrolimus", "prograf", 
                                            "tac", "tacro", 
                                            "fk", "fk506"),
                              drgunit = "mg",
                              windowlength = tune_pick$window_length[i],
                              max_edit_dist = tune_pick$max_edit_distance[i],
                              progress = FALSE)
}

# Functions to compute true positive, false positive, and false negatives
# number of true positives - how many annotations were correctly identified by extractMed
Tpos <- function(df){
  sum(df$annotation == df$expr, na.rm=TRUE)
}
# number of false positive (identified by extractMed but not annotated)
Fpos <- function(df){
  sum(is.na(df$annotation))
}
# number of false negatives (annotated but not identified by extractMed)
Fneg <- function(df){
  df_ann <- subset(df, !is.na(annotation)) #keep only rows with annotation
  sum(is.na(df$expr))
}

prf <- function(df){
  tp <- Tpos(df)
  fp <- Fpos(df)
  fn <- Fneg(df)
  
  precision <- tp/(tp + fp)
  recall <- tp/(tp + fn) 
  f1 <- (2*precision*recall)/(precision + recall)
  
  return(f1)
}

tune_pick$F1 <- sapply(mxr_tune, function(x){
  compare <- merge(x, tune_ann, 
                   by = c("filename", "entity", "pos"), all = TRUE)
  prf(compare)
})

# Plot F1-scores over defined grid
library(ggplot2)

ggplot(tune_pick) + geom_point(aes(max_edit_distance, window_length, size = F1)) + 
  scale_y_continuous(breaks=seq(30,120,30)) + 
  annotate("text", x = tune_pick$max_edit_distance+.2, y = tune_pick$window_length,
           label = round(tune_pick$F1, 2)) + 
  ggtitle("F1 for tuning parameter values")

##### Pro-Med-NLP Module: Generate longitudinal medication dose data from the output of MedExtractR #####

## Part I: Parse raw output and pairs entities together

# Load in dataset containing the same version of our output from the Extract-Med module into your global environment. Then, run the parseMedExtractR() function on your file
tac_mxr2 <- system.file("examples", "tac_mxr.csv", package = "EHR") #Load in dataset

tac_parsed <- parseMedExtractR(tac_mxr2) #Parse function. NOTE: you can only read in a file, not a dataframe

# Show first six rows of tac_parsed dataframe
head(tac_parsed)

# Output of parse function is a dataframe that transposes entities to column headers with corresponding expressions and their start:stop, with each row representing a separate mention of the drug of interest

# Run build dose function to pair entities together
tac_bd <- buildDose(dat = tac_parsed, #Output of parse function
                    checkForRare = TRUE) #Argument that, if set to TRUE, shows extracted expression with a proportion of occurrences less than 20% for each entity. Could be used to detect putatively inaccurate data before completely running Pro-Med-NLP Module. However, these values could still be correct extractions.

# There are additional arguments that can be included in the build dose function, including "dn", preserve", "dist_method", but these were not included here as they are either (a) not needed when processing medExtractR outputs or (b) the default is adequate. See 'EHR' documentation for more information on these arguments.

# Show first six rows of tac_bd dataframe
head(tac_bd)

## Part II: Form final datasets containing dosing information at the note and date level for each patient

# Create function for extracting meta data (e.g. pid/patient id, date) from each note
bmd <- function(x) {
  fns <- strsplit(x, '_')
  pid <- sapply(fns, `[`, 1)
  date <- as.Date(sapply(fns, `[`, 2), format = '%Y-%m-%d')
  note <- sapply(fns, `[`, 3)
  data.frame(filename = x, pid, date, note, stringsAsFactors = FALSE)
}

# Extract meta data from tac_bd filename column
tac_metadata <- bmd(tac_bd[['filename']])
head(tac_metadata)

# Handling 'lastdose'
data(tac_lab, package = 'EHR') #load in lab dataset containing times of lab drug measurements
tac_lab

tac_ld <- processLastDose(mxrData = tac_mxr, #output from extractMed function
                          noteMetaData = tac_metadata, #note meta data
                          labData = tac_lab) #dataframe containing lab dates and times associated with file names. Must contain patient ID (pid), date (in the same format as noteMetaData), and lab time

tac_ld #show dataframe

# Results of processLastDose() function converts time expressions (e.g. 10am, 7 last night, 14 hour level) into a processed and standardized time variable. 

tac_out_ld <- addLastDose(buildData = tac_bd, #output from buildDose()
                          lastdoseData = tac_ld) #dataset containing last dose times with corresponding file names in buildData()

# Since multiple times could be extracted from a clinical note, extracted times within 2 hours of each other are treated as equivalent. This rule was determined based on drugs administered every 12 hours and assuming a trough concentration drug. For time differences >2 hours, the last dose start position is used to pair the extracted time with the closest drug mention. You could also provide your own validated dataset with last dose times.

# Collapse data to date and note level to remove re-dundencies in data
tac_cD <- collapseDose(x = tac_out_ld, #Output from buildDose function or from addLastDose if last dose info included 
                       noteMetaData = tac_metadata, #Dataframe with columns filenmae, pid, date, and note
                       naFreq = 'most') #Method to use when assigning missing freq (default is the most common freq)

# If the drug you are processing has different formulations (e.g. extended release) that you want to separately take into consideration, they can be separated using a regular expression (e.g. "xr|er") after the naFreq argument in collapseDose(). Additional information on this can be found in the vignette's for the 'EHR' package (see Lamotrigine examples).

# Show first six rows of note level collapsing
head(tac_cD$note,10)

# Show first six rows of date level collapsing
head(tac_cD$date,10)

# Collapsing by date or note produces observations at the daily intake level. It is possible to collapse date even further to the daily level. You may want to drop the 'dose.intake' variable from your dataframe as it may potentially lose meaning.

# Collapse data to daily level using notes
x <- tac_cD[['note']] #create new dataframe with collapsed data at note level
k1 <- tac_metadata[match(x[,'filename'], tac_metadata[,'filename']), c('pid','date','note')] #retrieve metadata for each filename
k2 <- cbind(k1, x[,c('dose.daily','drugname_start')]) #select additional key data
chk <- do.call(paste, c(k2, sep = '|')) #turn keys into character string
tac_cD_note_dailyintake <- x[!duplicated(chk),] #keep first instance of each chk key
head(tac_cD_note_dailyintake[,c('filename','drugname','drugname_start','dose.daily')], 15) #show first 15 rows

# Collapse data to daily level using dates
x <- tac_cD[['date']] 
k1 <- tac_metadata[match(x[,'filename'], tac_metadata[,'filename']), c('pid','date')] #retrieve metadata for each filename, ignoring note for date level collapsing
k2 <- cbind(k1, x[,c('dose.daily','drugname_start')]) #select additional key data
chk <- do.call(paste, c(k2, sep = '|')) #turn keys into character string
tac_cD_date_dailyintake <- x[!duplicated(chk),] #keep first instance of each chk key
head(tac_cD_date_dailyintake[,c('filename','drugname','drugname_start','dose.daily')], 15) #show first 15 rows

##### Dose Building Using Example Vanderbilt EHR Data #####

## Part I: Parse raw output and pairs entities together

# Load in new dataset containing a larger version of the output from the Extract-Med module into your global environment. Then, run the parseMedExtractR() function on your file
tac_mxr_VU <- system.file("examples", "tac_mxr_out.csv", package = "EHR") #Load new dataset
tac_mxr_VU_notes <- read.csv(tac_mxr_VU, na = '')
head(tac_mxr_VU_notes)

tac_parsed_VU <- parseMedExtractR(tac_mxr_VU) #Parse function. NOTE: you can only read in a file, not a dataframe

# Show first six rows of tac_parsed dataframe
head(tac_parsed_VU)

# Output of parse function is a dataframe that transposes entities to column headers with corresponding expressions and their start:stop, with each row representing a separate mention of the drug of interest

# Run build dose function to pair entities together
tac_bd_VU <- buildDose(dat = tac_parsed_VU, #Output of parse function
                       checkForRare = TRUE) #Argument that, if set to TRUE, shows extracted expression with a proportion of occurrences less than 20% for each entity. Could be used to detect putatively inaccurate data before completely running Pro-Med-NLP Module. However, these values could still be correct extractions.

# There are additional arguments that can be included in the build dose function, including "dn", preserve", "dist_method", but these were not included here as they are either (a) not needed when processing medExtractR outputs or (b) the default is adequate. See 'EHR' documentation for more information on these arguments.

# Show first six rows of tac_bd dataframe
head(tac_bd_VU)

# Compare to Gold Standard
tac_gs_part1 <- read.csv(system.file("examples", "tac_gs_part1.csv", package = "EHR"),
                         stringsAsFactors = FALSE, na = '')

precall <- function(dat, gs) { #create function to provide precision and recall measures
  tp1 <- sum(dat %in% gs)
  fp1 <- sum(!(dat %in% gs))
  fn1 <- sum(!(gs %in% dat))
  r1 <- c(tp1, tp1 + fn1)
  p1 <- c(tp1, tp1 + fp1)
  r <- rbind(r1,p1)
  dimnames(r) <- list(c('recall','prec'), c('num','den'))
  cbind(r, prop = round(r[,1] / r[,2], 2))
}

colsToCompare <- c('filename','drugname','strength','dose','route','freq',
                   'dosestr','dosechange','drugname_start')
tac_bd_VU <- tac_bd_VU[,colsToCompare]
tac_gs_part1 <- tac_gs_part1[,colsToCompare]

tacxrrow <- do.call(paste, c(tac_bd_VU, sep = '|'))
gs.tacxrrow <- do.call(paste, c(tac_gs_part1, sep = '|'))

precall(tacxrrow, gs.tacxrrow)

## Part II: Form final datasets containing dosing information at the note and date level for each patient

# Create function for extracting meta data (e.g. pid/patient id, date) from each note
bmd <- function(x) {
  fns <- strsplit(x, '_')
  pid <- sapply(fns, `[`, 1)
  date <- as.Date(sapply(fns, `[`, 2), format = '%Y-%m-%d')
  note <- sapply(fns, `[`, 3)
  data.frame(filename = x, pid, date, note, stringsAsFactors = FALSE)
}

# Extract meta data from tac_bd filename column
tac_metadata_VU <- bmd(tac_bd_VU[['filename']])

# Collapse data to date and note level to remove redundencies in data
tac_cD_VU <- collapseDose(x = tac_bd_VU, #Output from buildDose function or from addLastDose if last dose info included 
                          noteMetaData = tac_metadata_VU, #Dataframe with columns filenmae, pid, date, and note
                          naFreq = 'most') #Method to use when assigning missing freq (default is the most common freq)

# If the drug you are processing has different formulations (e.g. extended release) that you want to separately take into consideration, they can be separated using a regular expression (e.g. "xr|er") after the naFreq argument in collapseDose(). Additional information on this can be found in the vignette's for the 'EHR' package (see Lamotrigine examples).

# Show first six rows of note level collapsing
head(tac_cD_VU$note,10)

# Show first six rows of date level collapsing
head(tac_cD_VU$date,10)

# Collapsing by date or note produces observations at the daily intake level. It is possible to collapse date even further to the daily level. You may want to drop the 'dose.intake' variable from your dataframe as it may potentially lose meaning.

# Collapse data to daily level using notes
x <- tac_cD_VU[['note']] #create new dataframe with collaposed data at note level
k1 <- tac_metadata_VU[match(x[,'filename'], tac_metadata_VU[,'filename']), c('pid','date','note')] #retrieve metadata for each filename
k2 <- cbind(k1, x[,c('dose.daily','drugname_start')]) #select additional key data
chk <- do.call(paste, c(k2, sep = '|')) #turn keys into character string
tac_cD_note_dailyintake_VU <- x[!duplicated(chk),] #keep first instance of each chk key
head(tac_cD_note_dailyintake_VU[,c('filename','drugname','drugname_start','dose.daily')], 15) #show first 15 rows

# Collapose data to daily level using dates
x <- tac_cD_VU[['date']] 
k1 <- tac_metadata_VU[match(x[,'filename'], tac_metadata_VU[,'filename']), c('pid','date')] #retrieve metadata for each filename, ignoring note for date level collapsing
k2 <- cbind(k1, x[,c('dose.daily','drugname_start')]) #select additional key data
chk <- do.call(paste, c(k2, sep = '|')) #turn keys into character string
tac_cD_date_dailyintake_VU <- x[!duplicated(chk),] #keep first instance of each chk key
head(tac_cD_date_dailyintake_VU[,c('filename','drugname','drugname_start','dose.daily')], 15) #show first 15 rows

# Compare to Gold Standard

tac_gs_part2_note <- read.csv(system.file("examples", "tac_gs_part2_note.csv", package = "EHR"), #Note level GS
                              stringsAsFactors = FALSE, na = '')

tac_gs_part2_date <- read.csv(system.file("examples", "tac_gs_part2_date.csv", package = "EHR"), #Date level GS
                              stringsAsFactors = FALSE, na = '')

precall <- function(dat, gs) {
  tp1 <- sum(dat %in% gs)
  fp1 <- sum(!(dat %in% gs))
  fn1 <- sum(!(gs %in% dat))
  r1 <- c(tp1, tp1 + fn1)
  p1 <- c(tp1, tp1 + fp1)
  r <- rbind(r1,p1)
  dimnames(r) <- list(c('recall','prec'), c('num','den'))
  cbind(r, prop = round(r[,1] / r[,2], 2))
}

metaData <- bmd(unique(tac_gs_part1$filename))
tacxr <- collapseDose(tac_gs_part1, metaData, 'bid')

tacxr.note <- tacxr[['note']]
tacxr.date <- tacxr[['date']]

tacxr.note$pid <- sub("_.*","",tacxr.note$filename)
tacxr.date$pid <- sub("_.*","",tacxr.date$filename)
tac_gs_part2_note$pid <- sub("_.*","",tac_gs_part2_note$filename)
tac_gs_part2_date$pid <- sub("_.*","",tac_gs_part2_date$filename)

tacxrrow.note.intake <- do.call(paste, c(tacxr.note[,c('pid','dose.intake',
                                                       'dosechange')],sep = '|'))
tacxrrow.note.daily <- do.call(paste, c(tacxr.note[,c('pid','intaketime','dose.daily',
                                                      'dosechange')], sep = '|'))
tacxrrow.date.intake <- do.call(paste, c(tacxr.date[,c('pid','dose.intake',
                                                       'dosechange')], sep = '|'))
tacxrrow.date.daily <- do.call(paste, c(tacxr.date[,c('pid','intaketime','dose.daily',
                                                      'dosechange')], sep = '|'))

gs.tacxrrow.note.intake <- do.call(paste, c(tac_gs_part2_note[,c('pid','doseintake',
                                                                 'dosechange')], sep = '|'))
gs.tacxrrow.note.daily <- do.call(paste, c(tac_gs_part2_note[,c('pid','intaketime','daily',
                                                                'dosechange')], sep = '|'))
gs.tacxrrow.date.intake <- do.call(paste, c(tac_gs_part2_date[,c('pid','doseintake',
                                                                 'dosechange')], sep = '|'))
gs.tacxrrow.date.daily <- do.call(paste, c(tac_gs_part2_date[,c('pid','intaketime','daily',
                                                                'dosechange')], sep = '|'))

precall(tacxrrow.note.intake, gs.tacxrrow.note.intake)

precall(tacxrrow.note.daily, gs.tacxrrow.note.daily)

precall(tacxrrow.date.intake, gs.tacxrrow.date.intake)

precall(tacxrrow.date.daily, gs.tacxrrow.date.daily)

##### Pre-Processing for Raw Extracted Data #####
#install.packages(pkdata)
#install.packages(lubridate)
library(pkdata)
library(lubridate)

# This examples shows the three main steps of pre-processing: (1) read and clean the raw data, (2) merge raw data to create new ID variables, and (3) make new data for use with modules. Each raw dataset contains a subject unique ID ('subject_uid') and subject visit ID ('subject_id').

# Define directories
td <- tempdir()
dir.create(file.path(td, 'data2'))
dir.create(file.path(td, 'check2'))
rawDataDir <- system.file("examples", "str_ex2", package="EHR")
dataDir <- file.path(td, 'data2')
checkDir <- file.path(td, 'check2')

# Demographics data
demo.in <- readTransform(file.path(rawDataDir, "Demographics_DATA.csv")) #read in raw data
head(demo.in)

# Concentration sampling times data
# Read in and transform data
samp.in <- readTransform(file.path(rawDataDir, "SampleTimes_DATA.csv"),
                         rename = c('Study.ID' = 'subject_id'), #rename subject ID so that it is consistent across datasets
                         modify = list(samp = expression(as.numeric(sub('Sample ', '', Event.Name))))) #create new variable called 'samp', which indexes the sample number
head(samp.in)

# helper function used to make subject_id
sampId <- function(x) {
  # remove leading zeroes or trailing periods
  subid <- gsub('(^0*|\\.$)', '', x)
  # change _ to .
  gsub('_([0-9]+[_].*)$', '.\\1', subid)
}

# Concentration sample values data
conc.in <- readTransform(file.path(rawDataDir, "SampleConcentration_DATA.csv"), #read in raw data
                         modify = list( #use helper function to make subject_id
                           subid = expression(sampId(name)), 
                           subject_id = expression(as.numeric(sub('[_].*', '', subid))), 
                           samp = expression(sub('[^_]*[_]', '', subid)),
                           name = NULL,
                           data_file = NULL,
                           subid = NULL
                         )
)
head(conc.in)

# FLOW dosing data
flow.in <- readTransform(file.path(rawDataDir, "FLOW_DATA.csv"), #read in raw data
                         rename = c('Subject.Id' = 'subject_id',
                                    'Subject.Uniq.Id' = 'subject_uid')) 
head(flow.in)

# MAR dosing data
mar.in0 <- read.csv(file.path(rawDataDir, "MAR_DATA.csv"), check.names = FALSE) #read in raw data
mar.in <- dataTransformation(mar.in0, rename = c('Uniq.Id' = 'subject_uid'))
head(mar.in)

# Serum creatinine lab data
creat.in <- readTransform(file.path(rawDataDir, "Creatinine_DATA.csv"), #read in raw data
                          rename = c('Subject.uniq' = 'subject_uid'))
head(creat.in)

# Albumin lab data
alb.in <- readTransform(file.path(rawDataDir, "Albumin_DATA.csv"), #read in raw data
                        rename = c('Subject.uniq' = 'subject_uid'))
head(alb.in)

# Merge all ID datasets
data <-  list(demo.in,
              samp.in,
              conc.in,
              flow.in,
              mar.in,
              creat.in,
              alb.in)

idcols <-  list(c('subject_id', 'subject_uid'), #id vars in demo.in
                'subject_id', #id var in samp.in
                'subject_id', #id var in conc.in
                c('subject_id', 'subject_uid'), #id vars in flow.in
                'subject_uid', #id var in mar.in
                'subject_uid', #id var in creat.in
                'subject_uid') #id var in creat.in

# Use idCrosswalk() to merge all cleaned input datasets with new IDs
mod.id <- idCrosswalk(data=data, #accepts list of input datasets
                      idcols=idcols, #accepts a list of vectors or character string that identify the ID var in corresponding input dataset
                      visit.id="subject_id", uniq.id="subject_uid")
saveRDS(mod.id, file=file.path(dataDir,"Fentanyl_module_id.rds"))

# Show output of idCrosswalk(). 'mod_id_visit' uniquely identifies each subjects' visit/course
mod.id

# Make new data for use with modules

# Using the pullFakeID() function we can replace the original IDs in each dataset with new IDs that can be used by the data processing modules

## demographics data
demo.cln <- pullFakeId(dat=demo.in, #cleaned demographics input dataframe 
                       xwalk=mod.id, #crosswalk dataframe produced using idCrosswalk()
                       firstCols = c('mod_id', 'mod_visit', 'mod_id_visit'), #controls which variables are in the first columns of output
                       uniq.id = 'subject_uid') #set which variable is the unique identifier for each subject
head(demo.cln)

saveRDS(demo.cln, file=file.path(dataDir,"Fentanyl_demo_mod_id.rds"))

## Drug level data
# sampling times
samp.cln <- pullFakeId(dat=samp.in, #cleaned samples input dataframe 
                       xwalk=mod.id, #crosswalk dataframe produced using idCrosswalk()
                       firstCols = c('mod_id', 'mod_visit', 'mod_id_visit', 'samp'), #controls which variables are in the first columns of output
                       orderBy = c('mod_id_visit','samp'), #sort order
                       uniq.id = 'subject_uid') #set which variable is the unique identifier for each subject
head(samp.cln)

saveRDS(samp.cln, file=file.path(dataDir,"Fentanyl_samp_mod_id.rds"))

# Sampling concentrations
conc.cln <- pullFakeId(dat=conc.in, #cleaned concentrations input dataframe 
                       xwalk=mod.id, #crosswalk dataframe produced using idCrosswalk()
                       firstCols = c('record_id', 'mod_id', 'mod_visit', 'mod_id_visit', 'samp'), #controls which variables are in the first columns of output
                       orderBy = 'record_id', #sort order
                       uniq.id = 'subject_uid') #set which variable is the unique identifier for each subject
head(conc.cln)

saveRDS(conc.cln, file=file.path(dataDir,"Fentanyl_conc_mod_id.rds"))

## Dosing data
# Flow
flow.cln <- pullFakeId(dat=flow.in, #cleaned flowsheet input dataframe 
                       xwalk=mod.id, #crosswalk dataframe produced using idCrosswalk()
                       firstCols = c('mod_id', 'mod_visit', 'mod_id_visit'), #controls which variables are in the first columns of output
                       uniq.id = 'subject_uid') #set which variable is the unique identifier for each subject
head(flow.cln)

saveRDS(flow.cln, file=file.path(dataDir,"Fentanyl_flow_mod_id.rds"))

# MAR
mar.cln <- pullFakeId(dat=mar.in, #cleaned MAR input dataframe 
                      xwalk=mod.id, #crosswalk dataframe produced using idCrosswalk()
                      firstCols = 'mod_id', #controls which variables are in the first columns of output
                      uniq.id = 'subject_uid') #set which variable is the unique identifier for each subject
head(mar.cln)

saveRDS(mar.cln, file=file.path(dataDir,"Fentanyl_mar_mod_id.rds"))

## laboratory data
creat.cln <- pullFakeId(dat=creat.in, #cleaned creatinine input dataframe 
                        xwalk=mod.id, #crosswalk dataframe produced using idCrosswalk()
                        firstCols='mod_id', #controls which variables are in the first columns of output
                        uniq.id = 'subject_uid') #set which variable is the unique identifier for each subject
head(creat.cln)

alb.cln <- pullFakeId(dat=alb.in, #cleaned albumin input dataframe 
                      xwalk=mod.id, #crosswalk dataframe produced using idCrosswalk()
                      firstCols='mod_id', #controls which variables are in the first columns of output
                      uniq.id = 'subject_uid') #set which variable is the unique identifier for each subject
head(alb.cln)

saveRDS(creat.cln, file=file.path(dataDir,"Fentanyl_creat_mod_id.rds"))
saveRDS(alb.cln, file=file.path(dataDir,"Fentanyl_alb_mod_id.rds"))

# set crosswalk option 
xwalk <- readRDS(file.path(dataDir, "Fentanyl_module_id.rds")) #this file is the same as our previously made mod.id dataframe
options(pkxwalk = 'xwalk') #allows modules to access the crosswalk file

# define parameters
drugname <- 'fent' #drugname stub
LLOQ <- 0.05 #define lower limit of quantification

##### Pro-Demographic #####

# Can be used to exclude patients based on specified criteria and create new variables

# Helper function
exclude_val <- function(x, val=1) 
{ !is.na(x) & x == val }

demo.out <- run_Demo(demo.path = file.path(dataDir, "Fentanyl_demo_mod_id.rds"), #read in fentanyl demographics file
                     toexclude = expression(exclude_val(in_hospital_mortality) | exclude_val(add_ecmo)), #exclude patients with a value of 1 in 'in_hospital_mortality' or 'add_ecmo'
                     demo.mod.list = list(length_of_icu_stay = expression(daysDiff(surgery_date, date_icu_dc)))) #create new variable 'length_of_icu_stay'

# Show first 6 rows of demo.out dataframe
head(demo.out$demo)

# Show mod_id_visit values that were excluded
demo.out$exclude

##### Pro-Med-Str #####

## Part I: IV Dose Data

# IV Dose data at Vanderbilt are located in the Flowsheet and medication administration records (MAR) data. Vanderbilt's flowsheet data contains information on record infusion rates and changes to all infusions for all inpatients outside of the operating room. The MAR data records all bolus doses of medications and infusions administered in the operating room. The function run_MedStrI() will generate several files to check potential data errors and allow for you to provide feedback via corrected 'fix' files that you provide. NOTE: IF CORRECTED FILES ARE PROVIDED, THIS FUNCTION SHOULD BE RE-RUN TO INCORPORATE THE CORRECTIONS. 

ivdose.out <- run_MedStrI(flow.path=file.path(dataDir,"Fentanyl_flow_mod_id.rds"), #file path where flowsheet data exists
                          flow.select = c('mod_id','mod_id_visit','Perform.Date','Final.Wt..kg.', #list of variables in the flowsheet data
                                          'Final.Rate..NFR.units.','Final.Units'),
                          flow.rename = c('mod_id','mod_id_visit', 'Perform.Date', 'weight', #rename variables defined by flow.select
                                          'rate', 'final.units'),
                          flow.mod.list = list( #list containing modifications to variables in the flowsheet data
                            date.time = expression(parse_dates(fixDates(Perform.Date))),
                            unit = expression(sub('.*[ ]', '', rate)),
                            rate = expression(as.numeric(sub('([0-9.]+).*', '\\1', rate)))),
                          medchk.path = file.path(rawDataDir, sprintf('medChecked-%s.csv', drugname)), #file path with list of medication names and variants in MAR data to be used
                          mar.path = file.path(dataDir,"Fentanyl_mar_mod_id.rds"), #file path where MAR data exists
                          demo.list = demo.out, #file name for processed demographics file used to exclude subjects based on exclusion criteria
                          check.path = checkDir, #file path where the generated files for data checking are stored and corresponding data files with fixed data exist
                          failflow_fn = 'FailFlow', #filename stub for invalid duplicate rows in the flowsheet data with rate of 0 check file. This will produce  a .csv file with the invalid data. Since we define this as 'FailFlow', the file 'failFailFlow.csv' with invalid duplicate rows will be created in the directory specific by check.path. The corrected version names 'fixFailFlow.csv' should be placed in the same directory.
                          failunit_fn = 'Unit', #filename stub for records with units other than those specified with infusion.unit and bolus.unit. If there are rows containing units other than those specified then a file will be generated similar to as described above. 
                          failnowgt_fn = 'NoWgt', #filename stub for records with missing weight in the flowsheet data and unit involving 'kg'. If there are rows containing missing weights, a file will be generated similar to as described above.
                          infusion.unit = 'mcg/kg/hr', #string specifying units for infusion doses.
                          bolus.unit = 'mcg', #string specifying units for bolus doses
                          bol.rate.thresh = Inf, #upper bound for retaining bolus doses (i.e. bolus units with a rate above the threshold are dropped). The default is 'Inf" which will keep all bolus doses
                          drugname = drugname) #drug name of interest (we previously defined this as 'fent')

head(ivdose.out)

## Part II: E-prescription Data

# For this module, all prescriptions must be for only one drug, but can handle different names (e.g. lamictal and lamotrigine). Data used in the module must include columns for ID, date, strength, dose amount, and frequency. The tasks this will handle include: (1) creating numeric variables for strength, dose, and frequency; (2) calculating daily dose; and (3) removing duplicate daily doses.

eRX <- read.csv(file.path(rawDataDir,"e-rx_DATA.csv"),stringsAsFactors = FALSE) #read in data
eRX

eRX.out <- run_MedStrII(file.path(rawDataDir,"e-rx_DATA.csv"), #file name of the prescription data
                        select = c('GRID','MED_NAME','RX_DOSE','FREQUENCY','ENTRY_DATE','STRENGTH_AMOUNT','DESCRIPTION'), #names of columns to select
                        rename = c('ID','MED_NAME','RX_DOSE','FREQUENCY','ENTRY_DATE','STRENGTH_AMOUNT','DESCRIPTION')) #rename columns

eRX.out

# In the above output, daily.dose was calculated by multiplying strength*dose*freq.num and a redundant daily dose was removed for the patient ID2. For patient ID3, the strength of 100 from the description was used because STRENGTH_AMOUNT was missing. For patient ID6, the dose amounts of 1.5, 1, and 1.5 are added together to equal a dose of 4 in the calculation of daily dose.

##### Pro-Drug Level #####

# Processes drug concentration data that can be merged with other datasets. Similar to above, run_DrugLevel() will generate several files to check potential data errors and allow for you to provide feedback via corrected 'fix' files that you provide.

conc.out <- run_DrugLevel(conc.path=file.path(dataDir,"Fentanyl_conc_mod_id.rds"), #file path of concentration data exists
                          conc.select=c('mod_id','mod_id_visit','samp','fentanyl_calc_conc'), #list of variables in drug conc data
                          conc.rename=c(fentanyl_calc_conc = 'conc.level', samp= 'event'), #rename variables in conc.select
                          conc.mod.list=list(mod_id_event = expression(paste(mod_id_visit, event, sep = '_'))), #list containing modifications to variables in the drug concentration data
                          samp.path=file.path(dataDir,"Fentanyl_samp_mod_id.rds"), #file path where sampling time data exists
                          samp.mod.list=list(mod_id_event = expression(paste(mod_id_visit, samp, sep = '_'))), #list of variables in sampling time data
                          check.path=checkDir, 
                          failmiss_fn = 'MissingConcDate-', #filename stub for missing concentration data check file
                          multsets_fn = 'multipleSetsConc-', #filename stub for multiple sets of concentration data check file
                          faildup_fn = 'DuplicateConc-', #filename stub for duplicate concentrations check file
                          drugname=drugname, #drug name of interest (we previously defined this as 'fent')
                          LLOQ=LLOQ, #lower limit of quantification (we previously defined this as 0.05 ng/mL)
                          demo.list=demo.out) #file name for processed demographics file use to exclude subjects based on exclusionary criteria

head(conc.out)

# The output in the console provides a message stating 3 rows are missing concentration date

# Load in file containing the 3 records with missing values for the date.time 
fail.miss.conc.date <- read.csv(file.path(checkDir,"failMissingConcDate-fent.csv"))
fail.miss.conc.date

# Correct missing concentration dates
fail.miss.conc.date[,"date.time"] <- c("9/30/2016 09:32","10/1/2016 19:20","10/2/2016 02:04")
fail.miss.conc.date

# Create file with missing concentration dates
write.csv(fail.miss.conc.date, file.path(checkDir,"fixMissingConcDate-fent.csv"))

# Re-run run_DrugLevel() to incorporate missing values
conc.out <- run_DrugLevel(conc.path=file.path(dataDir,"Fentanyl_conc_mod_id.rds"),
                          conc.select=c('mod_id','mod_id_visit','samp','fentanyl_calc_conc'),
                          conc.rename=c(fentanyl_calc_conc = 'conc.level', samp= 'event'),
                          conc.mod.list=list(mod_id_event = expression(paste(mod_id_visit, event, sep = '_'))),
                          samp.path=file.path(dataDir,"Fentanyl_samp_mod_id.rds"),
                          samp.mod.list=list(mod_id_event = expression(paste(mod_id_visit, samp, sep = '_'))),
                          check.path=checkDir,
                          failmiss_fn = 'MissingConcDate-',
                          multsets_fn = 'multipleSetsConc-',
                          faildup_fn = 'DuplicateConc-', 
                          drugname=drugname,
                          LLOQ=LLOQ,
                          demo.list=demo.out)

# Output is a dataset for processed drug concentration levels (conc.level) matched with data/time (date.time) with necessary identification numbers for further merging with other modules.

##### Pro-Laboratory #####

# Processes laboratory data that can be merged with other modules. 

creat.out <- run_Labs(lab.path=file.path(dataDir,"Fentanyl_creat_mod_id.rds"), #file path where creatinine data exists
                      lab.select = c('mod_id','date.time','creat'), #list of variables of interest
                      lab.mod.list = list(date.time = expression(parse_dates(fixDates(paste(date, time)))))) #list containing modifications to variables in the laboratory data

alb.out <- run_Labs(lab.path=file.path(dataDir,"Fentanyl_alb_mod_id.rds"), #file path where albumin data exists
                    lab.select = c('mod_id','date.time','alb'),
                    lab.mod.list = list(date.time = expression(parse_dates(fixDates(paste(date, time))))))

lab.out <- list(creat.out, alb.out)

str(lab.out)

##### Build PK-IV (NONMEM Format) #####

# Create PK Data in NONMEM format for IV medications. Dose data from Pro-Med-Str module and concentration data from the Pro-DrugLevel module are required. Data from other modules are optional. Similar to above, run_Build_PK_IV() will generate several files to check potential data errors and allow for you to provide feedback via corrected 'fix' files that you provide.

pk_dat <- run_Build_PK_IV(conc=conc.out, #concentration date output from run_DrugLevel()
                          dose=ivdose.out, #IV dose data output from run_MedStrI()
                          demo.list=demo.out, #file name for processed demographic file
                          demo.vars=c('weight', 'weight_demo', 'height', 'gender', #vector of demographic variables for final output
                                      'ageatsurgery', 'stat_sts', 'cpb_sts',
                                      'length_of_icu_stay'),
                          demo.abbr=c('wgt', 'wgt_demo', 'height', 'gender', #vector of abbreviations for demo.vars
                                      'age', 'stat', 'cpb', 'loi'),
                          lab.dat = lab.out, #list containing processed laboratory files
                          lab.vars = c('creat','alb'), #vector containing names for laboratory files
                          pk.vars=c('mod_id_visit', 'time', 'conc', 'dose', 'rate', 'event', #PK variables to include
                                    'other', 'multiple.record', 'date', 'mod_id'),
                          drugname=drugname, #drug name stub (fent)
                          check.path=checkDir, #file path where generated files for data checking are stored
                          missdemo_fn='-missing-demo', #file name stub for report of missingness frequency and percent for variables
                          faildupbol_fn='DuplicateBolus-', #filename stub for duplicate bolus dose records
                          date.format="%m/%d/%y %H:%M:%S", #date and time format
                          date.tz="America/Chicago") #date timezone

head(pk_dat)

# convert id back to original IDs
pk_dat <- pullRealId(pk_dat, remove.mod.id=TRUE)

head(pk_dat)

##### Build-PK-Oral (NONMEM Format) #####

# Create PK Data in NONMEM format for oral medications. Dose data from Pro-Med-Str module and concentration data from the Pro-DrugLevel module are required. Data from other modules are optional.

# Data generating function for examples
mkdat <- function() {
  npat=3
  visits <- floor(runif(npat, min=2, max=6))
  id <- rep(1:npat, visits)
  dt <- as.POSIXct(paste(as.Date(sort(sample(700, sum(visits))), 
                                 origin = '2019-01-01'), '10:00:00'), tz = 'UTC') 
  + rnorm(sum(visits), 0, 1*60*60)
  dose_morn <- sample(c(2.5,5,7.5,10), sum(visits), replace = TRUE)
  conc <- round(rnorm(sum(visits), 1.5*dose_morn, 1),1)
  ld <- dt - sample(10:16, sum(visits), replace = TRUE) * 3600
  ld[rnorm(sum(visits)) < .3] <- NA
  age <- rep(sample(40:75, npat), visits)
  weight <- rep(round(rnorm(npat, 180, 20)),visits)
  hgb <- round(rep(rnorm(npat, 10, 2), visits),1)
  cyp2c19 <- sample(c(1,2,3,4,5), sum(visits), replace=TRUE)
  cyp2d6 <- sample(c(1,2,3,4), sum(visits), replace=TRUE)
  data.frame(id, dt, dose_morn, conc, age, weight, hgb, cyp2c19, cyp2d6, ld)
}

# Make example data
set.seed(30)
dat <- mkdat()
dat

dat2 <- dat[,-10] #remove column with last dose ('ld') values
# Build PK data without last-dose times
run_Build_PK_Oral(x = dat2, #data
                  idCol = "id", #patient ID
                  dtCol = "dt", #time of concentration measurement
                  doseCol = "dose_morn", #dose
                  concCol = "conc", #drug concentration
                  ldCol = NULL, #last-dose time, with default set to NULL to ignore
                  first_interval_hours = 336, #hours of regular dosing leading up to the first concentration, with default set at 336 hours = 14 days. Specifying this argument assumes that all individuals were following a regular dosing schedule leading up to the first extracted concentration, and should therefore be long enough for trough concentrations to reach steady state.
                  imputeClosest = NULL) #vector of columns for imputation of missing data using last observation carried forward or, if unavailable, next observation propagated backward

# Build PK data with last-dose times
run_Build_PK_Oral(x = dat,
                  idCol = "id",
                  dtCol = "dt",
                  doseCol = "dose_morn",
                  concCol = "conc",
                  ldCol = "ld",
                  first_interval_hours = 336,
                  imputeClosest = NULL)

# In comparing run_Build_PK_Oral() with and without 'ldCol' specified, we can see several differences in rows 7-9. The measured concentration of 14.1 on data 2019-11-01 is associated with a last-dose time, addl drops from 69 to 68, and the extracted last-dose is added in row 8 with it's corresponding date. The number of doses leading up to the concentrations remain unchanged and the timing of the final dose was adjusted to reflect information in the EHR. NOTE: THIS DATASET STILL RELIES ON ASSUMPTIONS REGARDING DOSING, BUT SHOULD REFLECT THE ACTUAL DOSING SCHEDULE BETTER BY INCORPORATING LAST-DOSE TIME FROM THE EHR.



