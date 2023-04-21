library(EHR)
library(medExtractR)

##### Extract-Med Module #####

# Load in notes
tac_notes <- list(
  "tacronote1.txt",
  "tacronote2.txt",
  "tacronote3.txt",
  "tacronote4.txt"
)

# Execute Extract-Med Module
tac_mxr <- extractMed(note_fn = tac_notes, #File name or vector/list of notes
                      drugnames = c("tacrolimus", "prograf", #Names of drugs we want to extract
                                    "tac", "tacro", 
                                    "fk", "fk506"), 
                      drgunit = "mg", #Unit of drug
                      windowlength = 60, #Length of search window after drug name
                      max_edit_dist = 2, #Define edit distance to capture misspellings in drug names (e.g tcarolimius)
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