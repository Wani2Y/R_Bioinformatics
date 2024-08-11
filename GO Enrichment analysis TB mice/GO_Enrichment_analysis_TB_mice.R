# GO_Enrichment_analysis_TB_mice

#This notebook is designed to explore the molecular data published at DOI: 10.1128/iai.00263-23; This notebook is free to use; you can redistribute it and/or modify it under the terms of the GNU General Public License,version 3, as published by the Free Software Foundation. This notebook is distributed in the hope that it will serve educational and entertainment purposes, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. Use of this code may be subjected to licensing agreement and copyright requirement of the Journal, Infection and Immunity.

#For more details of the license, please refer to: www.gnu.org/licenses/gpl-3.0.en.html.

# Environment setup
# uncomment the following lines of this chunk if you have not installed all required packages. Modify the code if you already have some of the packages.
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("affy", "limma", "annotate", "oligo", "Biobase", "R.utils", "dplyr", "tidyverse", "clusterProfiler", "ggplot2", "enrichplot", "reshape2", "plotly"))

# Microarray data procurement and cleaning
library(affy)
library(limma)
library(preprocessCore)
library(Biobase)
library(R.utils)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(reshape2)
library(plotly)

## Downloading microarray data from GEO
cell_sample <- data.frame(cbind(
  Ex_group = c(rep("uninfected", 5),
               rep("asymptomatic_controller", 14), 
               rep("chronic_controller", 17), #on GEO, it is labelled as "Controller".
               rep("progressor", 10), # on GEO, it is labelled as "Supersusceptible Progressor".
               rep("uninfected", 20),
               rep("asymptomatic_controller", 39),
               rep("progressor", 29)
  ), 
  
  cell_series = c("GSM8250425", "GSM8250426", "GSM8250427", "GSM8250428", "GSM8250429", # uninfected_1
                  "GSM8250430", "GSM8250431", "GSM8250432", "GSM8250433", "GSM8250434", # asymptomatic controller_1
                  "GSM8250435", "GSM8250436", "GSM8250437", "GSM8250438", "GSM8250439", 
                  "GSM8250440", "GSM8250441", "GSM8250442", "GSM8250443", 
                  "GSM8250444", "GSM8250445", "GSM8250446", "GSM8250447", "GSM8250448", # chronic controller_1
                  "GSM8250449", "GSM8250450", "GSM8250451", "GSM8250452", "GSM8250453", 
                  "GSM8250454", "GSM8250455", "GSM8250456", "GSM8250457", "GSM8250458", 
                  "GSM8250459", "GSM8250460", 
                  "GSM8250461", "GSM8250462", "GSM8250463", "GSM8250464", "GSM8250465", # progressor_1
                  "GSM8250466", "GSM8250467", "GSM8250468", "GSM8250469", "GSM8250470", 
                  "GSM8250471", "GSM8250472", "GSM8250473", "GSM8250474", "GSM8250475", # uninfected_2
                  "GSM8250476", "GSM8250477", "GSM8250478", "GSM8250479", "GSM8250480", 
                  "GSM8250481", "GSM8250482", "GSM8250483", "GSM8250484", "GSM8250485",
                  "GSM8250486", "GSM8250487", "GSM8250488", "GSM8250489", "GSM8250490", 
                  "GSM8250491", "GSM8250492", "GSM8250493", "GSM8250494", "GSM8250495", # asymptomatic controller_2
                  "GSM8250496", "GSM8250497", "GSM8250498", "GSM8250499", "GSM8250500",
                  "GSM8250501", "GSM8250502", "GSM8250503", "GSM8250504", "GSM8250505",
                  "GSM8250506", "GSM8250507", "GSM8250508", "GSM8250509", "GSM8250510",
                  "GSM8250511", "GSM8250512", "GSM8250513", "GSM8250514", "GSM8250515",
                  "GSM8250516", "GSM8250517", "GSM8250518", "GSM8250519", "GSM8250520",
                  "GSM8250521", "GSM8250522", "GSM8250523", "GSM8250524", "GSM8250525",
                  "GSM8250526", "GSM8250527", "GSM8250528", "GSM8250529",
                  "GSM8250530", "GSM8250531", "GSM8250532", "GSM8250533", "GSM8250534", # progressor_2
                  "GSM8250535", "GSM8250536", "GSM8250537", "GSM8250538", "GSM8250539",
                  "GSM8250540", "GSM8250541", "GSM8250542", "GSM8250543", "GSM8250544",
                  "GSM8250545", "GSM8250546", "GSM8250547", "GSM8250548", "GSM8250549",
                  "GSM8250550", "GSM8250551", "GSM8250552", "GSM8250553", "GSM8250554",
                  "GSM8250555", "GSM8250556", "GSM8250557", "GSM8250558"
  ), 
  
  file_tag = c("5F161209", "5F161209", "5F161209", "5F161209", "5F161209", # uninfected_1
               "5F161209", "5F161209", "5F161209", "5F161209", "5F161209",  # asymptomatic controller_1
               "5F161209", "5F161209", "5F161209", "5F191212", "5F161209", 
               "5F161209", "5F161209", "5F161209", "5F191212",
               "5F161209", "5F161209", "5F161209", "5F191212", "5F161209", # chronic controller_1
               "5F161209", "5F161209", "5F191212", "5F161209", "5F161209", 
               "5F191212", "5F161209", "5F191212", "5F191212", "5F191212", 
               "5F191212", "5F191212", 
               "5F161209", "5F161209", "5F161209", "5F161209", "5F161209", # progressor_1
               "5F161209", "5F161209", "5F161209", "5F161209", "5F161209", 
               "5F161209", "5F161209", "5F191212", "5F170526", "5F170526", # uninfected_2
               "5F161209", "5F191212", "5F170526", "5F161209", "5F191212", 
               "5F170526", "5F161209", "5F191212", "5F170526", "5F191212", 
               "5F170526", "5F191212", "5F170526", "5F191212", "5F170526", 
               "5F170526", "5F170526", "5F170526", "5F170526", "5F170526", # asymptomatic controller_2
               "5F191212", "5F170526", "5F170526", "5F191212", "5F191212",
               "5F191212", "5F191212", "5F170526", "5F191212", "5F170526",
               
               "5F170526", "5F191212", "5F170526", "5F191212", "5F170526",
               "5F191212", "5F170526", "5F170526", "5F170526", "5F191212",
               "5F170526", "5F191212", "5F170526", "5F191212", "5F170526",
               "5F170526", "5F170526", "5F170526", "5F170526", "5F191212",
               "5F191212", "5F170526", "5F191212", "5F191212",
               "5F170526", "5F170526", "5F191212", "5F170526", "5F191212", # progressor_2
               "5F191212", "5F191212", "5F170526", "5F191212", "5F191212",
               "5F191212", "5F170526", "5F191212", "5F191212", "5F170526",
               "5F170526", "5F170526", "5F170526", "5F191212", "5F170526",
               "5F170526", "5F170526", "5F170526", "5F170526", "5F170526",
               "5F191212", "5F191212", "5F191212", "5F170526"
  ), 
  
  grp_tag = c("E1", "E2", "E3", "E4", "E5", # uninfected_1
              "A1", "A2", "A3", "A4", "A6",  # asymptomatic controller_1
              "A7", "B1", "B5", "A1", "B7", 
              "C3", "C7", "C8", "B2",
              "A5", "B3", "B4", "A2", "C1", # chronic controller_1
              "C4", "C5", "A3", "D1", "D3", 
              "A5", "D5", "A6", "A7", "A8", 
              "B1", "B3", 
              "A8", "B2", "B6", "B8", "C2", # progressor_1
              "C6", "D2", "D6", "D7", "D8", 
              "E6", "E7", "F2", "F1", "F2", # uninfected_2
              "E8", "F3", "F3", "F1", "F4", 
              "F4", "F2", "F5", "F5", "F6", 
              "F6", "F7", "F7", "F8", "F8", 
              "A2", "A4", "A7", "B1", "B2", # asymptomatic controller_2
              "C3", "B3", "B5", "C6", "C7",
              "C8", "D2", "C2", "D4", "C3",
              "C4", "D5", "C7", "D6", "C8",
              "D7", "D1", "D2", "D4", "D8",
              "D6", "E1", "D7", "E2", "E2",
              "E3", "E4", "E5", "E6", "E5",
              "E6", "E7", "E7", "F1",
              "A1", "A3", "B4", "A5", "B5", # progressor_2
              "B6", "B7", "A6", "B8", "C1",
              "C2", "A8", "C5", "D1", "B6",
              "B7", "B8", "C1", "D3", "C5",
              "C6", "D3", "D5", "D8", "E1",
              "E3", "E4", "E8", "E8"
  ), 
  
  file_name = c("F249", "F250", "F251", "F252", "F253", # uninfected_1
                "F194", "F195", "F196", "F197", "F199", # asymptomatic controller_1
                "F200", "F202", "F206", "F206", "F208", 
                "F213", "F217", "F218", "F236", 
                "F198", "F204", "F205", "F209", "F211", # chronic controller_1
                "F214", "F215", "F215", "F220", "F222", 
                "F225", "F227", "F229", "F232", "F233", 
                "F234", "F241", 
                "F201", "F203", "F207", "F210", "F212", # progressor_1
                "F216", "F221", "F231", "F239", "F240", 
                "F474", "F475", "F475", "F475", "F476", # uninfected_2
                "F477", "F477", "F477", "F478", "F478", 
                "F478", "F479", "F479", "F479", "F480", 
                "F480", "F481", "F481", "F482", "F482", 
                "F260", "F265", "F288", "F296", "F300", # asymptomatic controller_2
                "F303", "F305", "F311", "F318", "F329", 
                "F336", "F359", "F359", "F367", "F368",
                "F372", "F378", "F384", "F386", "F386",
                "F389", "F394", "F398", "F409", "F415",
                "F417", "F418", "F418", "F422", "F424",
                "F426", "F429", "F440", "F443", "F451",
                "F452", "F455", "F459", "F467",
                "F256", "F261", "F267", "F272", "F274", # progressor_2
                "F279", "F281", "F285", "F286", "F295",
                "F299", "F299", "F316", "F340", "F344",
                "F350", "F351", "F355", "F365", "F375",
                "F383", "F407", "F411", "F419", "F421",
                "F434", "F437", "F465", "F471"
  ))) 

# set row names for subsequent indexing.
rownames(cell_sample) <- c("GSM8250425.CEL", "GSM8250426.CEL", "GSM8250427.CEL", "GSM8250428.CEL", "GSM8250429.CEL", # uninfected_1
                           "GSM8250430.CEL", "GSM8250431.CEL", "GSM8250432.CEL", "GSM8250433.CEL", "GSM8250434.CEL", # asymptomatic controller_1
                           "GSM8250435.CEL", "GSM8250436.CEL", "GSM8250437.CEL", "GSM8250438.CEL", "GSM8250439.CEL", 
                           "GSM8250440.CEL", "GSM8250441.CEL", "GSM8250442.CEL", "GSM8250443.CEL", 
                           "GSM8250444.CEL", "GSM8250445.CEL", "GSM8250446.CEL", "GSM8250447.CEL", "GSM8250448.CEL", # chronic controller_1
                           "GSM8250449.CEL", "GSM8250450.CEL", "GSM8250451.CEL", "GSM8250452.CEL", "GSM8250453.CEL", 
                           "GSM8250454.CEL", "GSM8250455.CEL", "GSM8250456.CEL", "GSM8250457.CEL", "GSM8250458.CEL", 
                           "GSM8250459.CEL", "GSM8250460.CEL", 
                           "GSM8250461.CEL", "GSM8250462.CEL", "GSM8250463.CEL", "GSM8250464.CEL", "GSM8250465.CEL", # progressor_1
                           "GSM8250466.CEL", "GSM8250467.CEL", "GSM8250468.CEL", "GSM8250469.CEL", "GSM8250470.CEL", 
                           "GSM8250471.CEL", "GSM8250472.CEL", "GSM8250473.CEL", "GSM8250474.CEL", "GSM8250475.CEL", # uninfected_2
                           "GSM8250476.CEL", "GSM8250477.CEL", "GSM8250478.CEL", "GSM8250479.CEL", "GSM8250480.CEL", 
                           "GSM8250481.CEL", "GSM8250482.CEL", "GSM8250483.CEL", "GSM8250484.CEL", "GSM8250485.CEL", 
                           "GSM8250486.CEL", "GSM8250487.CEL", "GSM8250488.CEL", "GSM8250489.CEL", "GSM8250490.CEL", 
                           "GSM8250491.CEL", "GSM8250492.CEL", "GSM8250493.CEL", "GSM8250494.CEL", "GSM8250495.CEL", # asymptomatic controller_2
                           "GSM8250496.CEL", "GSM8250497.CEL", "GSM8250498.CEL", "GSM8250499.CEL", "GSM8250500.CEL",
                           "GSM8250501.CEL", "GSM8250502.CEL", "GSM8250503.CEL", "GSM8250504.CEL", "GSM8250505.CEL",
                           "GSM8250506.CEL", "GSM8250507.CEL", "GSM8250508.CEL", "GSM8250509.CEL", "GSM8250510.CEL",
                           "GSM8250511.CEL", "GSM8250512.CEL", "GSM8250513.CEL", "GSM8250514.CEL", "GSM8250515.CEL",
                           "GSM8250516.CEL", "GSM8250517.CEL", "GSM8250518.CEL", "GSM8250519.CEL", "GSM8250520.CEL",
                           "GSM8250521.CEL", "GSM8250522.CEL", "GSM8250523.CEL", "GSM8250524.CEL", "GSM8250525.CEL",
                           "GSM8250526.CEL", "GSM8250527.CEL", "GSM8250528.CEL", "GSM8250529.CEL",
                           "GSM8250530.CEL", "GSM8250531.CEL", "GSM8250532.CEL", "GSM8250533.CEL", "GSM8250534.CEL", # progressor_2
                           "GSM8250535.CEL", "GSM8250536.CEL", "GSM8250537.CEL", "GSM8250538.CEL", "GSM8250539.CEL",
                           "GSM8250540.CEL", "GSM8250541.CEL", "GSM8250542.CEL", "GSM8250543.CEL", "GSM8250544.CEL",
                           "GSM8250545.CEL", "GSM8250546.CEL", "GSM8250547.CEL", "GSM8250548.CEL", "GSM8250549.CEL",
                           "GSM8250550.CEL", "GSM8250551.CEL", "GSM8250552.CEL", "GSM8250553.CEL", "GSM8250554.CEL",
                           "GSM8250555.CEL", "GSM8250556.CEL", "GSM8250557.CEL", "GSM8250558.CEL"
) 

# download the gene data from GEO
geo <- "https://www.ncbi.nlm.nih.gov/geo/download/"

# iterate through all download links of each sample.
for (ID in cell_sample$cell_series){
  # paste a string for the url
  url <- paste0(geo, "?acc=", 
                ID, 
                "&format=file&file=", 
                ID, 
                "%5FGB%",
                cell_sample[ID, 'file_tag'],
                "%5F", 
                cell_sample[ID, 'grp_tag'], 
                "%5", 
                cell_sample[ID, 'file_name'], 
                "%2ECEL%2Egz")
  # create a file pointer
  fp <- paste0("gene_data/", 
               ID, 
               ".CEL.gz") 
  # check if files already exist
  if (!file.exists(fp)){
    download.file(url, fp)
    Sys.sleep(10) # wait for 5 - 10 seconds between downloads as recommended by GEO. Be nice to the server!
  }
}

# load the microarray data
cell_fp <- "gene_data/"
cell_name <- list.files(cell_fp, pattern = "\\.CEL\\.gz$", full.names = TRUE)
cell_file <- gsub("\\.gz$", "", cell_name)
for (i in seq_along(cell_name)) {
  if (!file.exists(cell_file[i])) {
    gunzip(cell_name[i], destname = cell_file[i])
  }
}

file_list <- list.celfiles(cell_fp, full.names = TRUE)

cell <- oligo::read.celfiles(cell_file)

# normalise with Robust Multi-Array Average expression(rma).
cell_norm <- rma(cell)

# modify cell_sample into an annotated dataframe.
pheno_data <- AnnotatedDataFrame(data = cell_sample)

# wrangle the data into an expression data.
exp_set <- ExpressionSet(assayData = exprs(cell_norm), 
                         phenoData = pheno_data)

# perform Differential Expressed Gene Analysis (DEGA)
## contrast of all four groups
exp_set_all <- exp_set
pData(exp_set_all)$Ex_group <- factor(pData(exp_set_all)$Ex_group, 
                                      levels = c("uninfected", 
                                                 "asymptomatic_controller", 
                                                 "chronic_controller", 
                                                 "progressor"))

mod_design <- model.matrix(~ Ex_group, data = pData(exp_set_all))
colnames(mod_design) <- c("(Intercept)",
                          "chronic_controller",
                          "progressor",
                          "uninfected" )
mod_fit_1 <- eBayes(lmFit(exp_set_all, design = mod_design))

# grab the p-value for all probes. Set a number if you cannot handle the size of dataset.
probe_results_all <- topTable(mod_fit_1, number = Inf)

# contrast between any two groups
two_grp_degs <- function(group_1, group_2){
  # convert input to character.
  grp_1 <- as.character(group_1)
  grp_2 <- as.character(group_2)
  
  # get the phenodata.
  pheno <- pData(exp_set)[pData(exp_set)$Ex_group %in% c(grp_1, grp_2), ]
  exp_subset <- exp_set[, rownames(pheno)]
  
  # fit a new model
  mod_design <- model.matrix(~ Ex_group, data = pheno)
  mod_fit <- eBayes(lmFit(exp_subset, design = mod_design))
  
  # the the p-value for relevant probes. Modify this line if you can not handle a lot of data.
  probe_results <- topTable(mod_fit, number = Inf)
  
  # make a volcano plot to see how many probes are differentially expressed and whether they are up or down regulated.
  p_title <- paste0(toupper(substring(grp_1, 1, 1)), tolower(substring(grp_1, 2)), "  VS  ",
                    toupper(substring(grp_2, 1, 1)), tolower(substring(grp_2, 2)))
  x <- probe_results$logFC
  y <- -log10(probe_results[, "adj.P.Val"])
  p <- ggplot(probe_results, aes(x, y)) +
    geom_point() +
    labs(
      x = "LogFC",
      y = "Negative log10 adjusted p-value",
      title =  p_title)+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
    theme_minimal()
  
  # put results together in a list
  result <- list(probe_results, p)
  names(result) <- c("deg_results", "plot")
  
  # return the data as list for subsequent cbind
  return(result)
}

# Uninfected vs Aymptomatic Controller
g <- c("uninfected", "asymptomatic_controller")
deg_ui_ac <- two_grp_degs(g[1], g[2])

# Uninfected vs Chronic Controller
g <- c("uninfected", "chronic_controller")
deg_ui_cc <- two_grp_degs(g[1], g[2])

# Uninfected vs Progressor
g <- c("uninfected", "progressor")
deg_ui_pr <- two_grp_degs(g[1], g[2])

# Aymptomatic Controller vs Chronic Controller
g <- c("asymptomatic_controller", "chronic_controller")
deg_ac_cc <- two_grp_degs(g[1], g[2])

# Aymptomatic Controller vs Progressor
g <- c("asymptomatic_controller", "progressor")
deg_ac_pr <- two_grp_degs(g[1], g[2])

# Chronic Controller vs Progressor
g <- c("chronic_controller", "progressor")
deg_cc_pr <- two_grp_degs(g[1], g[2])

# Enrichment/Overrepresentation Analysis
# download from GEO.
options(timeout = 600) # the file is about 240 mb.
fp <- "gene_data/GPL16570-1802.txt"
url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GPL16570&id=1802&db=GeoDb_blob100"

if (!file.exists(fp)){
  download.file(url, fp)}

# load the text file as a table, and omit the metadata
GPL16570 <- read.table("gene_data/GPL16570-1802.txt", header = TRUE, sep = "\t",
                       skip = 0, quote = "", row.names = 1,fill = TRUE)
GPL16570_lung <- GPL16570 %>% 
  filter(grepl("lung", unigene, ignore.case = TRUE)) 

# Create four sets of Gene Ontology (GO) terms
# make a function to grab the second and third items from gene_assignment.
parseGI <- function(gene_assignment){
  str_ful <- str_match(gene_assignment, "// (\\w+) // (.+?) //")
  gene_id <- str_ful[2]
  gene_title <- str_ful[3]
  # return the data as list for subsequent cbind actions.
  return(list(Gene_id = gene_id, Gene_title = gene_title))
}

# parsing GO terms
# biological process is chosen as the source article focus on the immunological aspect.
parsGO <- function(probe, ontology = "GO_biological_process"){
  go_segment <- strsplit(GPL16570_lung[probe, "GO_biological_process"], "///")[[1]]
  go_id <- trimws(unlist(lapply(go_segment, function(x) strsplit(x, "//")[[1]][2])))
  go_description <- trimws(unlist(lapply(go_segment, function(x) strsplit(x, "//")[[1]][3])))
  go_source <- trimws(unlist(lapply(go_segment, function(x) strsplit(x, "//")[[1]][4])))
  
  # in case of probe without id (a.k.a. errors in data wrangling)
  if (length(go_id) > 0){
    go_tab <- cbind(
      Id = go_id,
      Description = go_description,
      Source = go_source)
    rownames(go_tab) <- rep(probe, nrow(go_tab))
    return(go_tab)
    return() # return nothing if go_id is empty/NA
  }
}

# prepare AGO and GO data
prepGO <- function(degs){
  # match information between Array GO terms from DEG result and GPL16750_lung.
  AGO <- cbind(Gene_Symbol = GPL16570_lung[rownames(degs), c("probeset_id", "Gene_id", "Gene_title")],
               degs[, c("logFC", "P.Value", "adj.P.Val")])
  colnames(AGO) <- c("Probe_id", "Gene_id", "Gene_title", "logFC", "P.Val", "adj.P.Val")
  
  # parse all GO terms on Biological Process column from GPL16570. 
  GO <- c()
  
  for (probe in rownames(AGO)){
    GO <- rbind(GO, parsGO(probe))}
  
  # drop the rows with only na.
  GO <- GO[rowSums(is.na(GO)) != ncol(GO), ]
  
  # match AGO and GO, as DEGs without GO term is of little significant to our purpose.
  AGO <- AGO[which(rownames(AGO) %in% rownames(GO)), ]
  
  # combine the result into a list for subsequent enrichment analysis.
  result <- list(AGO, GO)
  names(result) <- c("AGO", "GO")
  
  return(result)
}

# modify GPL16570_lung in place.
GPL16570_lung <- GPL16570_lung %>%
  mutate(Gene_id = map_chr(gene_assignment, ~parseGI(.x)$Gene_id),
         Gene_title = map_chr(gene_assignment, ~parseGI(.x)$Gene_title))

# Uninfected vs Aymptomatic Controller
degs <- deg_ui_ac$deg_results
Genego_ui_ac <- prepGO(degs)

# Uninfected vs Chronic Controller
degs <- deg_ui_cc$deg_results
Genego_ui_cc <- prepGO(degs)

# Uninfected vs Progressor
degs <- deg_ui_pr$deg_results
Genego_ui_pr <- prepGO(degs)

# Aymptomatic Controller vs Chronic Controller
degs <- deg_ac_cc$deg_results
Genego_ac_cc <- prepGO(degs)

# Aymptomatic Controller vs Progressor
degs <- deg_ac_pr$deg_results
Genego_ac_pr <- prepGO(degs)

# Chronic Controller vs Progressor
degs <- deg_cc_pr$deg_results
Genego_cc_pr <- prepGO(degs)

# Perform the enrichment/overrespresentation analysis
EACompute <- function(ago, go){
  # assign the arguments.
  AGO <-  ago
  GOdata <- go
  
  # check the minimal adj.P.Val.
  if (min(AGO$adj.P.Val) < 0.05){
    # perform the EA computation.
    # compute for foreground.
    fg_probe <- rownames(AGO)[which(AGO[, "adj.P.Val"] < 0.05)]
    fg_go_id <- GOdata[fg_probe, "Id"]
    
    # compute for background.
    bg_probe <- rownames(AGO)
    bg_go_id <- GOdata[bg_probe, "Id"]
    
    # compute for the counts for relevant terms form the theoretical section for each gene ontology.
    EA <- c() # enrichment analysis
    
    # counting for each GO term.
    for (Id in unique(fg_go_id)){
      m <- length(unique(names(bg_go_id)[which(bg_go_id == Id)])) # get the count of probe names for each go
      n <- length(unique(names(bg_go_id))) - m # probes in the background lacking go
      q <- length(unique(names(fg_go_id)))
      x <- length(unique(names(fg_go_id)[which(fg_go_id == Id)]))
      p_val <- phyper(x - 1, m, n, q, lower.tail = FALSE) # x - 1 accounts for the discrete nature of the PMF when flipping 
      EA <- rbind(EA, c(m, m + n, q, x, p_val))
    }
    rownames(EA) = unique(fg_go_id)
    colnames(EA) = c("m", "o(m + n)", "q", "x", "P.Val.hygeo")
    
    if(EA[, "x"]> 0){
      # adjust the P.Val.hygeo for multiple testing across the entire array of probes. Finding true positive is more important to me in this analysis, but use one of the other adjusting methods as needed.
      EA <- cbind(EA, "adj.P.Val.hygeo" = p.adjust(EA[, "P.Val.hygeo"], "BH"))
      
      # filter again for only the significant results by adj.P.Val.hygeo. We use 0.05 here again as oup alpha level.
      GOsig <- rownames(EA)[which(EA[, "adj.P.Val.hygeo"] < 0.05)]
      
      # grab and add GO data.
      GOname <- c()
      
      for (go in GOsig){
        GOname <- c(GOname,
                    unique(GOdata[which(GOdata[, "Id"] == go), "Description"]))}
      
      EA <- cbind(Description = GOname, EA[GOsig, ])
      
      # we can format the enrichment analysis output as observed and expected, which may be easier to interpret.
      Ratio_fg <- as.numeric(EA[, "x"])/as.numeric(EA[, "q"])
      
      # background ratio.
      Ratio_bg <- as.numeric(EA[, "m"])/as.numeric(EA[, "o(m + n)"])
      
      # expected value of foreground based on the background ratio.
      Ratio_mu_bg <- round(Ratio_bg * as.numeric(EA[, "q"]), 5)
      
      # create a dataframe as a readable result.
      EA <- data.frame(cbind(`GO id` = rownames(EA),
                             Description = EA[, "Description"],
                             `Expected Frequency` = paste0(Ratio_mu_bg, "/", EA[, "q"]),
                             `Observed Frequency` = paste0(EA[, "m"], "/", EA[, "q"]),
                             `Background Ratio` = paste0(EA[, "m"], "/", EA[, "o(m + n)"]),
                             `FoldChange` = round(Ratio_fg/Ratio_bg, 5), # of GO term.
                             `P.Val` = as.numeric(EA[, "P.Val.hygeo"]),
                             `P.Val.adjusted` = as.numeric(EA[, "adj.P.Val.hygeo"])))
      
      # add gene information
      tempt <- EA %>%
        left_join(as.data.frame(cbind(probe_id = rownames(GOdata), GOdata)), by = c("GO.id" = "Id"))
      
      AGO$Probe_id <- as.character(AGO$Probe_id)
      
      EA_join <- tempt %>%
        left_join(AGO, by = c("probe_id" = "Probe_id"))
      
      # remove duplicates and unnecessary columns
      EA_join <- select(EA_join, c("probe_id", "GO.id", "Gene_id",  
                                   "logFC", "Description.x", "FoldChange", 
                                   "P.Val.x", "P.Val.adjusted", "Observed.Frequency", "Background.Ratio"))
      
      colnames(EA_join) <- c("probe_id", "GO.id", "Gene_id", 
                             "logFC", "Description", "FoldChange", 
                             "P.Val", "P.Val.adjusted", "Observed.Frequency", "Background.Ratio")
      
      EA_join <- unique(EA_join)
      
      # remove duplicate for a clean output table.
      tempt <- EA %>%
        select(c("GO.id", "Description", "Expected.Frequency", "Observed.Frequency", 
                 "Background.Ratio", "FoldChange", "P.Val", "P.Val.adjusted"))
      
      join_list <- EA_join %>%
        select(c("GO.id", "Gene_id"))
      
      tempt <- tempt %>%
        left_join(join_list, by = "GO.id") 
      
      tempt <- tempt %>%
        group_by(GO.id) %>%
        dplyr::reframe(
          Gene.id = paste(Gene_id, collapse = "/"))
      
      EA_tab <- tempt %>%
        left_join(EA, by = "GO.id")
      
      # create an enrichResult object for visualisation.
      # create dataframe for encichResult obj.
      enrich_df <- select(EA_join, c("probe_id", "GO.id", "Description", "Gene_id", "logFC", "P.Val", "P.Val.adjusted"))
      tempt <- GPL16570_lung[, c("probeset_id", "Gene_id")]
      tempt$probe_id <- rownames(tempt)
      
      enrich_df <- enrich_df %>%
        left_join(tempt, by = "probe_id") %>%
        select(c("probe_id", "GO.id", "Description", "Gene_id.x", "Gene_id.y", "logFC", "P.Val", "P.Val.adjusted")) %>% 
        na.omit()
      
      # set the column names to avoid confusions.
      colnames(enrich_df) <- c("probe_id", "GO.id", "Description", "AGene.id", "BgGene.id", "logFC", "P.Val", "P.Val.adjusted")
      
      # count the number of genes and concat the gene id column. The numerator of foreground and background gene ratios are supposed to be the same as we filtered probes by GO terms in earlier sections.
      enrich_df <- enrich_df %>%
        group_by(GO.id, Description, P.Val, P.Val.adjusted) %>%
        dplyr::summarize(
          AGene.id = paste(AGene.id, collapse = "/"),
          AGene.count = sapply(strsplit(AGene.id, "/"), length),
          BgGene.id = paste(BgGene.id, collapse = "/"),
          BgGene.count = sapply(strsplit(BgGene.id, "/"), length))
      
      array_gene_list <- unique(EA_join$Gene_id)
      gpl_gene_list <- unique(GPL16570_lung$Gene_id)
      
      # get denominator for computation of gene ratio.
      array_gene_count <- length(unique(EA_join$Gene_id))
      gpl_gene_count <- length(unique(GPL16570_lung$Gene_id))
      
      # wrangle data into a dataframe with correct format.
      ea_df <- data.frame(cbind(
        ID = enrich_df[, "GO.id"],
        Description = enrich_df[, "Description"],
        GeneRatio = paste0(enrich_df$AGene.count, "/", array_gene_count),
        BgRatio = paste0(enrich_df$BgGene.count, "/", gpl_gene_count),
        pvalue = as.numeric(unlist(enrich_df[, "P.Val"])),
        p.adjust = as.numeric(unlist(enrich_df[, "P.Val.adjusted"])),
        qvalue = rep(0.03, nrow(enrich_df)), # a dummy value for plotting, p.adjust are too small for calculating qvalue.
        geneID = enrich_df[, "AGene.id"],
        Count = enrich_df[, "AGene.count"],
        stringsAsFactors = FALSE))
      
      # ensure the column ID are exact.
      colnames(ea_df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
      
      # create a geneset list for mapping GO terms and gene ID.
      geneSets <- setNames(lapply(ea_df$geneID, function(x) unlist(strsplit(x, "/"))), ea_df$ID)
      
      # create enrichResult object.
      EA_obj <- new("enrichResult",
                    result = ea_df,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.02,
                    gene = array_gene_list,
                    universe = character(0),
                    gene2Symbol = character(0),
                    geneSets = geneSets,
                    organism = "Mus musculus",
                    keytype = "SYMBOL",
                    ontology = "BP", # biological process
                    readable = TRUE)
      
      result <- list(EA_tab, EA_obj)
      names(result) <- c("EA_tab", "EA_obj")
      return(result)
    }else print("Not enougn significant GO terms for computation. Increased sample size is required.")
  }else print("GO terms from contrast groups are not significant enough for computation. Set a larger alpha level to proceed.")
}

# Uninfected vs Aymptomatic Controller
EA_ui_ac <- EACompute(Genego_ui_ac$AGO, Genego_ui_ac$GO)
write.csv(as.data.frame(EA_ui_ac$EA_tab), "EA_Fin_ui_ac.csv", row.names = FALSE)

# Uninfected vs Chronic Controller
EA_ui_cc <- EACompute(Genego_ui_cc$AGO, Genego_ui_cc$GO)
write.csv(as.data.frame(EA_ui_cc$EA_tab), "EA_Fin_ui_cc.csv", row.names = FALSE)

# Uninfected vs Progressor
EA_ui_pr <- EACompute(Genego_ui_pr$AGO, Genego_ui_pr$GO)
write.csv(as.data.frame(EA_ui_pr$EA_tab), "EA_Fin_ui_pr.csv", row.names = FALSE)

# Aymptomatic Controller vs Chronic Controller
EA_ac_cc <- EACompute(Genego_ac_cc$AGO, Genego_ac_cc$GO)
write.csv(as.data.frame(EA_ac_cc$EA_tab), "EA_Fin_ac_cc.csv", row.names = FALSE)

# Aymptomatic Controller vs Progressor
EA_ac_pr <- EACompute(Genego_ac_pr$AGO, Genego_ac_pr$GO)
write.csv(as.data.frame(EA_ac_pr$EA_tab), "EA_Fin_ac_pr.csv", row.names = FALSE)

# Chronic Controller vs Progressor
EA_ccpr <- EACompute(Genego_cc_pr$AGO, Genego_cc_pr$GO)
write.csv(as.data.frame(EA_cc_pr$EA_tab), "EA_Fin_cc_pr.csv", row.names = FALSE)