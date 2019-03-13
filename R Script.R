library(curatedTCGAData)
library(TCGAutils)
library(TCGAbiolinks)
library(RCTGAtoolbox)

curatedTCGAData("OV", c("RNA*", "RPPA*", "Mutations"), dry.run = TRUE)

ov <- curatedTCGAData("OV", c("RNA*", "RPPA*", "Mutations"), dry.run = FALSE)
colData(ov)
getClinicalNames("OV")

#subtypes
TCGAutils::getSubtypeMap(ov)
head(metadata(colData(ov))[["subtypes"]]) 

#tcgabiolinks get ov data
ovtcgabio <- GDCquery_clinic("TCGA-OV", type = "clinical", save.csv = TRUE)
#subtypes
subtypes <- PanCancerAtlas_subtypes()
DT::datatable(subtypes,
                +               filter = 'top',
                +               options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
                +               rownames = FALSE)
ovtcgabio <- GDCquery_clinic("TCGA-OV", type = "clinical", save.csv = TRUE)

#RCTGCA
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("RTCGAToolbox") 
library(RTCGAToolbox)
ovrtcg <- getFirehoseData(dataset="OV",
                          forceDownload=TRUE, clinical=TRUE, Mutation=TRUE, RPPAArray = TRUE)


intersect(rownames(colData(ov)), )

# find Differences
discrep <- mapply(setdiff, ov, ovtcgabio)
?setdiff

## uco plots between two datasets (after isolating same patients)

## subset for patients in both subtype datasets (intersect)

## compare subtype variables in both datasets
library(curatedTCGAData)
library(TCGAutils)
library(TCGAbiolinks)
library(RTCGAtoolbox)

curatedTCGAData("OV", c("RNA*", "RPPA*", "Mutations"), dry.run = TRUE)

ov <- curatedTCGAData("OV", c("RNA*", "RPPA*", "Mutations"), dry.run = FALSE)
colData(ov)
getClinicalNames("OV")

#subtypes
TCGAutils::getSubtypeMap(ov)
head(metadata(colData(ov))[["subtypes"]]) 

#tcgabiolinks get ov data
ovtcgabio <- GDCquery_clinic("TCGA-OV", type = "clinical", save.csv = TRUE)
#subtypes
subtypes <- PanCancerAtlas_subtypes()
 DT::datatable(subtypes,
                +               filter = 'top',
                +               options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
                +               rownames = FALSE)
ovtcgabio <- GDCquery_clinic("TCGA-OV", type = "clinical", save.csv = TRUE)

#RCTA
ovrtcg <- getFirehoseData(dataset="OV",
                          forceDownload=TRUE, clinical=TRUE, Mutation=TRUE, RPPAArray = TRUE)
stddata <- getFirehoseRunningDates()
brcaov <- getFirehoseData(dataset="OV",
                          forceDownload=TRUE, clinical=TRUE, Mutation=TRUE)
clinov <- getData(brcaov, "clinical")
names(clinov)

intersect(rownames(colData(ov)), )

#firebrowser
require(FirebrowseR)
cohorts = Metadata.Cohorts(format = "csv")
cancer.Type = cohorts[grep("OV", cohorts$description, ignore.case = T), 1]
print(cancer.Type)
ov.pats = Samples.Clinical(cohort = cancer.Type, format="tsv")
dim(ov.pats)

#cbioportal
install.packages('cgdsr')
install.packages('tcgaretriever')
library('cgdsr')
library('TCGAretriever')
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
ov_study <- getCancerStudies(mycgds)[112,1]
ov_case_list <- getCaseLists(mycgds,ov_study)[1,1]
ov_clinic_data <- get_clinical_data("lusc_tcga_3way_complete")

# find Differences
discrep <- mapply(setdiff, ov, ovtcgabio)
?setdiff

## uco plots between two datasets (after isolating same patients)

## subset for patients in both subtype datasets (intersect)

## compare subtype variables in both datasets
