library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(genefilter)
library(survival)
library(gplots)

query = GDCquery(
                project = "TCGA-LGG", 
                 data.category = "Transcriptome Profiling", 
                 experimental.strategy = "RNA-Seq", 
                 workflow.type = "HTSeq - Counts") 

res = getResults(query)

res$sample_type <- as.factor(res$sample_type)

#Summarize Sample Types
summary(res$sample_type)



GDCdownload(query)

data = GDCprepare(query)

#View data in a separate R window
View(data@colData)

table(data@colData$vital_status)
table(data@colData$definition)
table(data@colData$tissue_or_organ_of_origin)
table(data@colData$gender)
table(data@colData$race)
table(data@colData$tissue_or_organ_of_origin)
table(data@colData$age_at_diagnosis)

cl = data@colData

#Keep the data-frame to primary solid tumor
clin = cl[cl$definition == "Primary solid Tumor", 
          c("patient", 
            "vital_status", 
            "days_to_death", 
            "days_to_last_follow_up", 
            "gender",
            "race",
            "tissue_or_organ_of_origin",
            "age_at_diagnosis")]


#Separate column for "if dead"
clin$dead = clin$vital_status == "Dead"

#New variable to create days left until end result
clin$overall_survival = ifelse(clin$dead,
                               clin$days_to_death,
                               clin$days_to_last_follow_up)

#New variable to describe whether the cancer is in the cerebrum or not
clin$cerebrum = ifelse(clin$tissue_or_organ_of_origin == "Cerebrum",
                       clin$tissue_or_organ_of_origin,
                       "Not Cerebrum")

#New variable to describe whether the patient is above the age of 40 or not
clin$aboveForty = ifelse(clin$age_at_diagnosis > 14600,
                         TRUE,
                         FALSE)



#Model with "aboveForty" strata
Fortyfit = survfit(Surv(clin$overall_survival, clin$dead) ~ clin$aboveForty)

print(Fortyfit)

#Kaplan-Meier plot for aboveForty
ggsurvplot(Fortyfit, data=clin, pval=T)



#Model with "Gender" strata
Genderfit = survfit(Surv(clin$overall_survival, clin$dead) ~ clin$gender)

print(Genderfit)

#Kaplan-Meier plot for Gender
ggsurvplot(Genderfit, data=clin, pval=T)




#Model with "Race" strata
Racefit = survfit(Surv(clin$overall_survival, clin$dead) ~ clin$race)

print(Racefit)

#Kaplan-Meier plot for Race
ggsurvplot(Racefit, data=clin, pval=T)



#Model with "Cerebrum" strata
Cerebrumfit = survfit(Surv(clin$overall_survival, clin$dead) ~ clin$cerebrum)

print(Cerebrumfit)

#Kaplan-Meier plot for Cerebrum
ggsurvplot(Cerebrumfit, data=clin, pval=T)




