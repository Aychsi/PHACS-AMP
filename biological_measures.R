#install.packages("haven")
library(haven)
data <- read_sas("/Users/hansoochang/Drexel/PHACS:AMP/data/lbw0070.sas7bdat")
# Organ functioning (e.g., kidney(creatinine, blood urea nitrogen), liver(AST, ALT))
# Metabolic markers (fasting LDL, HDL, total cholesterol, Triglycerides, fasting glucose)
subdata <- data[ , names(data) %in% c("publicID", "study", "week", "specage", "sgptval", "sgotval", "bunval", 
                 "cretval", "ftcval", "fldlval", "fhdlval", "ftgval", "fbsval") ]
# average starting time 12.03 year, Min: 6.87, Max: 16.38
summary(subdata[which(subdata$week==0),  ]$specage )
# Cardiovascular indicator(blood pressure) 
bloodpressure <- read_sas("evw0173.sas7bdat")
## the above biological measures have at least 3 data points


#not enough longitudinal
GlucoseToleranceTest <- read_sas("lbw0074.sas7bdat")
GlucoseToleranceTest$hba1cval
#"hba1cval" in the data set LBW0074