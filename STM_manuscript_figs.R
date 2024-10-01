###R code for manuscript: Targeting mosquito hydroxyphenylpyruvate dioxygenase (HPPD) as a novel strategy for preventing malaria transmission

##Erase Everything and Set the working directory to R file location
rm(list=ls())

### install required packages
list.of.packages <- c("gridExtra", "grid", "ggplot2", "ggpubr", "patchwork", "readr", "plotrix", "rstudioapi", "tidyr", "reshape2", "dplyr", "ggtext", "survival", "survminer", "readxl", "stringr", "plyr", "RColorBrewer","GGally", "unikn", "viridis")
###Check if any of these packages is not installed and add those not already installed to object new.packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
##Install the packages not already installed
if(length(new.packages)) install.packages(new.packages)
##Load all the packages including the installed ones
lapply(list.of.packages, require, character.only = TRUE)

#Set the working directory to file location
setwd(dirname(getActiveDocumentContext()$path))

#Set colour themes for figures
IVMcolor= "#689ed4"
NTBCcolor=    "#c07142"
bgcolor="#fdf8e7"
legendcolor="#f4eed4"
legendline="grey93"


setwd("/Users/annatrett/Box/STM Manuscript Sept 2024/Data")


PlotTHEME=theme(axis.text.x = element_text(family = "Myriad Pro", size = 9, colour="black"),
                axis.text.y = element_text(family = "Myriad Pro",size=9, colour="black"),
                axis.title.y = element_text( family = "Myriad Pro",size=9, vjust= 1),
                axis.title.x = element_text( family = "Myriad Pro",size=9, vjust = 1),
                text=element_text(family= "Myriad Pro Bold"),
                legend.title = element_blank(),legend.spacing.y = unit(1, 'mm'),
                axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
                plot.margin = margin(5,5,5,5),
                
                legend.direction = "horizontal",
                legend.position = "top",
                legend.box.background=element_rect(fill=legendcolor, size=0.5),
                legend.background= element_rect(fill=legendcolor),
                legend.text = element_text(family = "Myriad Pro",size=7),
                plot.title = element_text(family = "Myriad Pro",size=9),
                legend.key.height =unit(0.35, 'cm'),
                legend.margin = margin(0.06,0.06,0.06,0.06, unit="cm"),
                panel.background = element_rect(fill=bgcolor),
                panel.grid.major = element_line(colour = legendline, linetype=1, size=0.25),
                panel.grid.minor = element_blank(),
                plot.tag = element_text(family = "Myriad Pro Bold", size = 9, face= "bold", hjust=0, vjust=0),
                legend.key = element_rect(fill = alpha("white", 0.0)),
                legend.key.width = unit(0.3, "cm"))


#Set working directory

#Set function to convert survival data into Kaplan ready analysis data
Kaplan_Func = function (kaplan) {
  
  
  kaplan=kaplan[c(1,3:nrow(kaplan)),]
  
  kaplan[2,]=kaplan[1,]
  kaplan=kaplan[-1,]
  
  
  names(kaplan)[1]="Time"
  kaplan$Time = as.numeric(kaplan$Time)
  kaplan[1,1]=0
  
  kaplan$`...11`=NULL
  
  
  kaplanS=kaplan
  
  for(i in 2:nrow(kaplan)){
    
    for (j in 2:length(kaplan)){
      kaplanS[i,j]=kaplan[1,j]-kaplan[i,j]
    }
  }
  kaplanS=kaplanS[1:15,]
  
  #Melt the data into a column format (easier to work with for R)
  MeltedKaplan<-reshape2::melt(data = kaplanS, id.vars = 1)
  
  ###split the data into elements representing a treatment each
  KaplanList<-split(MeltedKaplan, f = as.factor(MeltedKaplan$variable))
  
  x<-KaplanList[[2]]
  ##Function to convert data to COX/Kaplan Meyer Format
  COXIT <-function (x) { ##you can text the function by defining x as x<-KaplanList[[1]] - so then you can run it step by step and see what it's doing
    ##Step 1 - Change variable column from factor to character
    x$variable<-as.character(x$variable)
    
    ##Step -2 Create a column to count events for each time point (The reduction of mosquitoes from the last measurement at each time point)
    x$Events<-NA
    for (i in (1:nrow(x))) {
      x[i,4]<-ifelse(i==1, 0, x[i-1,3]-x[i,3])
    }
    
    ##Add a row for remaining non-events ### This is basically the number of mosquitoes that remain alive at the end of the experiment. The last value
    x[nrow(x)+1,]<-c(x[nrow(x),1] ,unique(x$variable),x[nrow(x),3], x[nrow(x),3])
    
    
    ##Add a column which defines whether Event is life or death (0 or 1)
    x$Nature<-c(rep(0, nrow(x)-1),1)
    
    
    ##Change Time, value and Events columns to numeric
    x$Time<-as.numeric(x$Time)
    x$value<-as.numeric(x$value)
    x$Events<-as.numeric(x$Events)
    
    ##Remove time points not containing events. They're not needed for a COX regression
    x<-subset(x, !x$Events==0)
    
    ##Repeat Times and Events by a factor of number events. If you have 3 mosquitoes that die at time 12 you want three entries for time 12, if you have 10 more mosquitoes that die at time 24 you want 10 entries for time 24 and so on
    Times<-rep(x$Time, x$Events)
    
    ##Same as for times. Events are coded with 0. We want to record total number of events (deaths) at the last time point and then code all reamining living mosquiotes, quantified by the last added row in the events column with   1
    
    if (1 %in% x$Nature) { ##Check if there are any living mosquitoes left at the end of the experiment , if there are then rep 0 for all events and 1 for remaining living mosquitoes
      Events<-rep(c(rep(0,(nrow(x)-1)),1), x$Events)
    } else {  #if no living mosquitoes remaining use this: Only repeat 0 (for events) for all, no 1s because no living mosquitoes
      Events<-rep(rep(0,nrow(x)), x$Events)
    }
    ##0= Event of death, 1=censored at end of experiment/no event
    variable<-rep(unique(x$variable), sum(x$Events)) ##repeat the variable to insert with the data
    
    ##Create the data
    COX<-data.frame("Times"=Times, "Events"=Events, "Treatment"= variable)
    return(COX)
  }
  
  
  #Apply the function on the new list
  newlist<-lapply(KaplanList, FUN= COXIT)
  
  ##
  COXdata<-bind_rows(newlist)
  names(COXdata)
  COXdata$Treatment<-as.factor(COXdata$Treatment)
  
  
  #Switch 0s and 1s
  COXdata$Events<-ifelse(COXdata$Events==1,0,1)
  
  
  return (COXdata)
  
}
#Same function as above but includes the CMax column, not in the other experiments
Kaplan_Func_max = function (kaplan) {
  
  kaplan=kaplan[c(1,3:nrow(kaplan)),]
  
  kaplan[2,]=kaplan[1,]
  kaplan=kaplan[-1,]
  
  
  # names(kaplan)=gsub(x=names(kaplan), pattern=" ", replacement="")
  names(kaplan)[1]="Time"
  kaplan$Time = as.numeric(kaplan$Time)
  kaplan[1,1]=0
  
  kaplan$`...11`=NULL
  
  
  
  kaplanS=kaplan
  
  for(i in 2:nrow(kaplan)){
    
    for (j in 2:length(kaplan)){
      kaplanS[i,j]=kaplan[1,j]-kaplan[i,j]
    }
  }
  kaplanS=kaplanS[1:15,]
  
  #Melt the data into a column format (easier to work with for R)
  MeltedKaplan<-reshape2::melt(data = kaplanS, id.vars = 1)
  
  ###split the data into elements representing a treatment each
  KaplanList<-split(MeltedKaplan, f = as.factor(MeltedKaplan$variable))
  
  x<-KaplanList[[2]]
  ##Function to convert data to COX/Kaplan Meyer Format
  COXIT <-function (x) { ##you can text the function by defining x as x<-KaplanList[[1]] - so then you can run it step by step and see what it's doing
    ##Step 1 - Change variable column from factor to character
    x$variable<-as.character(x$variable)
    
    ##Step -2 Create a column to count events for each time point (The reduction of mosquitoes from the last measurement at each time point)
    x$Events<-NA
    for (i in (1:nrow(x))) {
      x[i,4]<-ifelse(i==1, 0, x[i-1,3]-x[i,3])
    }
    
    ##Add a row for remaining non-events ### This is basically the number of mosquitoes that remain alive at the end of the experiment. The last value
    x[nrow(x)+1,]<-c(x[nrow(x),1] ,unique(x$variable),x[nrow(x),3], x[nrow(x),3])
    
    
    ##Add a column which defines whether Event is life or death (0 or 1)
    x$Nature<-c(rep(0, nrow(x)-1),1)
    
    
    ##Change Time, value and Events columns to numeric
    x$Time<-as.numeric(x$Time)
    x$value<-as.numeric(x$value)
    x$Events<-as.numeric(x$Events)
    
    ##Remove time points not containing events. They're not needed for a COX regression
    x<-subset(x, !x$Events==0)
    ##Repeat Times and Events by a factor of number events. If you have 3 mosquitoes that die at time 12 you want three entries for time 12, if you have 10 more mosquitoes that die at time 24 you want 10 entries for time 24 and so on
    Times<-rep(x$Time, x$Events)
    
    ##Same as for times. Events are coded with 0. We want to record total number of events (deaths) at the last time point and then code all reamining living mosquiotes, quantified by the last added row in the events column with   1
    
    if (1 %in% x$Nature) { ##Check if there are any living mosquitoes left at the end of the experiment , if there are then rep 0 for all events and 1 for remaining living mosquitoes
      Events<-rep(c(rep(0,(nrow(x)-1)),1), x$Events)
    } else {  #if no living mosquitoes remaining use this: Only repeat 0 (for events) for all, no 1s because no living mosquitoes
      Events<-rep(rep(0,nrow(x)), x$Events)
    }
    ##0= Event of death, 1=censored at end of experiment/no event
    variable<-rep(unique(x$variable), sum(x$Events)) ##repeat the variable to insert with the data
    
    ##Create the data
    COX<-data.frame("Times"=Times, "Events"=Events, "Treatment"= variable)
    return(COX)
  }
  
  
  #Apply the function on the new list
  newlist<-lapply(KaplanList, FUN= COXIT)
  
  ##
  COXdata<-bind_rows(newlist)
  names(COXdata)
  COXdata$Treatment<-as.factor(COXdata$Treatment)
  
  #Switch 0s and 1s
  COXdata$Events<-ifelse(COXdata$Events==1,0,1)
  
  
  return (COXdata)
  
}



#Read the data from excel

#Experiment 1 data, (Experiment 1 of young mosquitoes experiment)
#This is already in Kaplan ready format
COX_1 = read_csv("00-OriginalData/NTBC_IVM_Exp1.csv")
COX_1$Treatment = str_replace_all(COX_1$Treatment, "Ethanol.Control", "Control 0"	)

COX_1$Treatment=revalue(COX_1$Treatment, c("NTBC_10.000.ng.ul" ="NTBC 10000","IVM_1000.ng.ul"="IVM 1000","IVM_5000.ng.ul"="IVM 5000",
                                           "IVM_100.ng.ul"="IVM 100",  "NTBC_250.ng.ul" = "NTBC 250" ,  "IVM_15.ng.ul"="IVM 15",
                                           "NTBC_100.ng.ul"="NTBC 100", "NTBC_50.ng.ul"="NTBC 50"))


#Experiment 2 data, (Experiment 2 of young mosquitoes experiment)

exp_2_young<-read_excel(path="00-OriginalData/NTBC_IVM_Exp2_rep1.xlsx",
                        sheet = "Sheet1", range = "A2:K25", col_names = TRUE)

names(exp_2_young)=c("Time",  "Control 0",  "IVM 125", "IVM 20","IVM 8.5" , "IVM 3.5", "IVM 1.5",
                     "NTBC 3500", "NTBC 350", "NTBC 37", "NTBC 3.8")


exp_22_young<-read_excel(path="00-OriginalData/NTBC_IVM_Exp2_rep2.xlsx",
                         sheet = "Sheet1", range = "A2:K19", col_names = TRUE)

names(exp_22_young)=c("Time",   "IVM 125", "IVM 20","IVM 8.5" , "IVM 3.5", "IVM 1.5",
                      "NTBC 3500", "NTBC 350", "NTBC 37", "NTBC 3.8",  "Control 0")


exp_23_young<-read_excel(path="00-OriginalData/NTBC_IVM_Exp2_rep3.xlsx",
                         sheet = "Sheet1", range = "A2:K19", col_names = TRUE)

names(exp_23_young)=c("Time",   "IVM 125", "IVM 20","IVM 8.5" , "IVM 3.5", "IVM 1.5",
                      "NTBC 3500", "NTBC 350", "NTBC 37", "NTBC 3.8",  "Control 0")


#Convert to Kaplan format
COX_2 = Kaplan_Func_max (exp_2_young)
COX_22 = Kaplan_Func_max (exp_22_young)
COX_23 = Kaplan_Func_max (exp_23_young)

COX_2$Times = COX_2$Times/24
COX_22$Times = COX_22$Times/24
COX_23$Times = COX_23$Times/24

COX_2_cont = COX_2 %>%
  filter(Treatment== "Control 0")
COX_22_cont = COX_22 %>%
  filter(Treatment== "Control 0")
COX_23_cont = COX_23 %>%
  filter(Treatment== "Control 0")

COX_2_drug= COX_2 %>%
  filter(Treatment!= "Control 0")
COX_22_drug = COX_22 %>%
  filter(Treatment!= "Control 0")
COX_23_drug = COX_23 %>%
  filter(Treatment!= "Control 0")

#Experiment 3 data, (Experiment 3 of young mosquitoes experiment, only including extra concentrations of NTBC)

conc_exp = read_excel(path= "00-OriginalData/NTBC_Exp3.xlsx",
                      sheet = "Sheet1", range = "A2:G30", col_names = TRUE)


names(conc_exp)=c("Time",  "Control 0",   "NTBC 250", "NTBC 200", "NTBC 150", "NTBC 100",
                  "IVM 15")

conc_COX = Kaplan_Func (conc_exp)
COX_3 = conc_COX


######### Figure 1 ########## 
#Survival curve - for Experiment 1 and 2
Paste_6 = c("black", "#FF6666", "#CC9933", "#00CC33" , "#00CCFF", "#FF66FF")
Paste_5 = c("black", "#FF6666", "#CC9933", "#00CCFF", "#FF66FF")
Paste_4 = c("black", "#FF6666", "#CC9933","#00CC33")


#Survival curve plot function
Survival_curv_Func = function (data, legend.labs, Title, paste, lines) {
  data$surv = with(data, Surv(Times, Events == 1))
  
  curve <- survfit(surv ~ Treatment, data = data, conf.type = "log-log")
  survplot = ggsurvplot(curve, data,
                        conf.int = FALSE,
                        palette =  paste,
                        linetype = lines,
                        size=0.7,
                        title = Title ,
                        # font.title = c(9, "bold"),
                        
                        pval = TRUE,
                        pval.size = 3,
                        pval.coord = c(10, 0.2),
                        ylab="Mosquito survival (%)", xlab="Time (Days)",
                        risk.table = FALSE,
                        risk.table.title = "",
                        
                        xlim=c(0,14),
                        legend.title = "Conc. [ng/mL]",
                        risk.table.height = 0.2,
                        
                        #legend = c(0.75,0.8),
                        ggtheme = PlotTHEME,
                        break.time.by = 4,
                        
                        legend.labs = legend.labs) +
    guides(colour = guide_legend(nrow=1))
  
  
  
  return(survplot)
  
}
Survival_curv_Func_1_4 = function (data, legend.labs, Title, paste, lines) {
  data$surv = with(data, Surv(Times, Events == 1))
  
  curve <- survfit(surv ~ Treatment, data = data, conf.type = "log-log")
  survplot = ggsurvplot(curve, data,
                        conf.int = FALSE,
                        palette =  paste,
                        linetype = lines,
                        size=0.7,
                        title = Title ,
                        # font.title = c(9, "bold"),
                        
                        pval = TRUE,
                        pval.size = 3,
                        pval.coord = c(10, 0.2),
                        ylab="Mosquito survival (%)", xlab="",
                        risk.table = FALSE,
                        risk.table.title = "",
                        
                        xlim=c(0,14),
                        legend.title = "Conc. [ng/mL]",
                        risk.table.height = 0.2,
                        
                        #legend = c(0.75,0.8),
                        ggtheme = PlotTHEME,
                        break.time.by = 4,
                        
                        legend.labs = legend.labs) +
    guides(colour = guide_legend(nrow=1))
  
  
  
  return(survplot)
  
}

Survival_curv_Func_6 = function (data, legend.labs, Title, paste, lines) {
  data$surv = with(data, Surv(Times, Events == 1))
  
  curve <- survfit(surv ~ Treatment, data = data, conf.type = "log-log")
  survplot = ggsurvplot(curve, data,
                        conf.int = FALSE,
                        palette =  paste,
                        linetype = lines,
                        size=0.7,
                        title = Title ,
                        # font.title = c(9, "bold"),
                        
                        pval = TRUE,
                        pval.size = 3,
                        pval.coord = c(10, 0.2),
                        ylab="", xlab="Time (Days)",
                        risk.table = FALSE,
                        risk.table.title = "",
                        
                        xlim=c(0,14),
                        legend.title = "Conc. [ng/mL]",
                        risk.table.height = 0.2,
                        
                        #legend = c(0.75,0.8),
                        ggtheme = PlotTHEME,
                        break.time.by = 4,
                        
                        legend.labs = legend.labs) +
    guides(colour = guide_legend(nrow=1))
  
  
  
  return(survplot)
  
}

Survival_curv_Func_2_4 = function (data, legend.labs, Title, paste, lines) {
  data$surv = with(data, Surv(Times, Events == 1))
  
  curve <- survfit(surv ~ Treatment, data = data, conf.type = "log-log")
  survplot = ggsurvplot(curve, data,
                        conf.int = FALSE,
                        palette =  paste,
                        linetype = lines,
                        size=0.7,
                        title = Title ,
                        # font.title = c(9, "bold"),
                        
                        pval = TRUE,
                        pval.size = 3,
                        pval.coord = c(10, 0.2),
                        ylab="", xlab="",
                        risk.table = FALSE,
                        risk.table.title = "",
                        
                        xlim=c(0,14),
                        legend.title = "Conc. [ng/mL]",
                        risk.table.height = 0.2,
                        
                        #legend = c(0.75,0.8),
                        ggtheme = PlotTHEME,
                        break.time.by = 4,
                        
                        legend.labs = legend.labs) +
    guides(colour = guide_legend(nrow=1))
  
  
  
  return(survplot)
  
}



Conc_Func_IVM = function (NTBC_IV, Control) {
  
  NTBC_IV$Treatment = as.character(NTBC_IV$Treatment)
  str(NTBC_IV)
  NTBC_IV$Drug <- sapply(strsplit(NTBC_IV$Treatment, " "), "[", 1)
  
  NTBC_IV$Conc <- sapply(strsplit(NTBC_IV$Treatment, " "), "[", 2)
  NTBC_IV = NTBC_IV [ ,-3]
  colnames(NTBC_IV)[4] <- "Treatment"
  NTBC_IV$Treatment = as.character(NTBC_IV$Treatment)
  
  IVM_1=subset(NTBC_IV, NTBC_IV$Drug=="IVM")
  NTBC_1=subset(NTBC_IV, NTBC_IV$Drug=="NTBC")
  
  NTBC_1 = NTBC_1 [ ,-3]
  IVM_1 = IVM_1 [ ,-3]
  
  
  Control$Treatment = "Control"
  
  IVM_1 = rbind (Control, IVM_1)
  NTBC_1 = rbind (Control, NTBC_1)
  #Control = NTBC_IV [c(1:89),]
  
  
  return(IVM_1)
  
}
Conc_Func_NTBC = function (NTBC_IV, Control) {
  
  NTBC_IV$Treatment = as.character(NTBC_IV$Treatment)
  str(NTBC_IV)
  NTBC_IV$Drug <- sapply(strsplit(NTBC_IV$Treatment, " "), "[", 1)
  
  NTBC_IV$Conc <- sapply(strsplit(NTBC_IV$Treatment, " "), "[", 2)
  NTBC_IV = NTBC_IV [ ,-3]
  colnames(NTBC_IV)[4] <- "Treatment"
  NTBC_IV$Treatment = as.character(NTBC_IV$Treatment)
  
  IVM_1=subset(NTBC_IV, NTBC_IV$Drug=="IVM")
  NTBC_1=subset(NTBC_IV, NTBC_IV$Drug=="NTBC")
  
  NTBC_1 = NTBC_1 [ ,-3]
  IVM_1 = IVM_1 [ ,-3]
  
  
  Control$Treatment = "Control"
  
  IVM_1 = rbind (Control, IVM_1)
  NTBC_1 = rbind (Control, NTBC_1)
  #Control = NTBC_IV [c(1:89),]
  
  
  return(NTBC_1)
  
}


IVM_1 = Conc_Func_IVM (COX_2_drug, COX_2_cont)
NTBC_1 = Conc_Func_NTBC (COX_2_drug, COX_2_cont)

IVM_2 = Conc_Func_IVM (COX_22_drug, COX_22_cont)
NTBC_2 = Conc_Func_NTBC (COX_22_drug, COX_22_cont)

IVM_3 = Conc_Func_IVM (COX_23_drug, COX_23_cont)
NTBC_3 = Conc_Func_NTBC (COX_23_drug, COX_23_cont)

NTBC_2$Treatment <- factor(NTBC_2$Treatment, levels = c("Control", "3.8", "37", "350", "3500"))
IVM_2$Treatment <- factor(IVM_2$Treatment, levels = c("Control",  "1.5", "3.5", "8.5", "20", "125"))

NTBC_1$Treatment <- factor(NTBC_1$Treatment, levels = c("Control", "3.8", "37", "350", "3500"))
IVM_1$Treatment <- factor(IVM_1$Treatment, levels = c("Control",  "1.5", "3.5", "8.5", "20", "125"))

NTBC_3$Treatment <- factor(NTBC_3$Treatment, levels = c("Control", "3.8", "37", "350", "3500"))
IVM_3$Treatment <- factor(IVM_3$Treatment, levels = c("Control",  "1.5", "3.5", "8.5", "20", "125"))


NTBC_curve_1 = Survival_curv_Func_1_4 (NTBC_1, legend.labs = c("Control", "3.8", "37",  "350", "3500"), Title =   "NTBC - Replicate 1", paste = Paste_5, lines=c(1,1,1,1,2))
NTBC_curve_1$plot <-NTBC_curve_1$plot + labs(tag='A');

IVM_curve_1 = Survival_curv_Func_2_4 (IVM_1, legend.labs = c("Control", "1.5", "3.5","8.5",  "20", "125"), Title =    "IVM - Replicate 1", paste = Paste_6, lines=c(1,1,1,1,1,2) )
IVM_curve_1$plot <-IVM_curve_1$plot + labs(tag='B');

NTBC_curve_2 = Survival_curv_Func_1_4 (NTBC_2, legend.labs = c("Control", "3.8", "37",  "350", "3500"), Title =    "NTBC - Replicate 2", paste = Paste_5, lines=c(1,1,1,1,2))
NTBC_curve_2$plot <-NTBC_curve_2$plot + labs(tag='C');

IVM_curve_2 = Survival_curv_Func_2_4 (IVM_2, legend.labs = c("Control", "1.5", "3.5","8.5",  "20", "125"), Title =    "IVM - Replicate 2", paste = Paste_6, lines=c(1,1,1,1,1,2) )
IVM_curve_2$plot <-IVM_curve_2$plot + labs(tag='D');

NTBC_curve_3 = Survival_curv_Func (NTBC_3, legend.labs = c("Control", "3.8", "37",  "350", "3500"), Title =    "NTBC - Replicate 3", paste = Paste_5, lines=c(1,1,1,1,2))
NTBC_curve_3$plot <-NTBC_curve_3$plot + labs(tag='E');

IVM_curve_3 = Survival_curv_Func_6 (IVM_3, legend.labs = c("Control", "1.5", "3.5","8.5",  "20", "125"), Title =    "IVM - Replicate 3", paste = Paste_6, lines=c(1,1,1,1,1,2) )
IVM_curve_3$plot <-IVM_curve_3$plot + labs(tag='F');

plots = list(NTBC_curve_1, NTBC_curve_2, NTBC_curve_3, IVM_curve_1, IVM_curve_2,IVM_curve_3)

Figure_1_TIFF = arrange_ggsurvplots (plots, nrow=3, ncol=2)

tiff("MS_Fig_1.tiff", width = 7.3, height = 11, units = "in", res = 600)
print(Figure_1_TIFF)  # Use print() to ensure the plot is rendered
dev.off() 



######### Figure 2 ########## 



#Split data into drugs
NTBC_COX_Func = function (final) {
  final$Treatment = as.character(final$Treatment)
  
  final$Drug=gsub("Treatment", "", final$Treatment)
  final$Drug <- sapply(strsplit(final$Drug, " "), "[", 1)
  final$Conc <- sapply(strsplit(final$Treatment, " "), "[", 2)
  
  Control<-subset(final, Drug=="Control") ###NTBC data
  IVM<-subset(final, Drug=="IVM") ###NTBC data
  NTBC<-subset(final, Drug=="NTBC") ###NTBC data
  
  NTBC = rbind (NTBC, Control)
  IVM = rbind (IVM, Control)
  
  
  NTBC = NTBC [,-3]
  return (NTBC)
}
IVM_COX_Func = function (final) {
  final$Treatment = as.character(final$Treatment)
  
  final$Drug=gsub("Treatment", "", final$Treatment)
  final$Drug <- sapply(strsplit(final$Drug, " "), "[", 1)
  final$Conc <- sapply(strsplit(final$Treatment, " "), "[", 2)
  
  Control<-subset(final, Drug=="Control") ###NTBC data
  IVM<-subset(final, Drug=="IVM") ###NTBC data
  NTBC<-subset(final, Drug=="NTBC") ###NTBC data
  
  NTBC = rbind (NTBC, Control)
  IVM = rbind (IVM, Control)
  IVM = IVM [,-3]
  
  return (IVM)
}


NTBC_1 = NTBC_COX_Func (COX_1)
IVM_1 = IVM_COX_Func (COX_1)

IVM_1$Conc <- factor(IVM_1$Conc, levels = c("0", "15", "100", "1000", "5000"))
NTBC_1$Conc <- factor(NTBC_1$Conc, levels = c("0", "50", "100", "250", "10000"))

NTBC_2 = NTBC_COX_Func (COX_2)
IVM_2 = IVM_COX_Func (COX_2)

IVM_2$Conc <- factor(IVM_2$Conc, levels = c("0", "1.5", "3.5", "8.5", "20", "125"))
NTBC_2$Conc <- factor(NTBC_2$Conc, levels = c("0", "3.8", "37", "350", "3500"))


NTBC_22 = NTBC_COX_Func (COX_22)
IVM_22 = IVM_COX_Func (COX_22)

IVM_22$Conc <- factor(IVM_22$Conc, levels = c("0", "1.5", "3.5", "8.5", "20", "125"))
NTBC_22$Conc <- factor(NTBC_22$Conc, levels = c("0", "3.8", "37", "350", "3500"))


NTBC_23 = NTBC_COX_Func (COX_23)
IVM_23 = IVM_COX_Func (COX_23)

IVM_23$Conc <- factor(IVM_23$Conc, levels = c("0", "1.5", "3.5", "8.5", "20", "125"))
NTBC_23$Conc <- factor(NTBC_23$Conc, levels = c("0", "3.8", "37", "350", "3500"))



NTBC_3 = NTBC_COX_Func (COX_3)

NTBC_3$Conc <- factor(NTBC_3$Conc, levels = c("0", "100", "150", "200", "250"))

#Cox regressions
cox_IVM_1 <- coxph(Surv(Times, Events) ~ Conc, data = IVM_1)
cox_NTBC_1 <- coxph(Surv(Times, Events) ~ Conc, data = NTBC_1)

cox_IVM_2 <- coxph(Surv(Times, Events) ~ Conc, data = IVM_2)
cox_NTBC_2 <- coxph(Surv(Times, Events) ~ Conc, data = NTBC_2)

cox_IVM_22 <- coxph(Surv(Times, Events) ~ Conc, data = IVM_22)
cox_NTBC_22 <- coxph(Surv(Times, Events) ~ Conc, data = NTBC_22)

cox_IVM_23 <- coxph(Surv(Times, Events) ~ Conc, data = IVM_23)
cox_NTBC_23 <- coxph(Surv(Times, Events) ~ Conc, data = NTBC_23)

cox_NTBC_3 <- coxph(Surv(Times, Events) ~ Conc, data = NTBC_3)


#Extract the Cox regression out come into table format
res_COX_Func = function (res.cox) {
  
  ##extract the main findings
  final<-as.data.frame(summary(res.cox)[7])
  
  ##Add a column for treatment names taken from rownames
  final$treatment=row.names(final)
  
  ##Delete the row names
  row.names(final)=NULL
  
  #extract Treatment name, Coefficient, SE and the P value
  final=final[,c(6,2,3,5)]
  
  ##Change the names of the data frame
  names(final)=c("Treatment", "HR", "SE", "P")
  
  ##take the exp of SE
  final$SE=exp(final$SE)
  
  data=as.data.frame(summary(res.cox)$conf.int)
  
  
  final$P95=data$`upper .95`
  final$P05=data$`lower .95`
  final=final[match(str_sort(as.character(final$Treatment), numeric = TRUE),final$Treatment),]
  
  final$Conc=gsub("Conc", "", final$Treatment)
  
  
  
  #final$Treatment <- sapply(strsplit(final$Treatment, "Treatment"), "[", 2)
  
  
  
  return (final)
}


IVM_1_HR = res_COX_Func(cox_IVM_1)
NTBC_1_HR = res_COX_Func(cox_NTBC_1)

IVM_2_HR = res_COX_Func(cox_IVM_2)
NTBC_2_HR = res_COX_Func(cox_NTBC_2)

IVM_22_HR = res_COX_Func(cox_IVM_22)
NTBC_22_HR = res_COX_Func(cox_NTBC_22)

IVM_23_HR = res_COX_Func(cox_IVM_23)
NTBC_23_HR = res_COX_Func(cox_NTBC_23)

NTBC_3_HR = res_COX_Func(cox_NTBC_3)

#Combine HRs from all three experiments into one table
IVM_HRs = rbind (IVM_1_HR, IVM_2_HR,IVM_22_HR, IVM_23_HR)
NTBC_HRs = rbind (NTBC_1_HR, NTBC_2_HR, NTBC_22_HR, NTBC_23_HR, NTBC_3_HR)

IVM_HRs = as.data.frame(IVM_HRs)
NTBC_HRs = as.data.frame(NTBC_HRs)


#Concentration curves
#New + Original NTBC data - same control


detach(package:plyr)
library(dplyr)
NTBC_HRs = NTBC_HRs %>%
  group_by (Conc) %>%
  summarise (HR=mean(HR),
             SE=mean(SE),
             P=mean(P),
             P05=min(P05),
             P95=max(P95))

NTBCdata = NTBC_HRs

IVM_HRs = IVM_HRs %>%
  group_by (Conc) %>%
  summarise (HR=mean(HR),
             SE=mean(SE),
             P=mean(P),
             P05=min(P05),
             P95=max(P95))

IVMdata = IVM_HRs

IVMdata = as.data.frame(IVMdata)
NTBCdata = as.data.frame(NTBCdata)

IVMdata$HR[IVMdata$HR >= 20] <- 20
NTBCdata$HR[NTBCdata$HR >= 20] <- 20

IVMdata$Conc = as.numeric (IVMdata$Conc)
NTBCdata$Conc = as.numeric (NTBCdata$Conc)

h= 12
  var=NTBCdata$HR
  conc= c(seq(1,100000, 1))
  param4<-nls(var ~ ((19 * Conc^h)/(ec^h+Conc^h))+1, data=NTBCdata, start=list(ec=205))
  Coef4param<-as.data.frame(coef(param4))
  ec<-Coef4param[which(row.names(Coef4param)=="ec"),]
  Concens<-conc
  Effect<-((19 * Concens^h)/(ec^h+Concens^h))+1
  params_data_ntbc<-data.frame("Conc"=Concens, "Effect"=Effect)



  
  
  timepoint = "EC[50]"
  
  second_ntbc = "  = 205.31ng/mL"

  NTBC_HRs = rbind (NTBC_1_HR, NTBC_2_HR, NTBC_22_HR, NTBC_23_HR, NTBC_3_HR)
  
  NTBC_HRs = as.data.frame(NTBC_HRs)
  
  NTBC_HRs$Conc = as.numeric(NTBC_HRs$Conc)
  NTBC_HRs$HR[NTBC_HRs$HR >= 20] <- 20
  
ntbc_plot = ggplot()+
  geom_point(data=NTBC_HRs, aes(x=Conc, y=HR), size =4.5, alpha=0.6, col=NTBCcolor) +
  geom_line(data=params_data_ntbc, aes(x=Conc, y= Effect), size=1.3)   +
  geom_vline(xintercept = 200, size=0.4, linetype=2) +
  scale_y_continuous("Hazard Ratio Compared to Control", breaks=c(1,5,10,15,20))+
  
  geom_hline(yintercept = 10, size=0.4, linetype=2) +

  geom_text(aes(3000, 10.8, label = timepoint), size=3, col=NTBCcolor, parse=TRUE) +
  geom_text(aes(25000,10.8,label = second_ntbc), size=3, color = NTBCcolor)  +
  
  
    labs(title =  " NTBC Concentration curve" , tag= "A") +
  
  scale_x_log10("Predicted blood concentration (ng/mL)", breaks=c(1,10,100,1000,10000, 100000),
                labels=c(expression(paste("10"^0)), expression(paste("10"^1)),expression(paste("10"^2)),expression(paste("10"^3)),expression(paste("10"^4)),
                         expression(paste("10"^5)))) +

  
  PlotTHEME
ntbc_plot


h=3
var=IVMdata$HR
conc= c(seq(1,100000, 1))
param4<-nls(var ~ ((19 * Conc^h)/(ec^h+Conc^h))+1, data=IVMdata, start=list(ec=3))
Coef4param<-as.data.frame(coef(param4))
ec<-Coef4param[which(row.names(Coef4param)=="ec"),]
Concens<-conc
Effect<-((19 * Concens^h)/(ec^h+Concens^h))+1
params_data_ivm<-data.frame("Conc"=Concens, "Effect"=Effect)

IVM_HRs = rbind (IVM_1_HR, IVM_2_HR,IVM_22_HR, IVM_23_HR)

IVM_HRs = as.data.frame(IVM_HRs)

IVM_HRs$Conc = as.numeric(IVM_HRs$Conc)
IVM_HRs$HR[IVM_HRs$HR >= 20] <- 20

library(lubridate)


second_ivm  = " =13.43 ng/mL"
timepoint = "EC[50]"

ivm_plot = ggplot()+
  geom_point(data=IVM_HRs, aes(x=Conc, y=HR), size =4.5, alpha=0.6, col=IVMcolor) +
  geom_line(data=params_data_ivm, aes(x=Conc, y= Effect), size=1.3)   +
  geom_vline(xintercept = 13.43, size=0.4, linetype=2) +
  scale_y_continuous("", breaks=c(1,5,10,15,20))+
  
  geom_hline(yintercept = 10, size=0.4, linetype=2) +
  geom_text(aes(3500, 10.8, label = timepoint), size=3, col=IVMcolor, parse=TRUE) +
  geom_text(aes(25000,10.8,label = second_ivm), size=3, color = IVMcolor)  +
  
  labs(title = " IVM Concentration curve", tag = expression(paste(bold("B")))) +

  scale_x_log10("Predicted blood concentration (ng/mL)", breaks=c(1,10,100,1000,10000, 100000),
                labels=c(expression(paste("10"^0)), expression(paste("10"^1)),expression(paste("10"^2)),expression(paste("10"^3)),expression(paste("10"^4)),
                         expression(paste("10"^5)))) +
  
  PlotTHEME
ivm_plot

library(cowplot)


Fig_2 = plot_grid(ntbc_plot, ivm_plot)



tiff("MS_Fig_2.tiff", width = 7.3, height = 3.5, units = "in", res = 600)
print(Fig_2)  # Use print() to ensure the plot is rendered
dev.off() 


######### Figure 3 ########## 

IVM_3_600_Cc = readRDS("PKProfiles/IVM_3_600_Cc.RDS")
IVM_3_600_HR = readRDS("PKProfiles/IVM_3_600_HR.RDS")

NTBC_3_1mg_Cc = readRDS("PKProfiles/NTBC_3_1mg_Cc.RDS")
NTBC_3_1mg_HR = readRDS("PKProfiles/NTBC_3_1mg_HR.RDS")

IVM_3_300_Cc = readRDS("PKProfiles/IVM_3_300_Cc.RDS")
IVM_3_300_HR = readRDS("PKProfiles/IVM_3_300_HR.RDS")

IVM_1_400_Cc = readRDS("PKProfiles/IVM_1_400_Cc.RDS")
IVM_1_400_HR = readRDS("PKProfiles/IVM_1_400_HR.RDS")


##Plotting PK/PD profiles of specified dosing regimens

label1 ="EC[50]"
label2=  "EC[50]"

Cc_Func = function (Cc_profile_1, Cc_profile_2) {
  
  plot =  ggplot() +
    
    geom_line(data=Cc_profile_1, aes(x=TIME/24, y= clin_med, col=IVMcolor), size=0.5)  +
    geom_line(data=Cc_profile_2, aes(x=TIME/24, y= clin_med, col=NTBCcolor), size=0.5) +
    geom_hline(yintercept = 13, size=0.4, linetype=2, col=IVMcolor) +
    
    geom_hline(yintercept = 205, size=0.4, linetype=2, col=NTBCcolor) +
    
    annotate("text", x = 26.8, y = 17,  label = label1, parse = TRUE, size = 2.5, col= IVMcolor, fontface = "bold")+
    
    annotate("text", x = 26.8, y = 266,  label = label2, parse = TRUE, size = 2.5, col= NTBCcolor, fontface = "bold")+
    
    labs(title = " Predicted PK profile", tag=expression(paste(bold("A")))) +
    scale_y_log10("Predicted blood conc. (ng/ml)", breaks=c(1,10,100,1000,10000, 100000))+
    
    scale_color_manual(values=c(IVMcolor, NTBCcolor), labels=c("IVM (0.6mg/kg qd 3 d)", "NTBC (1.0 mg/kg qd 3 d)"))+
    
    scale_x_continuous("Time (Days)" , limits=c(0,28), breaks=c(0,7,14,21,28)) +
    
    PlotTHEME +
    theme(legend.position = "none")
  
  return (plot)
}

PK_plot = Cc_Func (IVM_3_600_Cc, NTBC_3_1mg_Cc)

Cc_profile_1 = IVM_3_600_Cc
Cc_profile_2 = NTBC_3_1mg_Cc

HR_Func = function (HR_profile_1, HR_profile_2) {
  
  plot =  ggplot()+
    # geom_smooth(data=dataset, aes(x=conc, y=mortality, col=red))
    geom_line(data=HR_profile_1, aes(x=TIME/24, y= clin_med, col=IVMcolor), size=0.5)   +
    geom_line(data=HR_profile_2, aes(x=TIME/24, y= clin_med,col=NTBCcolor), size=0.5)   +
    
    geom_hline(yintercept = 4, size=0.4, linetype=2) +
    
    labs(title = " Predicted PD profile", tag=expression(paste(bold("B")))) +
    
    
    scale_x_continuous("Time (Days)", breaks=c(0,7,14,21,28))+
    
    geom_text(aes(26,5.6, label = "HR=4", vjust = 2), size=3, col="black") +
    scale_color_manual(values=c(IVMcolor, NTBCcolor), labels=c("IVM (0.6mg/kg qd 3 d)", "NTBC (1.0 mg/kg qd 3 d)"))+
    
    scale_y_continuous("Predicted HR", breaks=c(0,5,10,15,20)) +
    PlotTHEME  +
    theme(legend.position = "none")
  return(plot)
}

PD_plot = HR_Func (IVM_3_600_HR, NTBC_3_1mg_HR)



Fig_3 = plot_grid(PK_plot, PD_plot)



tiff("MS_Fig_3.tiff", width = 7.3, height = 3.5, units = "in", res = 600)
print(Fig_3)  # Use print() to ensure the plot is rendered
dev.off() 


############# Figure 4 ##########



#Cox regression - Adult results

#Experiment 3 data, (Experiment 1 of adult mosquitoes experiment)


#Experiment 4 data, (Experiment 2 of adult mosquitoes experiment)
Kaplan_Func_max = function (kaplan) {
  
  kaplan = kaplan [-2, ]
  kaplan = kaplan [-c(17:20), ]
  
   kaplan =  kaplan
  kaplan=kaplan[c(1,3:nrow(kaplan)),]
  
  kaplan[2,]=kaplan[1,]
  kaplan=kaplan[-1,]
  
  
  # names(kaplan)=gsub(x=names(kaplan), pattern=" ", replacement="")
  names(kaplan)[1]="Time"
  kaplan$Time = as.numeric(kaplan$Time)
  kaplan[1,1]=0
  
  kaplan$`...11`=NULL
  
  
  
  kaplanS=kaplan
  
  for(i in 2:nrow(kaplan)){
    
    for (j in 2:length(kaplan)){
      kaplanS[i,j]=kaplan[1,j]-kaplan[i,j]
    }
  }
  kaplanS=kaplanS[1:14,]
  #Melt the data into a column format (easier to work with for R)
  MeltedKaplan<-reshape2::melt(data = kaplanS, id.vars = 1)
  
  ###split the data into elements representing a treatment each
  KaplanList<-split(MeltedKaplan, f = as.factor(MeltedKaplan$variable))
  
  x<-KaplanList[[2]]
  ##Function to convert data to COX/Kaplan Meyer Format
  COXIT <-function (x) { ##you can text the function by defining x as x<-KaplanList[[1]] - so then you can run it step by step and see what it's doing
    ##Step 1 - Change variable column from factor to character
    x$variable<-as.character(x$variable)
    
    ##Step -2 Create a column to count events for each time point (The reduction of mosquitoes from the last measurement at each time point)
    x$Events<-NA
    for (i in (1:nrow(x))) {
      x[i,4]<-ifelse(i==1, 0, x[i-1,3]-x[i,3])
    }
    
    ##Add a row for remaining non-events ### This is basically the number of mosquitoes that remain alive at the end of the experiment. The last value
    x[nrow(x)+1,]<-c(x[nrow(x),1] ,unique(x$variable),x[nrow(x),3], x[nrow(x),3])
    
    
    ##Add a column which defines whether Event is life or death (0 or 1)
    x$Nature<-c(rep(0, nrow(x)-1),1)
    
    
    ##Change Time, value and Events columns to numeric
    x$Time<-as.numeric(x$Time)
    x$value<-as.numeric(x$value)
    x$Events<-as.numeric(x$Events)
    
    ##Remove time points not containing events. They're not needed for a COX regression
    x<-subset(x, !x$Events==0)
    ##Repeat Times and Events by a factor of number events. If you have 3 mosquitoes that die at time 12 you want three entries for time 12, if you have 10 more mosquitoes that die at time 24 you want 10 entries for time 24 and so on
    Times<-rep(x$Time, x$Events)
    
    ##Same as for times. Events are coded with 0. We want to record total number of events (deaths) at the last time point and then code all reamining living mosquiotes, quantified by the last added row in the events column with   1
    
    if (1 %in% x$Nature) { ##Check if there are any living mosquitoes left at the end of the experiment , if there are then rep 0 for all events and 1 for remaining living mosquitoes
      Events<-rep(c(rep(0,(nrow(x)-1)),1), x$Events)
    } else {  #if no living mosquitoes remaining use this: Only repeat 0 (for events) for all, no 1s because no living mosquitoes
      Events<-rep(rep(0,nrow(x)), x$Events)
    }
    ##0= Event of death, 1=censored at end of experiment/no event
    variable<-rep(unique(x$variable), sum(x$Events)) ##repeat the variable to insert with the data
    
    ##Create the data
    COX<-data.frame("Times"=Times, "Events"=Events, "Treatment"= variable)
    return(COX)
  }
  
  
  #Apply the function on the new list
  newlist<-lapply(KaplanList, FUN= COXIT)
  
  ##
  COXdata<-bind_rows(newlist)
  names(COXdata)
  COXdata$Treatment<-as.factor(COXdata$Treatment)
  
  #Switch 0s and 1s
  COXdata$Events<-ifelse(COXdata$Events==1,0,1)
  
  
  return (COXdata)
  
}

exp_1_adult<-read_excel(path="Adult_experiments/Old_NTBC_IVM_Exp1_R1.xlsx",
                        sheet = "Sheet1", range = "A2:K23", col_names = TRUE)

names(exp_1_adult)=c("Time",  "Control",  "Cmax IVM", "IVM Day 7", "IVM Day 14", "IVM Day 21", "IVM Day 28",
                     "NTBC Day 7", "NTBC Day 14", "NTBC Day 21", "NTBC Day 28")

COX_1_old = Kaplan_Func_max (exp_1_adult)



exp_2_adult<-read_excel(path="Adult_experiments/Old_NTBC_IVM_Exp1_R2.xlsx",
                        sheet = "Sheet1", range = "A2:K19", col_names = TRUE)
names(exp_2_adult)=c("Time",    "Cmax IVM", "IVM Day 7", "IVM Day 14", "IVM Day 21", "IVM Day 28",
                     "NTBC Day 7", "NTBC Day 14", "NTBC Day 21", "NTBC Day 28", "Control")


COX_2_old = Kaplan_Func_max (exp_2_adult)


Kaplan_Func_max = function (kaplan) {
  kaplan=kaplan[c(1,3:nrow(kaplan)),]
  
  kaplan[2,]=kaplan[1,]
  kaplan=kaplan[-1,]
  
  
  # names(kaplan)=gsub(x=names(kaplan), pattern=" ", replacement="")
  names(kaplan)[1]="Time"
  kaplan$Time = as.numeric(kaplan$Time)
  kaplan[1,1]=0
  
  kaplan$`...11`=NULL
  
  
  
  kaplanS=kaplan
  
  for(i in 2:nrow(kaplan)){
    
    for (j in 2:length(kaplan)){
      kaplanS[i,j]=kaplan[1,j]-kaplan[i,j]
    }
  }
  kaplanS=kaplanS[1:14,]
  
  #Melt the data into a column format (easier to work with for R)
  MeltedKaplan<-reshape2::melt(data = kaplanS, id.vars = 1)
  
  ###split the data into elements representing a treatment each
  KaplanList<-split(MeltedKaplan, f = as.factor(MeltedKaplan$variable))
  
  x<-KaplanList[[2]]
  ##Function to convert data to COX/Kaplan Meyer Format
  COXIT <-function (x) { ##you can text the function by defining x as x<-KaplanList[[1]] - so then you can run it step by step and see what it's doing
    ##Step 1 - Change variable column from factor to character
    x$variable<-as.character(x$variable)
    
    ##Step -2 Create a column to count events for each time point (The reduction of mosquitoes from the last measurement at each time point)
    x$Events<-NA
    for (i in (1:nrow(x))) {
      x[i,4]<-ifelse(i==1, 0, x[i-1,3]-x[i,3])
    }
    
    ##Add a row for remaining non-events ### This is basically the number of mosquitoes that remain alive at the end of the experiment. The last value
    x[nrow(x)+1,]<-c(x[nrow(x),1] ,unique(x$variable),x[nrow(x),3], x[nrow(x),3])
    
    
    ##Add a column which defines whether Event is life or death (0 or 1)
    x$Nature<-c(rep(0, nrow(x)-1),1)
    
    
    ##Change Time, value and Events columns to numeric
    x$Time<-as.numeric(x$Time)
    x$value<-as.numeric(x$value)
    x$Events<-as.numeric(x$Events)
    
    ##Remove time points not containing events. They're not needed for a COX regression
    x<-subset(x, !x$Events==0)
    ##Repeat Times and Events by a factor of number events. If you have 3 mosquitoes that die at time 12 you want three entries for time 12, if you have 10 more mosquitoes that die at time 24 you want 10 entries for time 24 and so on
    Times<-rep(x$Time, x$Events)
    
    ##Same as for times. Events are coded with 0. We want to record total number of events (deaths) at the last time point and then code all reamining living mosquiotes, quantified by the last added row in the events column with   1
    
    if (1 %in% x$Nature) { ##Check if there are any living mosquitoes left at the end of the experiment , if there are then rep 0 for all events and 1 for remaining living mosquitoes
      Events<-rep(c(rep(0,(nrow(x)-1)),1), x$Events)
    } else {  #if no living mosquitoes remaining use this: Only repeat 0 (for events) for all, no 1s because no living mosquitoes
      Events<-rep(rep(0,nrow(x)), x$Events)
    }
    ##0= Event of death, 1=censored at end of experiment/no event
    variable<-rep(unique(x$variable), sum(x$Events)) ##repeat the variable to insert with the data
    
    ##Create the data
    COX<-data.frame("Times"=Times, "Events"=Events, "Treatment"= variable)
    return(COX)
  }
  
  
  #Apply the function on the new list
  newlist<-lapply(KaplanList, FUN= COXIT)
  
  ##
  COXdata<-bind_rows(newlist)
  names(COXdata)
  COXdata$Treatment<-as.factor(COXdata$Treatment)
  
  #Switch 0s and 1s
  COXdata$Events<-ifelse(COXdata$Events==1,0,1)
  
  
  return (COXdata)
  
}

exp_3_adult<-read_excel(path="Adult_experiments/Old_NTBC_IVM_Exp1_R3.xlsx",
                        sheet = "Sheet1", range = "A2:K18", col_names = TRUE)

names(exp_3_adult)=c("Time",    "Cmax IVM", "IVM Day 7", "IVM Day 14", "IVM Day 21", "IVM Day 28",
                     "NTBC Day 7", "NTBC Day 14", "NTBC Day 21", "NTBC Day 28", "Control")


COX_3_old = Kaplan_Func_max (exp_3_adult)


COXdata = rbind (COX_1_old, COX_2_old, COX_3_old)
COXdata

COXdata$Treatment <- factor(COXdata$Treatment, levels = c("Control", "IVM Day 28", "IVM Day 21","NTBC Day 28", "NTBC Day 21", "IVM Day 14","NTBC Day 14", "IVM Day 7" ,"Cmax IVM","NTBC Day 7"))

res.cox <- coxph(Surv(Times, Events) ~ Treatment, data = COXdata)

##extract the main findings
final<-as.data.frame(summary(res.cox)[7])

##Add a column for treatment names taken from rownames
final$treatment=row.names(final)

##Delete the row names
row.names(final)=NULL

#extract Treatment name, Coefficient, SE and the P value
final=final[,c(6,2,3,5)]

##Change the names of the data frame
names(final)=c("Treatment", "HR", "SE", "P")

##take the exp of SE
final$SE=exp(final$SE)
#read the 5th and 95th percentiles
data=as.data.frame(summary(res.cox)$conf.int)

##Add 5th and 95th percentiles to the main dataframe
final$P95=data$`upper .95`
final$P05=data$`lower .95`

##Sort the data alphebatically (this by passes the issue of double digits regarded as lower than single digits for example IVM7 is after IVM28
final=final[match(str_sort(as.character(final$Treatment), numeric = TRUE),final$Treatment),]



##Add a column for the Day
final$Day=gsub("TreatmentIVM|TreatmentNTBC", "", final$Treatment)
final$Day=gsub("Day", "Day ", final$Day)


##Add a column for the drug used
final$Drug=gsub("Treatment|Day|[0-9]", "", final$Treatment)


library(plyr)
#Convert Day 7 or Day07
final$Day=revalue(final$Day, c(" Day  7"=" Day  07"))
final$Drug=revalue(final$Drug, c("NTBC  "="NTBC"))
final$Drug=revalue(final$Drug, c("IVM  "="IVM"))

#Reorder Drug factors so NTBC comes first in plots
final$Drug=factor(final$Drug, levels=c("NTBC", "IVM"))


#Reorder by day
final=final[order(final$Day),]

COXdata = as.data.frame(COXdata)

####Figure 5 Code -Faceted surival plots, PK curve and HR comparison between drug concs/days
DrugCompare=COXdata
##Add a time point column
DrugCompare$Timepoint=gsub("IVM|NTBC", "", DrugCompare$Treatment)
##Add a drug column
DrugCompare$DRUG=gsub("Day|[0-9]", "", DrugCompare$Treatment)
#remove the control column
DrugCompare=subset(DrugCompare, !DrugCompare$Treatment=="Control")
DrugCompare=subset(DrugCompare, !DrugCompare$Treatment=="Cmax IVM")
##split into a list based on the time point
DrugList<-split(DrugCompare, f=DrugCompare$Timepoint)
try=DrugList[[2]]

##Function to compare drugs for each time point
comparedrugs<-function(try) {
  res.cox <- coxph(Surv(Times, Events) ~ DRUG, data = try)
  
  #Extract the hazard ratio and P value
  summary(res.cox)$conf.int
  
  Hazard_P=data.frame("Parameter"=c("Hazard", "P"), Value=summary(res.cox)$coefficients[c(2,5)], P05=summary(res.cox)$conf.int[1,3], P95=summary(res.cox)$conf.int[1,4])
  Hazard_P$Significant=ifelse(Hazard_P$Value[which(Hazard_P$Parameter=="P")]<0.05, "YES", "NO")
  Hazard_P$TimePoint=unique(try$Timepoint)
  return(Hazard_P)
}
library(plyr)
###Apply the function on all elements of the DrugList
stats=bind_rows(lapply(DrugList, FUN=comparedrugs))
#Convert Day 7 or Day07
stats$TimePoint=revalue(stats$TimePoint, c(" Day 7"="Day 07"))

stats=stats[order(stats$TimePoint),]

###Calculate hazard ratios between drugs for each group
#Convert times to days
COXdata$Times=COXdata$Times/24
COXdata=subset(COXdata, !COXdata$Treatment=="Cmax IVM")
##rename the treatments

##Reorder the treatments
COXdata$Treatment <- factor(COXdata$Treatment, levels = c("Control", "NTBC Day 7","IVM Day 7",  "NTBC Day 14", "IVM Day 14", "NTBC Day 21","IVM Day 21", "NTBC Day 28", "IVM Day 28" ))
# COXdata = COXdata[,-5]

#Compare the drugs at  each time point
fitall <- survfit(Surv(Times, Events) ~ Treatment, data = COXdata)
####Make a facet surv plot
COXdata$Drug=gsub(" Day |[0-9]", "", COXdata$Treatment)
COXdata$Day=gsub("NTBC|IVM", "", COXdata$Treatment)

COXdata$Day=factor(COXdata$Day, levels=c("Control", " Day 7",  " Day 14", " Day 21",  " Day 28"))

##Create a list where the data is split by Day
DayList=split(COXdata, f=COXdata$Day)

##Add the control to each element of the list
control=DayList[[1]]
addcontrol=function(x) {
  control1=control
  control1$Day=unique(x$Day)
  all=rbind(control,x)
  return(all)
}
##New list with control added
DayList1=lapply(DayList, FUN=addcontrol)


##Remove the control from Day List1
DayList1[[1]]=NULL

##Create a survfit for each element
fit_t=function(x){
  fit <- survfit(Surv(Times, Events) ~ Treatment,
                 data = x)
  return(fit)
}


final = final [-9, ]
final$Day1=as.character(rep(seq(7,28,7), each=2))
final$PP=ifelse(final$P<0.0001, "p<0.0001",
                ifelse(final$P<0.001, "p<0.001",
                       ifelse(final$P<0.01, "p<0.01",
                              ifelse(final$P<0.05, "p<0.05", paste0("p=", round(final$P,3))))))

DayList2=lapply(DayList1, FUN=fit_t)
fit=DayList2[[3]]

Panel6Theme=theme(plot.margin = margin(10,10,10,10),legend.position = c(0.75,0.8), legend.box.margin = margin(t=0.01,r=0.01,b=0.01,l=0.01, unit="mm"))
Panel4Theme=theme(plot.margin = margin(20,20,20,20),legend.position = c(0.75,0.75), legend.box.margin = margin(t=0.01,r=0.01,b=0.01,l=0.01, unit="cm"))



Panel6Theme=theme(plot.margin = margin(10,10,10,10),legend.position = c(0.74,0.81), legend.margin=margin(0.1, 0.1, 0.1, 0.1  , unit ="mm"))


Panel6Theme=theme(plot.margin = margin(10,10,10,10),legend.position = c(0.74,0.81), legend.box.margin = margin(0.001,0.001,0.001,0.001, unit="cm"))
Panel4Theme=theme(plot.margin = margin(20,20,20,20),legend.position = c(0.74,0.81), legend.box.margin = margin(t=0.01,r=0.01,b=0.01,l=0.01, unit="cm"))


HR = final$HR
P05  = final$P05
P95  = final$P95
PP  = final$P


CONC = c(20, 3500, 8.5, 350, 3.5, 37, 1.5, 3.8)

final$CONC = CONC

plottem=function (fit) {
  DayE=gsub("Treatment=|NTBC |IVM |", "", names(fit$strata[2]))
  DayN=gsub("Treatment=|NTBC|IVM|Day| ", "", names(fit$strata[2]))
  
  IVMCONC=final$CONC[which(final$Drug=="IVM"&final$Day1==DayN)]
  NTBCCONC=final$CONC[which(final$Drug=="NTBC"&final$Day1==DayN)]
  data=data.frame("IVM"=as.character(c("20.0",	"8.50",	"3.50",	"1.50")), "NTBC"=as.character(c("3500",	"360.0",	"37.0",	"4.00")), "Day"=as.character(seq(7,28,7)))
  melt(data,id.vars=c(3))
  final=merge(final, data, by.x = 9, by.y = 3)
  
  #titles for a 6_pnanel figure
  panel=ifelse(DayN=="7", "B ", ifelse(DayN=="14", "C ", ifelse(DayN=="21", "D ", ifelse(DayN=="28", "E "))))
  #Titles for a 4-panel figure
  panel_4=ifelse(DayN=="7", "(a) ", ifelse(DayN=="14", "(b) ", ifelse(DayN=="21", "(c) ", ifelse(DayN=="28", "(d) "))))
  
  

  plot=ggsurv(fit,
              
              size.est=0.5, ) +
    scale_color_manual(values=c("grey30", IVMcolor, NTBCcolor),#this line changed.
                       labels=c("Control  " ,
                                paste0("IVM (", IVMCONC,   " ng/mL)"),
                                paste0("NTBC (", NTBCCONC,   " ng/mL)")))+
    
  
    
    scale_x_continuous("Time [Days]", limits=c(0,14.9), breaks=seq(0,14.9,2)) +
    coord_cartesian(xlim = c(0,14.1)) +
    scale_y_continuous("Mosquito survival", limits=c(0,1)) +
    

    PlotTHEME + 
    theme(legend.position = c(0.82,0.87),
          legend.direction = "vertical") +
    guides (linetype = FALSE) 
  
  plot4panel=plot+Panel4Theme+labs(title = substitute(paste(bold(a), "Mosquitocidal activity: ", b), list(a=panel_4,b=DayE)))

  plot6panel=plot+labs(title = substitute( paste( "Mosquitocidal activity: ", b), list(b=DayE)), tag = substitute(paste(bold(a)), list (a=panel)))
                                           
                                           
                                           
  #     tag = substitute(paste(bold(a)), list (a=panel))))
  
  

  
  return(list(plot4panel,plot6panel))
}
AllPlots=lapply(DayList2, FUN=plottem)
dev.off()
##Generate lists for the 4 panel and 6 panel figures
AllPlots4 <- lapply(AllPlots, `[[`, 1)
AllPlots6 <- lapply(AllPlots, `[[`, 2)
AllPlots=lapply(DayList2, FUN=plottem)

##Generate lists for the 4 panel and 6 panel figures
AllPlots6 <- lapply(AllPlots, `[[`, 2)
fitall <- survfit(Surv(Times, Events) ~ Treatment, data = COXdata)
dev.off()
dev.list()
  

AllPlots6[[1]]


MED_NTBC = plyr::ddply(.data = NTBC_3_1mg_Cc, .variables = "TIME", summarise, "clin_med" = median(clin_med)) #take the median HR profile
# MED_IVM_300 = plyr::ddply(.data = IVM_3_300_Cc, .variables = "TIME", summarise, "clin_med" = median(clin_med))#take the median HR profile
MED_IVM_600 = plyr::ddply(.data = IVM_3_600_Cc, .variables = "TIME", summarise, "clin_med" = median(clin_med)) #take the median HR profile


MED_N = MED_NTBC [c(7*24, 14*24, 21*24, 28*24), ]
MED_I = MED_IVM_600 [c(2.3*24,7*24, 14*24, 21*24, 28*24), ]
#  MED_I = MED_IVM_300 [c(7*24, 14*24, 21*24, 28*24), ]

# plot_PK =

MED_I$med = MED_I$clin_med+0.6
MED_IVM_600$med = MED_IVM_600$clin_med+0.6
# MED_IVM_300$med = MED_IVM_300$clin_med+0.6


plot_PK =
  ggplot() +
  geom_line(data=MED_NTBC, aes(TIME/24, clin_med, col=IVMcolor),  size=1) +
  geom_line(data=MED_IVM_600, aes(TIME/24, med, col=NTBCcolor),  size=1) +
  geom_point(data=MED_I, aes(TIME/24, med, col=NTBCcolor),  size=4.5)+
  geom_point(data=MED_N, aes(TIME/24, clin_med, col=IVMcolor), size=4.5)+
  
  scale_y_log10("Predicted plasma conc. [ng/mL]", limits = c(1,20000)) +
  scale_x_continuous("Time [Days]", breaks=seq(0,31,7))+
  scale_color_manual(values=c( NTBCcolor, IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_fill_manual(values=c(NTBCcolor,IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_linetype_manual(values=c(1,1),            labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  guides(colour = guide_legend(override.aes = list(fatten = 15))) +
  
  labs(title =  " Predicted PK profiles", tag = expression(paste(bold("A")))) +
  

  annotate("text", x = 10, y= 6000 ,  colour = NTBCcolor, label="3467.6 ng/mL", size=3)+
  annotate("text", x = 16.5, y= 780 ,  colour = NTBCcolor, label="358.9 ng/mL", size=3)+
  annotate("text", x = 23.5, y= 70 ,  colour = NTBCcolor, label="37.1 ng/mL", size=3)+
  annotate("text", x = 27.2, y= 14.5 ,  colour = NTBCcolor, label="3.8 ng/mL", size=3)+
  
  annotate("text", x = 3.4, y= 220 ,  colour = IVMcolor, label="125 ng/mL", size=3)+
  annotate("text", x = 9, y= 39.75 ,  colour = IVMcolor, label="20 ng/mL", size=3)+
  annotate("text", x = 16.3, y= 11.5 ,  colour = IVMcolor, label="8.3 ng/mL", size=3)+
  annotate("text", x = 23, y= 4.5 ,  colour = IVMcolor, label="3.4 ng/mL", size=3)+
  annotate("text", x = 27.2, y= 2.1 ,  colour = IVMcolor, label="1.4 ng/mL", size=3) +
  
  
  PlotTHEME +
  theme(legend.position = c(0.75,0.87),
        legend.direction = "vertical")



plot_PK



dev.off()
comp<-ggplot(final, aes(Day, HR, group=as.character(Drug)))+
  
  geom_pointrange(
    aes(ymin = HR, ymax = HR, color = Drug), size=1.5,
    position = position_dodge(0.3), show.legend=FALSE
  )+
  geom_point(aes(color=Drug), x=1, y=500, size=4.5)+
  geom_line(aes(color=Drug, linetype=Drug), position=position_dodge(width=0.3), size=1.5)+
  
  geom_errorbar(aes(ymin = P05, ymax = P95, color=Drug),  width = 0.25,  position=position_dodge(width=0.3))+
  
  
  geom_hline(yintercept = 1, lty=2)+
  
  scale_y_log10("Mosquitocidal activity [HR]", limits = c(NA,100)) +
  scale_x_discrete("Predicted drug concentration") +
  scale_color_manual(values=c(NTBCcolor, IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_fill_manual(values=c(NTBCcolor, IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_linetype_manual(values=c(1,1),            labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  guides(colour = guide_legend(override.aes = list(fatten = 15)))+
  labs(title =  " Hazard ratio comparison", tag = expression(paste(bold("F")))) +
  
PlotTHEME +
  theme(legend.position = c(0.75,0.87),
        legend.direction = "vertical")


AllPlots[[9]]=comp

AllPlots[[8]] = plot_PK
dev.off()

panel6final=ggarrange(AllPlots[[8]], AllPlots6[[1]], AllPlots6[[2]], AllPlots6[[3]], AllPlots6[[4]], AllPlots[[9]], nrow = 2, ncol=3, align="hv")


tiff("MS_Fig_4_new.tiff", width = 14, height = 7.3, units = "in", res = 600)
print(panel6final )  # Use print() to ensure the plot is rendered
dev.off() 



############Figure 5  ############
#- kisumu vs tiassale strain

Tiassale_one<-read_excel(path="Tiassale_experiments/Tia_NTBC_Exp1_rep1.xlsx",
                     sheet = "Sheet1", range = "A2:K18", col_names = TRUE)


Tiassale_one = Tiassale_one [-2 , ]


Tiassale_two<-read_excel(path="Tiassale_experiments/Tia_NTBC_Exp1_rep2.xlsx",
                         sheet = "Sheet1", range = "A2:K18", col_names = TRUE)


Tiassale_two = Tiassale_two [-2 , ]

Tiassale_three<-read_excel(path="Tiassale_experiments/Tia_NTBC_Exp1_rep3.xlsx",
                         sheet = "Sheet1", range = "A2:K19", col_names = TRUE)

Tiassale_three = Tiassale_three [-2 , ]


kaplan_IR_Func = function (data) {
  

kaplan=data[c(1,3:nrow(data)),]

kaplan[2,]=kaplan[1,]
kaplan=kaplan[-1,]


names(kaplan)=gsub(x=names(kaplan), pattern=" ", replacement="")
names(kaplan)[1]="Time"
kaplan$Time = as.numeric(kaplan$Time)
kaplan[1,1]=0

kaplan$`...11`=NULL

return(kaplan) }



Tiassale_one_kap = kaplan_IR_Func (Tiassale_one)
Tiassale_one_kap = Tiassale_one_kap [-14, ]

Tiassale_two_kap = kaplan_IR_Func (Tiassale_two)
Tiassale_two_kap = Tiassale_two_kap [-14, ]
Tiassale_three_kap = kaplan_IR_Func (Tiassale_three)

Kisumu_split_Func = function(kaplan) {

  
kisumu=kaplan [,c(1:6)]
Tiassale=kaplan [,c(1,7,8,9,10,11)]

names(Tiassale)=c("Time",  "3500", "350", "37", "3.8","Control")
Tiassale$Time = as.numeric(Tiassale$Time)
names(kisumu)=c("Time",  "3500", "350", "37", "3.8","Control")
kisumu$Time = as.numeric(kisumu$Time)

return(kisumu) }


Tiassale_split_Func = function(kaplan) {
  
  kisumu=kaplan [,c(1:6)]
  Tiassale=kaplan [,c(1,7,8,9,10,11)]
  
  names(Tiassale)=c("Time",  "3500", "350", "37", "3.8","Control")
  Tiassale$Time = as.numeric(Tiassale$Time)
  names(kisumu)=c("Time",  "3500", "350", "37", "3.8","Control")
  kisumu$Time = as.numeric(kisumu$Time)
  
  return(Tiassale) }

Kisumu_one = Kisumu_split_Func (Tiassale_one_kap)
Tiassale_one = Tiassale_split_Func (Tiassale_one_kap)
Kisumu_two = Kisumu_split_Func (Tiassale_two_kap)
Tiassale_two = Tiassale_split_Func (Tiassale_two_kap)
Kisumu_three = Kisumu_split_Func (Tiassale_three_kap)
Tiassale_three = Tiassale_split_Func (Tiassale_three_kap)

kaplan_Tiassale_Func = function (kaplan, kaplanS) {
  

  for(i in 2:nrow(kaplan)){
    
    for (j in 2:length(kaplan)){
      kaplanS[i,j]=kaplan[1,j]-kaplan[i,j]
    }
  }
  
  
  
  #Melt the data into a column format (easier to work with for R)
  MeltedKaplan<-reshape2::melt(data = kaplanS, id.vars = 1)
  
  ###split the data into elements representing a treatment each
  KaplanList<-split(MeltedKaplan, f = as.factor(MeltedKaplan$variable))
  
  x<-KaplanList[[2]]
  ##Function to convert data to COX/Kaplan Meyer Format
  ##Function to convert data to COX/Kaplan Meyer Format
  COXIT <-function (x) { ##you can text the function by defining x as x<-KaplanList[[1]] - so then you can run it step by step and see what it's doing
    ##Step 1 - Change variable column from factor to character
    x$variable<-as.character(x$variable)
    
    ##Step -2 Create a column to count events for each time point (The reduction of mosquitoes from the last measurement at each time point)
    x$Events<-NA
    for (i in (1:nrow(x))) {
      x[i,4]<-ifelse(i==1, 0, x[i-1,3]-x[i,3])
    }
    
    ##Add a row for remaining non-events ### This is basically the number of mosquitoes that remain alive at the end of the experiment. The last value
    x[nrow(x)+1,]<-c(x[nrow(x),1] ,unique(x$variable),x[nrow(x),3], x[nrow(x),3])
    
    
    ##Add a column which defines whether Event is life or death (0 or 1)
    x$Nature<-c(rep(0, nrow(x)-1),1)
    
    
    ##Change Time, value and Events columns to numeric
    x$Time<-as.numeric(x$Time)
    x$value<-as.numeric(x$value)
    x$Events<-as.numeric(x$Events)
    
    ##Remove time points not containing events. They're not needed for a COX regression
    x<-subset(x, !x$Events==0)
    
    ##Repeat Times and Events by a factor of number events. If you have 3 mosquitoes that die at time 12 you want three entries for time 12, if you have 10 more mosquitoes that die at time 24 you want 10 entries for time 24 and so on
    Times<-rep(x$Time, x$Events)
    
    ##Same as for times. Events are coded with 0. We want to record total number of events (deaths) at the last time point and then code all reamining living mosquiotes, quantified by the last added row in the events column with   1
    
    if (1 %in% x$Nature) { ##Check if there are any living mosquitoes left at the end of the experiment , if there are then rep 0 for all events and 1 for remaining living mosquitoes
      Events<-rep(c(rep(0,(nrow(x)-1)),1), x$Events)
    } else {  #if no living mosquitoes remaining use this: Only repeat 0 (for events) for all, no 1s because no living mosquitoes
      Events<-rep(rep(0,nrow(x)), x$Events)
    }
    ##0= Event of death, 1=censored at end of experiment/no event
    variable<-rep(unique(x$variable), sum(x$Events)) ##repeat the variable to insert with the data
    
    ##Create the data
    COX<-data.frame("Times"=Times, "Events"=Events, "Treatment"= variable)
    
    COX$Times = COX$Times/24
    return(COX)
  }
  
  
  #Apply the function on the new list
  newlist<-lapply(KaplanList, FUN= COXIT)
  
  ##
  COXdata<-bind_rows(newlist)
  names(COXdata)
  COXdata$Treatment<-as.factor(COXdata$Treatment)
  ##rename the treatments"CmaxIVM" = "Cmax IVM",
  COXdata$Treatment=revalue(COXdata$Treatment, c("Control"="Control"))
  
  ##Reorder the treatments"Cmax IVM",
  COXdata$Treatment <- factor(COXdata$Treatment, levels = c("Control", "3.8", "37",  "350", "3500" ))
  
  
  #Switch 0s and 1s
  COXdata$Events<-ifelse(COXdata$Events==1,0,1)
  
  return(COXdata)
}

Kisumu_1 = kaplan_Tiassale_Func (Kisumu_one, Kisumu_one)
Tiassale_1 = kaplan_Tiassale_Func (Tiassale_one, Tiassale_one)

Kisumu_2 = kaplan_Tiassale_Func (Kisumu_two, Kisumu_two)
Tiassale_2 = kaplan_Tiassale_Func (Tiassale_two, Tiassale_two)

Kisumu_3 = kaplan_Tiassale_Func (Kisumu_three, Kisumu_three)
Tiassale_3 = kaplan_Tiassale_Func (Tiassale_three, Tiassale_three)




Paste_6 = c("black", "#FF6666", "#CC9933",  "#00CCFF", "#FF66FF")




Survival_curv_Func = function (data, legend.labs, Title) {
  data$surv = with(data, Surv(Times, Events == 1))
  
  curve <- survfit(surv ~ Treatment, data = data, conf.type = "log-log")
  survplot = ggsurvplot(curve, data,
                        conf.int = FALSE,
                        palette =  Paste_6,
                        linetype = c(1,1,1,1,2),
                        size=0.7,
                        title = Title ,
                        # font.title = c(9, "bold"),
                        
                        pval = TRUE,
                        pval.size = 3,
                        pval.coord = c(10, 0.2),
                        ylab="Mosquito survival (%)", xlab="Time (Days)",
                        risk.table = FALSE,
                        risk.table.title = "",
                        
                        xlim=c(0,14),
                        legend.title = "Conc. [ng/mL]",
                        risk.table.height = 0.2,
                        
                        #legend = c(0.75,0.8),
                        ggtheme = PlotTHEME,
                        break.time.by = 4,
                        
                        legend.labs = legend.labs) +
    guides(colour = guide_legend(nrow=1))
  
  
  
  return(survplot)
  
}



survival_curve_kisumu_1 = Survival_curv_Func  (Kisumu_1,
                                             legend.labs = c("Control", "3.8", "37",
                                                             "350", "3500"),
                                             Title = " Kisumu Replicate 1")


survival_curve_tiassale_1= Survival_curv_Func  (Tiassale_1,
                                             legend.labs = c("Control", "3.8", "37",
                                                             "350", "3500"),
                                             Title = " Tiassalé Replicate 1")

survival_curve_kisumu_2 = Survival_curv_Func  (Kisumu_2,
                                               legend.labs = c("Control", "3.8", "37",
                                                               "350", "3500"),
                                               Title = " Kisumu Replicate 2")


survival_curve_tiassale_2= Survival_curv_Func  (Tiassale_2,
                                                legend.labs = c("Control", "3.8", "37",
                                                                "350", "3500"),
                                                Title = " Tiassalé Replicate 2")


survival_curve_kisumu_3 = Survival_curv_Func  (Kisumu_3,
                                               legend.labs = c("Control", "3.8", "37",
                                                               "350", "3500"),
                                               Title = " Kisumu Replicate 3")


survival_curve_tiassale_3= Survival_curv_Func  (Tiassale_3,
                                                legend.labs = c("Control", "3.8", "37",
                                                                "350", "3500"),
                                                Title = " Tiassalé Replicate 3")



survival_curve_kisumu_1$plot <-survival_curve_kisumu_1$plot + labs(tag='A');
survival_curve_tiassale_1$plot <-survival_curve_tiassale_1$plot + labs(tag='B');
survival_curve_kisumu_2$plot <-survival_curve_kisumu_2$plot + labs(tag='C');
survival_curve_tiassale_2$plot <-survival_curve_tiassale_2$plot + labs(tag='D');
survival_curve_kisumu_3$plot <-survival_curve_kisumu_3$plot + labs(tag='E');
survival_curve_tiassale_3$plot <-survival_curve_tiassale_3$plot + labs(tag='F');


plots = list(survival_curve_kisumu_1, survival_curve_kisumu_2,survival_curve_kisumu_3, survival_curve_tiassale_1, survival_curve_tiassale_2, survival_curve_tiassale_3)

Figure_5_TIFF = arrange_ggsurvplots (plots, nrow=3, ncol=2)

tiff("MS_Fig_5.tiff", width = 7.3, height = 11, units = "in", res = 600)
print(Figure_5_TIFF)  # Use print() to ensure the plot is rendered
dev.off() 


############# Figure 6  ##########


single_HRs = readRDS("PKProfiles/single_HRs.RDS")  ##Generated in R script 'VirtualSimulation'
single_Ccs = readRDS("PKProfiles/single_Ccs.RDS")  ##Generated in R script 'VirtualSimulation'


IVM = readRDS ("PKProfiles/IVM_300_profile.RDS")

Cc_Func = function (GROUPEDHR) {
  
  
  data_HR_med = GROUPEDHR %>%
    filter (median>136)
  data_HR_05 = GROUPEDHR %>%
    filter (percen5>136)
  data_HR_95 = GROUPEDHR %>%
    filter (percen95>136)
  data_HR_25 = GROUPEDHR %>%
    filter (percen25>136)
  data_HR_75 = GROUPEDHR %>%
    filter (percen75>136)
  
  THR_med = max(data_HR_med$TIME)/24
  THR_05 = max(data_HR_05$TIME)/24
  THR_95 = max(data_HR_95$TIME)/24
  THR_25 = max(data_HR_25$TIME)/24
  THR_75 = max(data_HR_75$TIME)/24
  
  THR_med = unlist(THR_med)
  THR_05 = unlist(THR_05)
  THR_95 = unlist(THR_95)
  THR_25 = unlist(THR_25)
  THR_75 = unlist(THR_75)
  
  THR_med = as.data.frame(THR_med)
  THR_05 = as.data.frame(THR_05)
  THR_95 = as.data.frame(THR_95)
  THR_25 = as.data.frame(THR_25)
  THR_75 = as.data.frame(THR_75)
  
  time_above = cbind (THR_med, THR_05, THR_95, THR_75, THR_25)
  return(time_above)
  
}
Single_Cc<-lapply(single_Ccs, FUN= Cc_Func)
Cc_sim <- plyr::ldply(Single_Cc, data.frame)
Cc_sim$ID = rownames(Cc_sim)
Cc_sim$ID = as.numeric(Cc_sim$ID)

HR_Func = function (GROUPEDHR) {
  
  data_HR_med = GROUPEDHR %>%
    filter (median>4)
  data_HR_05 = GROUPEDHR %>%
    filter (percen5>4)
  data_HR_95 = GROUPEDHR %>%
    filter (percen95>4)
  data_HR_25 = GROUPEDHR %>%
    filter (percen25>4)
  data_HR_75 = GROUPEDHR %>%
    filter (percen75>4)
  
  THR_med = max(data_HR_med$TIME)/24
  THR_05 = max(data_HR_05$TIME)/24
  THR_95 = max(data_HR_95$TIME)/24
  THR_25 = max(data_HR_25$TIME)/24
  THR_75 = max(data_HR_75$TIME)/24
  
  THR_med = unlist(THR_med)
  THR_05 = unlist(THR_05)
  THR_95 = unlist(THR_95)
  THR_25 = unlist(THR_25)
  THR_75 = unlist(THR_75)
  
  THR_med = as.data.frame(THR_med)
  THR_05 = as.data.frame(THR_05)
  THR_95 = as.data.frame(THR_95)
  THR_25 = as.data.frame(THR_25)
  THR_75 = as.data.frame(THR_75)
  
  time_above = cbind (THR_med, THR_05, THR_95, THR_25, THR_75)
  return(time_above)
  
}
Single_HR<-lapply(single_HRs, FUN= HR_Func)
HR_sim <- plyr::ldply(Single_HR, data.frame)
HR_sim$ID = rownames(HR_sim)
HR_sim$ID = as.numeric(HR_sim$ID)



HR4_plot =
  ggplot(HR_sim) +
  geom_line(aes(ID*0.01, THR_med),col=NTBCcolor, size=0.5) +
  geom_ribbon(aes(ymin=THR_05, ymax=THR_95,x=ID*0.01), fill=NTBCcolor, alpha=0.4) +
  geom_ribbon(aes(ymin=THR_25, ymax=THR_75,x=ID*0.01), fill=NTBCcolor, alpha=0.4) +
  
  geom_hline(yintercept = 11, linetype=2, col="dodgerblue3", size=0.5)+
  annotate("text", x = 0.68, y= 10 ,  colour = "dodgerblue3", label="3x0.6mg/kg IVM: 11 days", size=3)+
  
  scale_y_continuous("Time above HR=4 (Days)", limits = c(0,28), breaks = c(0,5,10,15,20,25)) +
  
  scale_x_continuous("Single dose of NTBC (mg/kg)", limits = c(0,1), breaks = c(0, 1, 0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  
  scale_color_manual(values=c(NTBCcolor)) +
  scale_fill_manual(values=c(NTBCcolor))  +
  labs(title = expression(paste( " Time above HR=4")), tag =  expression(paste(bold("B")))) +
  

PlotTHEME + 
  
  theme( axis.title.y = element_text(family="Myriad Pro", size=9, vjust=2))

HR4_plot

LC50_plot =
  ggplot(Cc_sim) +
  geom_line(aes(ID*0.01, THR_med),col=NTBCcolor, size=0.5) +
  geom_ribbon(aes(ymin=THR_05, ymax=THR_95,x=ID*0.01), fill=NTBCcolor, alpha=0.4) +
  geom_ribbon(aes(ymin=THR_25, ymax=THR_75,x=ID*0.01), fill=NTBCcolor, alpha=0.4) +
  
  geom_hline(yintercept = 8, linetype=2, col="dodgerblue3",size=0.5)+
  annotate("text", x = 0.68, y= 9 ,  colour = "dodgerblue3",  label="3x0.6mg/kg IVM: 8 days", size=3)+
 
  scale_y_continuous(expression(paste("Time above LC"["50"], "(Days)")), limits = c(0,28), breaks = c(0,10,20)) +
  scale_x_continuous("Single dose of NTBC (mg/kg)", limits = c(0,1), breaks = c(0, 1, 0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  
  scale_color_manual(values=c(NTBCcolor)) +
  scale_fill_manual(values=c(NTBCcolor))  +
  labs(title = expression(paste( " Time above LC"["50"])), tag =  expression(paste(bold("A")))) +

  #labs(y=expression(bold(paste("L",og["10"]," cell concentration (cells ",ml^"-1",")",sep=""))))
  PlotTHEME + 
  
theme( axis.title.y = element_text(family="Myriad Pro", size=9, vjust=2))

LC50_plot

Fig_6 = plot_grid(LC50_plot, HR4_plot)



tiff("MS_Fig_6.tiff", width = 7.3, height = 3, units = "in", res = 600)
print(Fig_6)  # Use print() to ensure the plot is rendered
dev.off() 




############# Figure 7 ##########


simulation_2_31 = readRDS("PKProfiles/simulation_2_31.rds")

label1 ="EC[50]"
label2=  " = 205.31ng/mL"




PK_plot =
  ggplot() +
  geom_line(data=simulation_2_31, aes(TIME/24, Cc, col=NTBCcolor),  size=0.5) +

  scale_y_log10("Predicted plasma conc. (ng/mL)", limits = c(1,1000)) +
  scale_x_continuous("Time (Days)", breaks=seq(0,31,7))+
  scale_color_manual(values=c( NTBCcolor, IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_fill_manual(values=c(NTBCcolor,IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_linetype_manual(values=c(1,1),            labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  guides(colour = guide_legend(override.aes = list(fatten = 15))) +
  geom_hline(yintercept = 205, size=0.4, linetype=2, col=NTBCcolor) +
  
  labs(title =  " Predicted PK profile", tag = expression(paste(bold("A")))) +

  geom_text(aes( 17.3, 260, label = label1), size=3, col=NTBCcolor, parse=TRUE) +
  geom_text(aes(23, 260, label = label2), size=3, color = NTBCcolor)  +
  
  
  PlotTHEME +
  theme(legend.position = "none")



simulation_2_31$HR[simulation_2_31$HR >= 20] <- 20

label3 = "HR = 4"
  
PD_plot =
  ggplot() +
  geom_line(data=simulation_2_31, aes(TIME/24, HR, col=NTBCcolor),  size=0.5) +
  
  scale_y_continuous("Predicted HR",breaks=seq(0,20,5)) +
  scale_x_continuous("Time (Days)", breaks=seq(0,31,7))+
  scale_color_manual(values=c( NTBCcolor, IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_fill_manual(values=c(NTBCcolor,IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_linetype_manual(values=c(1,1),            labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  guides(colour = guide_legend(override.aes = list(fatten = 15))) +
  geom_hline(yintercept = 4, size=0.4, linetype=2, col=NTBCcolor) +
  
  labs(title =  " Predicted PD profile", tag = expression(paste(bold("D")))) +
  
  geom_text(aes( 22, 4.7, label = label3), size=3, col=NTBCcolor, parse=FALSE) +

  
  PlotTHEME +
  theme(legend.position = "none")


library(cowplot)

Fig_7 = plot_grid(PK_plot, PD_plot)


tiff("MS_Fig_7.tiff", width = 7.3, height = 3, units = "in", res = 600)
print(Fig_7)  # Use print() to ensure the plot is rendered
dev.off() 



############# Figure 8 ##########

##Feb
Exp_1<-read_excel(path="AKU_experiments/AKU_Exp1.xlsx",
                  sheet = "Sheet1", range = "A2:Q19", col_names = TRUE)

names(Exp_1)=c("Time", "S1_1",  "S1_2",  "S1_5", "S1_10", "S2_1",  "S2_2",  "S2_5", "S2_10",
               "S3_1",  "S3_2",  "S3_5","S3_10" , "S4", "Sub-lethal", "Lethal", "Control")


Survival_curv_Func = function (data, legend.labs, Title, paste, lines) {
  data$surv = with(data, Surv(Times, Events == 1))
  
  curve <- survfit(surv ~ Treatment, data = data, conf.type = "log-log")
  survplot = ggsurvplot(curve, data,
                        conf.int = FALSE,
                        palette =  paste,
                        linetype = lines,
                        size=0.7,
                        title = Title ,
                        # font.title = c(9, "bold"),
                        
                        pval = TRUE,
                        pval.size = 3,
                        pval.coord = c(10, 0.2),
                        ylab="Mosquito survival (%)", xlab="Time (Days)",
                        risk.table = FALSE,
                        risk.table.title = "",
                        
                        xlim=c(0,15),
                        legend.title = "Conc. [ng/mL]",
                        risk.table.height = 0.2,
                        
                        #legend = c(0.75,0.8),
                        ggtheme = PlotTHEME,
                        break.time.by = 4,
                        
                        legend.labs = legend.labs) +
    guides(colour = guide_legend(nrow=1))
  
  
  
  return(survplot)
  
}



#Convert to Kaplan format
Exp_1_suvr = Kaplan_Func (Exp_1)

Exp_suvr = Exp_1_suvr

Exp_suvr$Treatment =as.character(Exp_suvr$Treatment)
Exp_suvr$Patient <- sapply(strsplit(Exp_suvr$Treatment, "_"), "[", 1)

Exp_suvr$Times = Exp_suvr$Times/24


S_1 = Exp_suvr %>%
  filter(Patient=="S1")

S_2 = Exp_suvr %>%
  filter(Patient=="S2")

S_3 = Exp_suvr %>%
  filter(Patient=="S3")

S_4 = Exp_suvr %>%
  filter(Patient=="S4")

S_5 = Exp_suvr %>%
  filter(Patient=="S5")

Control = Exp_suvr %>%
  filter(Patient=="Control")

S_1 = rbind (S_1, Control)
S_2 = rbind (S_2, Control)
S_3 = rbind (S_3, Control)
S_4 = rbind (S_4, Control)
S_5 = rbind (S_5, Control)



#Survival curve - for Experiment 1 and 2


Paste_5_cyan = c("#008B8B",  "#00CDCD", "#00EEEE", "#00FFFF", "#97FFFF")
Paste_5_turq = c("#00868B",  "#00C5CD", "#00F5FF", "#40E0D0", "#97FFFF")

#green
Paste_5_green = c("black",  "#00868B", "#5CACEE", "#8B5A2B", "#FFA54F")
Paste_1 =  c("black", "#663399")

#purple
Paste_5_purple = c("#999999",  "#663399", "#9966CC", "#CC99FF")



#Experiment 1 - Feb
S_1$Treatment <- factor(S_1$Treatment, levels = c("Control","S1_1",  "S1_2", "S1_5", "S1_10"))
S_2$Treatment <- factor(S_2$Treatment, levels = c("Control","S2_1",  "S2_2", "S2_5", "S2_10"))
S_3$Treatment <- factor(S_3$Treatment, levels = c("Control","S3_1",  "S3_2", "S3_5", "S3_10"))
S_4$Treatment <- factor(S_4$Treatment, levels = c("Control","S4"))

#Experiment 1 - Feb
S1_curve = Survival_curv_Func (S_1, legend.labs = c("Control", "1:1", "1:2", "1:5",  "1:10"), Title =    " Patient 1", paste = Paste_5_green, lines=c(1,1,1,1,1))
S2_curve = Survival_curv_Func (S_2, legend.labs = c("Control","1:1", "1:2", "1:5",  "1:10"), Title =     " Patient 2", paste = Paste_5_green, lines=c(1,1,1,1,1))
S3_curve = Survival_curv_Func (S_3, legend.labs = c("Control","1:1", "1:2", "1:5",  "1:10"), Title =    " Patient 3",   paste = Paste_5_green, lines=c(1,1,1,1,1))
S4_curve = Survival_curv_Func (S_4, legend.labs = c("Control","1:1"), Title =   " Patient 4", paste = Paste_1, lines=c(1,1,1))


S1_curve$plot <-S1_curve$plot + labs(tag='A');
S2_curve$plot <-S2_curve$plot + labs(tag='B');
S3_curve$plot <-S3_curve$plot + labs(tag='C');
S4_curve$plot <-S4_curve$plot + labs(tag='D');

plots = list(S1_curve, S3_curve, S2_curve, S4_curve)

Figure_8_TIFF = arrange_ggsurvplots (plots, nrow=2, ncol=2)

tiff("MS_Fig_8.tiff", width = 7.3, height = 7.3, units = "in", res = 600)
print(Figure_8_TIFF)  # Use print() to ensure the plot is rendered
dev.off() 



############# Figure 9 ##########

##Need to open r script: ntbc_tyr_coxreg



Tyr_Sept<-read_excel(path="Tyrosine_experiments/Tyrosine_exp1_rep1.xlsx",
                     sheet = "Sheet1", range = "A2:K18", col_names = TRUE)

names(Tyr_Sept)=c("Time", "low_tyr", "high_tyr", "lowtyr_3.8", "hightyr_3.8", "lowtyr_37", "hightyr_37",
                  "NTBC250", "NTBC37", "NTBC3.8", "Control")

Tyr_Oct<-read_excel(path="Tyrosine_experiments/Tyrosine_exp1_rep2.xlsx",
                    sheet = "Sheet1", range = "A2:K18", col_names = TRUE)

names(Tyr_Oct)=c("Time", "low_tyr", "high_tyr", "lowtyr_3.8", "hightyr_3.8", "lowtyr_37", "hightyr_37",
                 "NTBC250", "NTBC37", "NTBC3.8", "Control")


Tyr_Nov<-read_excel(path="Tyrosine_experiments/Tyrosine_exp1_rep3.xlsx",
                    sheet = "Sheet1", range = "A2:L18", col_names = TRUE)

names(Tyr_Nov)=c("Time", "low_tyr", "high_tyr", "lowtyr_3.8", "hightyr_3.8", "lowtyr_37", "hightyr_37",
                 "NTBC250", "NTBC37", "NTBC3.8", "Control", "hightyr_150")





Kaplan_Func = function (kaplan) {
  
  kaplan=kaplan[c(1,3:nrow(kaplan)),]
  
  names(kaplan)[1]="Time"
  kaplan$Time = as.numeric(kaplan$Time)
  kaplan[1,1]=0
  
  
  kaplanS=kaplan
  
  for(i in 2:nrow(kaplan)){
    
    for (j in 2:length(kaplan)){
      kaplanS[i,j]=kaplan[1,j]-kaplan[i,j]
    }
  }
  
  
  #Melt the data into a column format (easier to work with for R)
  MeltedKaplan<-reshape2::melt(data = kaplanS, id.vars = 1)
  
  ###split the data into elements representing a treatment each
  KaplanList<-split(MeltedKaplan, f = as.factor(MeltedKaplan$variable))
  
  x<-KaplanList[[2]]
  ##Function to convert data to COX/Kaplan Meyer Format
  COXIT <-function (x) { ##you can text the function by defining x as x<-KaplanList[[1]] - so then you can run it step by step and see what it's doing
    ##Step 1 - Change variable column from factor to character
    x$variable<-as.character(x$variable)
    
    ##Step -2 Create a column to count events for each time point (The reduction of mosquitoes from the last measurement at each time point)
    x$Events<-NA
    for (i in (1:nrow(x))) {
      x[i,4]<-ifelse(i==1, 0, x[i-1,3]-x[i,3])
    }
    
    ##Add a row for remaining non-events ### This is basically the number of mosquitoes that remain alive at the end of the experiment. The last value
    x[nrow(x)+1,]<-c(x[nrow(x),1] ,unique(x$variable),x[nrow(x),3], x[nrow(x),3])
    
    
    ##Add a column which defines whether Event is life or death (0 or 1)
    x$Nature<-c(rep(0, nrow(x)-1),1)
    
    
    ##Change Time, value and Events columns to numeric
    x$Time<-as.numeric(x$Time)
    x$value<-as.numeric(x$value)
    x$Events<-as.numeric(x$Events)
    
    ##Remove time points not containing events. They're not needed for a COX regression
    x<-subset(x, !x$Events==0)
    
    ##Repeat Times and Events by a factor of number events. If you have 3 mosquitoes that die at time 12 you want three entries for time 12, if you have 10 more mosquitoes that die at time 24 you want 10 entries for time 24 and so on
    Times<-rep(x$Time, x$Events)
    
    ##Same as for times. Events are coded with 0. We want to record total number of events (deaths) at the last time point and then code all reamining living mosquiotes, quantified by the last added row in the events column with   1
    
    if (1 %in% x$Nature) { ##Check if there are any living mosquitoes left at the end of the experiment , if there are then rep 0 for all events and 1 for remaining living mosquitoes
      Events<-rep(c(rep(0,(nrow(x)-1)),1), x$Events)
    } else {  #if no living mosquitoes remaining use this: Only repeat 0 (for events) for all, no 1s because no living mosquitoes
      Events<-rep(rep(0,nrow(x)), x$Events)
    }
    ##0= Event of death, 1=censored at end of experiment/no event
    variable<-rep(unique(x$variable), sum(x$Events)) ##repeat the variable to insert with the data
    
    ##Create the data
    COX<-data.frame("Times"=Times, "Events"=Events, "Treatment"= variable)
    return(COX)
  }
  
  
  #Apply the function on the new list
  newlist<-lapply(KaplanList, FUN= COXIT)
  
  ##
  COXdata<-bind_rows(newlist)
  names(COXdata)
  COXdata$Treatment<-as.factor(COXdata$Treatment)
  
  
  #Switch 0s and 1s
  COXdata$Events<-ifelse(COXdata$Events==1,0,1)
  
  return (COXdata)
  
}

Kaplan_Sept = Kaplan_Func(Tyr_Sept)
Kaplan_Oct = Kaplan_Func(Tyr_Oct)
Kaplan_Nov = Kaplan_Func(Tyr_Nov)

Kaplan_all = rbind(Kaplan_Sept, Kaplan_Oct, Kaplan_Nov)

  Kaplan_all$Times=Kaplan_all$Times/24
  final = Kaplan_all
  
  NTBC <- final[grep("NTBC", final$Treatment), ]
  Control <- final[grep("Control", final$Treatment), ]
  NTBC3.8 <- final[grep("3.8", final$Treatment), ]
  NTBC37 <- final[grep("37", final$Treatment), ]
  low_tyr <- final[grep("low_tyr", final$Treatment), ]
  high_tyr <- final[grep("high_tyr", final$Treatment), ]
  
  tyr = rbind (Control, low_tyr, high_tyr)
  NTBC = rbind (Control, NTBC)
  
  
  
  ###Survival analysis, survival curves
  Survival_curv_Func = function (data, legend.labs, Title, paste, lines) {
    data$surv = with(data, Surv(Times, Events == 1))
    
    curve <- survfit(surv ~ Treatment, data = data, conf.type = "log-log")
    survplot = ggsurvplot(curve, data,
                          conf.int = FALSE,
                          palette =  paste,
                          linetype = lines,
                          size=0.5,
                          title = Title ,
                          # font.title = c(9, "bold"),
                          
                          pval = TRUE,
                          pval.size = 3,
                          pval.coord = c(10, 0.2),
                          ylab="Mosquito survival (%)", xlab="Time (Days)",
                          risk.table = FALSE,
                          risk.table.title = "",
                          
                          xlim=c(0,14),
                          legend.title = "Conc. [ng/mL]",
                          risk.table.height = 0.2,
                          
                          #legend = c(0.75,0.8),
                          ggtheme = PlotTHEME,
                          break.time.by = 4,
                          
                          legend.labs = legend.labs) +
      guides(colour = guide_legend(nrow=1))
    
    
    
    return(survplot)
    
  }
  

  NTBC$Treatment <- factor(NTBC$Treatment, levels = c("Control",  "NTBC3.8", "NTBC37", "NTBC250"))
  tyr$Treatment <- factor(tyr$Treatment, levels = c("Control",   "low_tyr", "high_tyr"))
  NTBC3.8$Treatment <- factor(NTBC3.8$Treatment, levels = c("NTBC3.8",  "lowtyr_3.8", "hightyr_3.8"))
  NTBC37$Treatment <- factor(NTBC37$Treatment, levels =  c("NTBC37",  "lowtyr_37", "hightyr_37"))
  
  
  ##Code for tyrosine experiment data results
  NTBC_concs_surv = Survival_curv_Func (NTBC, legend.labs = c("Control",  "3.8 NTBC", "37 NTBC", "250 NTBC"), Title =  " NTBC concentrations  (ng/mL)" , paste = Paste_4, lines=c(1,1,1,1))
  Tyr_concs_surv = Survival_curv_Func (tyr, legend.labs = c("Control", "Low tyr", "High tyr"), Title =    " Tyrosine concentrations" , paste = Paste_4, lines=c(1,1,1))
  NTBC3.8_concs_surv = Survival_curv_Func (NTBC3.8, legend.labs = c("3.8 NTBC",  "NTBC+low tyr", "NTBC+high tyr"), Title =   " Non-lethal NTBC 3.7ng/mL" , paste = Paste_4, lines=c(1,1,1))
  NTBC37_concs_surv = Survival_curv_Func (NTBC37, legend.labs = c("37 NTBC",  "NTBC+low tyr", "NTBC+high tyr"), Title =  " Sub-lethal NTBC 37ng/mL" , paste = Paste_4, lines=c(1,1,1))
  
  
  
  NTBC_concs_surv$plot <-NTBC_concs_surv$plot + labs(tag='A');
  NTBC3.8_concs_surv$plot <-NTBC3.8_concs_surv$plot + labs(tag='C');
  Tyr_concs_surv$plot <-Tyr_concs_surv$plot + labs(tag='B');
  NTBC37_concs_surv$plot <-NTBC37_concs_surv$plot + labs(tag='D');
  
  plots = list(NTBC_concs_surv, NTBC3.8_concs_surv, Tyr_concs_surv, NTBC37_concs_surv)

Figure_9_TIFF = arrange_ggsurvplots (plots, nrow=2, ncol=2)

tiff("MS_Fig_9.tiff", width = 7.3, height = 7.3, units = "in", res = 600)
print(Figure_9_TIFF)  # Use print() to ensure the plot is rendered





############# Figure S1 ##########


#Set working directory


#Experiment 1 data, (Experiment 1 of young mosquitoes experiment)
#This is already in Kaplan ready format
COX_1 = read_csv("00-OriginalData/NTBC_IVM_Exp1.csv")
COX_1$Treatment = str_replace_all(COX_1$Treatment, "Ethanol.Control", "Control 0"	)

COX_1$Treatment=revalue(COX_1$Treatment, c("NTBC_10.000.ng.ul" ="NTBC 10000","IVM_1000.ng.ul"="IVM 1000","IVM_5000.ng.ul"="IVM 5000",
                                           "IVM_100.ng.ul"="IVM 100",  "NTBC_250.ng.ul" = "NTBC 250" ,  "IVM_15.ng.ul"="IVM 15",
                                           "NTBC_100.ng.ul"="NTBC 100", "NTBC_50.ng.ul"="NTBC 50"))


COX_1$Times = COX_1$Times/24


#Split data into drugs
NTBC_COX_Func = function (final) {
  final$Treatment = as.character(final$Treatment)
  
  final$Drug=gsub("Treatment", "", final$Treatment)
  final$Drug <- sapply(strsplit(final$Drug, " "), "[", 1)
  final$Conc <- sapply(strsplit(final$Treatment, " "), "[", 2)
  
  Control<-subset(final, Drug=="Control") ###NTBC data
  IVM<-subset(final, Drug=="IVM") ###NTBC data
  NTBC<-subset(final, Drug=="NTBC") ###NTBC data
  
  NTBC = rbind (NTBC, Control)
  IVM = rbind (IVM, Control)
  
  
  NTBC = NTBC [,-3]
  return (NTBC)
}
IVM_COX_Func = function (final) {
  final$Treatment = as.character(final$Treatment)
  
  final$Drug=gsub("Treatment", "", final$Treatment)
  final$Drug <- sapply(strsplit(final$Drug, " "), "[", 1)
  final$Conc <- sapply(strsplit(final$Treatment, " "), "[", 2)
  
  Control<-subset(final, Drug=="Control") ###NTBC data
  IVM<-subset(final, Drug=="IVM") ###NTBC data
  NTBC<-subset(final, Drug=="NTBC") ###NTBC data
  
  NTBC = rbind (NTBC, Control)
  IVM = rbind (IVM, Control)
  IVM = IVM [,-3]
  
  return (IVM)
}

NTBC_1 = NTBC_COX_Func (COX_1)
IVM_1 = IVM_COX_Func (COX_1)

IVM_1$Conc <- factor(IVM_1$Conc, levels = c("0", "15", "100", "1000", "5000"))
NTBC_1$Conc <- factor(NTBC_1$Conc, levels = c("0", "50", "100", "250", "10000"))
colnames(NTBC_1) [4] = "Treatment"
colnames(IVM_1) [4] = "Treatment"

NTBC_curve_1 = Survival_curv_Func (data = NTBC_1, legend.labs = c("Control", "50", "100", "250", "10,000"), Title =   "NTBC - Experiment 1", paste = Paste_5, lines=c(1,1,1,1,2))

IVM_curve_1 = Survival_curv_Func (data = IVM_1, legend.labs = c("Control", "15", "100", "1000", "5,000"), Title =   "IVM - Experiment 1", paste = Paste_5, lines=c(1,1,1,1,2))


NTBC_curve_1$plot <-NTBC_curve_1$plot + labs(tag='A');
IVM_curve_1$plot <-IVM_curve_1$plot + labs(tag='B');

plots = list(NTBC_curve_1, IVM_curve_1)

Figure_S1_TIFF = arrange_ggsurvplots (plots, nrow=1, ncol=2)

tiff("MS_Fig_S1.tiff", width = 7.3, height = 3.5, units = "in", res = 600)
print(Figure_S1_TIFF)  # Use print() to ensure the plot is rendered
dev.off() 






############# Figure S2 ##########


#Experiment 3 data, (Experiment 3 of young mosquitoes experiment, only including extra concentrations of NTBC)


conc_exp = read_excel(path= "00-OriginalData/NTBC_Exp3.xlsx",
                      sheet = "Sheet1", range = "A2:G30", col_names = TRUE)


names(conc_exp)=c("Time",  "Control 0",   "NTBC 250", "NTBC 200", "NTBC 150", "NTBC 100",
                  "IVM 15")


Kaplan_Func_3 = function (kaplan) {
  
  kaplan = conc_exp
  kaplan = kaplan [-2 , ]
  kaplan=kaplan[c(1,3:nrow(kaplan)),]
  
  kaplan[2,]=kaplan[1,]
  kaplan=kaplan[-1,]
  
  
  names(kaplan)[1]="Time"
  kaplan$Time = as.numeric(kaplan$Time)
  kaplan[1,1]=0
  
  kaplan$`...11`=NULL
  
  
  kaplanS=kaplan
  
  for(i in 2:nrow(kaplan)){
    
    for (j in 2:length(kaplan)){
      kaplanS[i,j]=kaplan[1,j]-kaplan[i,j]
    }
  }
  kaplanS=kaplanS[1:25,]
  

  #Melt the data into a column format (easier to work with for R)
  MeltedKaplan<-reshape2::melt(data = kaplanS, id.vars = 1)
  
  ###split the data into elements representing a treatment each
  KaplanList<-split(MeltedKaplan, f = as.factor(MeltedKaplan$variable))
  
  x<-KaplanList[[2]]
  ##Function to convert data to COX/Kaplan Meyer Format
  COXIT <-function (x) { ##you can text the function by defining x as x<-KaplanList[[1]] - so then you can run it step by step and see what it's doing
    ##Step 1 - Change variable column from factor to character
    x$variable<-as.character(x$variable)
    
    ##Step -2 Create a column to count events for each time point (The reduction of mosquitoes from the last measurement at each time point)
    x$Events<-NA
    for (i in (1:nrow(x))) {
      x[i,4]<-ifelse(i==1, 0, x[i-1,3]-x[i,3])
    }
    
    ##Add a row for remaining non-events ### This is basically the number of mosquitoes that remain alive at the end of the experiment. The last value
    x[nrow(x)+1,]<-c(x[nrow(x),1] ,unique(x$variable),x[nrow(x),3], x[nrow(x),3])
    
    
    ##Add a column which defines whether Event is life or death (0 or 1)
    x$Nature<-c(rep(0, nrow(x)-1),1)
    
    
    ##Change Time, value and Events columns to numeric
    x$Time<-as.numeric(x$Time)
    x$value<-as.numeric(x$value)
    x$Events<-as.numeric(x$Events)
    
    ##Remove time points not containing events. They're not needed for a COX regression
    x<-subset(x, !x$Events==0)
    
    ##Repeat Times and Events by a factor of number events. If you have 3 mosquitoes that die at time 12 you want three entries for time 12, if you have 10 more mosquitoes that die at time 24 you want 10 entries for time 24 and so on
    Times<-rep(x$Time, x$Events)
    
    ##Same as for times. Events are coded with 0. We want to record total number of events (deaths) at the last time point and then code all reamining living mosquiotes, quantified by the last added row in the events column with   1
    
    if (1 %in% x$Nature) { ##Check if there are any living mosquitoes left at the end of the experiment , if there are then rep 0 for all events and 1 for remaining living mosquitoes
      Events<-rep(c(rep(0,(nrow(x)-1)),1), x$Events)
    } else {  #if no living mosquitoes remaining use this: Only repeat 0 (for events) for all, no 1s because no living mosquitoes
      Events<-rep(rep(0,nrow(x)), x$Events)
    }
    ##0= Event of death, 1=censored at end of experiment/no event
    variable<-rep(unique(x$variable), sum(x$Events)) ##repeat the variable to insert with the data
    
    ##Create the data
    COX<-data.frame("Times"=Times, "Events"=Events, "Treatment"= variable)
    return(COX)
  }
  
  
  #Apply the function on the new list
  newlist<-lapply(KaplanList, FUN= COXIT)
  
  ##
  COXdata<-bind_rows(newlist)
  names(COXdata)
  COXdata$Treatment<-as.factor(COXdata$Treatment)
  
  
  #Switch 0s and 1s
  COXdata$Events<-ifelse(COXdata$Events==1,0,1)
  
  
  return (COXdata)
  
}

conc_COX = COXdata

conc_COX = Kaplan_Func_3 (conc_exp)
COX_3 = conc_COX

COX_3$Times = COX_3$Times/24
NTBC_3 = NTBC_COX_Func (COX_3)

NTBC_3$Conc <- factor(NTBC_3$Conc, levels = c("0", "100", "150", "200", "250"))


colnames(NTBC_3) [4] = "Treatment"


Survival_curv_Func = function (data, legend.labs, Title, paste, lines) {
  data$surv = with(data, Surv(Times, Events == 1))
  
  curve <- survfit(surv ~ Treatment, data = data, conf.type = "log-log")
  survplot = ggsurvplot(curve, data,
                        conf.int = FALSE,
                        palette =  paste,
                        linetype = lines,
                        size=0.7,
                        title = Title ,
                        # font.title = c(9, "bold"),
                        
                        pval = TRUE,
                        pval.size = 3,
                        pval.coord = c(10, 0.4),
                        ylab="Mosquito survival (%)", xlab="Time (Days)",
                        risk.table = FALSE,
                        risk.table.title = "",
                        
                        xlim=c(0,15),
                        legend.title = "Conc. [ng/mL]",
                        risk.table.height = 0.2,
                        
                        #legend = c(0.75,0.8),
                        ggtheme = PlotTHEME,
                        break.time.by = 4,
                        
                        legend.labs = legend.labs) +
    guides(colour = guide_legend(nrow=1))
  
  
  
  return(survplot)
  
}



NTBC_curve_3 = Survival_curv_Func (data = NTBC_3, legend.labs = c("Control", "100", "150", "200", "250"), Title =   "NTBC - Experiment 3", paste = Paste_5, lines=c(1,1,1,1,2))


tiff("MS_Fig_S2.tiff", width = 3.5, height = 3.5, units = "in", res = 600)
print(NTBC_curve_3)  # Use print() to ensure the plot is rendered
dev.off() 





############# Figure S3 ##########





HT_NTBC_1_01 = readRDS("PKProfiles/HT_NTBC_1_01.RDS")
HT_NTBC_1_01 = as.data.frame(HT_NTBC_1_01)
label1 ="EC[50]"

TIME = 0 
HR = 0 
Cc = 0

one_row = data.frame (TIME, HR, Cc)
HT_NTBC_1_01 = rbind (one_row ,HT_NTBC_1_01)


PK_plot =
  ggplot() +
  geom_line(data=HT_NTBC_1_01, aes(TIME/24, Cc, col=NTBCcolor),  size=0.5) +
  
  scale_y_log10("Predicted plasma conc. (ng/mL)", limits = c(1,1000)) +
  scale_x_continuous("Time (Days)", breaks=seq(0,31,7))+
  scale_color_manual(values=c( NTBCcolor, IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_fill_manual(values=c(NTBCcolor,IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_linetype_manual(values=c(1,1),            labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  guides(colour = guide_legend(override.aes = list(fatten = 15))) +
  geom_hline(yintercept = 205, size=0.4, linetype=2, col=NTBCcolor) +
  
  labs(title =  " Predicted PK profile", tag = expression(paste(bold("A")))) +
  
  geom_text(aes( 23, 260, label = label1), size=3, col=NTBCcolor, parse=TRUE) +

  
  PlotTHEME +
  theme(legend.position = "none")



HT_NTBC_1_01$HR[HT_NTBC_1_01$HR >= 20] <- 20

label3 = "HR = 4"

PD_plot =
  ggplot() +
  geom_line(data=HT_NTBC_1_01, aes(TIME/24, HR, col=NTBCcolor),  size=0.5) +
  
  scale_y_continuous("Predicted HR",breaks=seq(0,20,5)) +
  scale_x_continuous("Time (Days)", breaks=seq(0,31,7))+
  scale_color_manual(values=c( NTBCcolor, IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_fill_manual(values=c(NTBCcolor,IVMcolor), labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  scale_linetype_manual(values=c(1,1),            labels=c("NTBC (1.0 mg/kg qd 3 d)", "IVM (0.6 mg/kg qd 3 d)"))+
  guides(colour = guide_legend(override.aes = list(fatten = 15))) +
  geom_hline(yintercept = 4, size=0.4, linetype=2, col=NTBCcolor) +
  
  labs(title =  " Predicted PD profile", tag = expression(paste(bold("B")))) +
  
  geom_text(aes( 22, 4.7, label = label3), size=3, col=NTBCcolor, parse=FALSE) +
  
  
  PlotTHEME +
  theme(legend.position = "none")

Fig_S3 = plot_grid(PK_plot, PD_plot)



tiff("MS_Fig_S3.tiff", width = 7.3, height = 3, units = "in", res = 600)
print(Fig_S3)  # Use print() to ensure the plot is rendered
dev.off() 






############# Figure S4

getwd()
##Read the file

new_data<-read_excel(path="Metabolites/Metabolite_intensities_14May2021.xlsx",
                         sheet = "Sheet1", range = "A2:O22", col_names = T)

new_data = as.data.frame(new_data)

keeps = c(names(new_data)[which(grepl("A1|A2|A3|A4|B1|B2|B3|B4|C1|C2|C3|C4", names(new_data)))], "AminoAcid")

dataNEW = new_data[keeps]

groups = ifelse(grepl("A", names(dataNEW)), "A", ifelse(grepl("B", names(dataNEW)), "B", ifelse(grepl("C", names(dataNEW)), "C", "AminoAcid")))

MeltedData = melt(dataNEW, "AminoAcid")

MeltedData$Unique = paste0(MeltedData$AminoAcid, ifelse(grepl("A", MeltedData$variable), "A", ifelse(grepl("B", MeltedData$variable), "B", ifelse(grepl("C", MeltedData$variable), "C", NA))))


str(MeltedData)
MeltedDataNEW = ddply(.data = MeltedData, .variables = c("Unique"), summarise, "MeanValue" = mean(as.numeric(value), na.rm=T), "SD" = sd(value, na.rm=T))

MeltedDataNEW$AminoAcid = substr(MeltedDataNEW$Unique, 1, nchar(MeltedDataNEW$Unique)-1)

MeltedDataNEW$Group = substr(MeltedDataNEW$Unique, nchar(MeltedDataNEW$Unique) , nchar(MeltedDataNEW$Unique))
MeltedDataNEW$Names = factor(ifelse(MeltedDataNEW$Group == "A", "Complete", ifelse(MeltedDataNEW$Group == "B", "Partial", "Control")), levels = c("Control", "Partial", "Complete"))


MeltedDataNEW = as.data.frame(MeltedDataNEW)
MeltedDataNEW$AminoAcid = str_replace_all(MeltedDataNEW$AminoAcid, "\r\n3-Hydroxyphenylacetate", "3-Hydroxyphenylacetate"	)


Tyrosine = filter(MeltedDataNEW, AminoAcid %in%  c("3-Hydroxyphenylacetate", "L-tyrosine",  "L-Phenylalanine", "3-Hydroxyphenyllactate", "3-Hydroxyphenylpyruvate"))

Citrate = filter(MeltedDataNEW, AminoAcid %in%  c("Citrate", "2-Oxoglutarate",  "Cis-Aconitate", "Malate", "Succinate"))

Lipid = filter(MeltedDataNEW, AminoAcid %in%  c("Dodecenedioic acid", "Decenedioic acid",  "(E)-10-Oxo-8-decenoic acid",  "Hydroxybutyrylcarnitine"  ,  "3-Hydroxydodecanedioic acid", "2-Hydroxydecanedioic acid", "Butenylcarnitine" , "O-Butanoylcarnitine" ,  "Propenoylcarnitine" , "Suberic acid"))

library(scales)


Tyrosine_fig = 
  
  ggplot(Tyrosine, aes(x=Names, y=MeanValue, fill=Names)) + 
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  geom_errorbar(mapping=aes(ymin=MeanValue-SD, ymax=MeanValue+SD), size=0.25, width=0.5, position=position_dodge(.9)) +
  scale_y_continuous("Relative Abundance", labels = scales::scientific)+ 
  scale_fill_manual(values=alpha(c("lightblue1", "lightblue2", "lightblue3"), c(0.8, 0.8, 1))) +
  guides(fill=guide_legend(title = "Knockdown")) +
  
  facet_wrap(~AminoAcid, scales="free", ncol=3,  labeller = labeller(cluster = facet_labels)) +
  PlotTHEME +
  theme(legend.background= element_rect(fill="white"),  
        axis.title.x = element_blank(),
        legend.box.background= element_rect(fill="white"),
        panel.background = element_rect(fill="white"), 
        strip.text.x = element_text(family = "Myriad Pro", size = 9, colour="black"),
        
        panel.grid.major.y =  element_line(color="white", size=0.2),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.placement = "outside") 


tiff("MS_Fig_S4.tiff", width = 7.3, height = 4, units = "in", res = 600)
print(Tyrosine_fig)  # Use print() to ensure the plot is rendered
dev.off() 


Fig_S5 = 
  ggplot(Citrate, aes(x=Names, y=MeanValue, fill=Names)) + 
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  geom_errorbar(mapping=aes(ymin=MeanValue-SD, ymax=MeanValue+SD), size=0.25, width=0.5, position=position_dodge(.9)) +
  scale_y_continuous("Relative Abundance", labels = scales::scientific)+ 
  scale_fill_manual(values=alpha(c("#FFE5B4", "#FDC89D", "#F8AE8F"), c(0.8, 0.8, 1))) +
  guides(fill=guide_legend(title = "Knockdown")) +
  
  facet_wrap(~AminoAcid, ncol=3, scales="free", labeller = labeller(cluster = facet_labels)) +
  PlotTHEME +
  theme(legend.background= element_rect(fill="white"),  
        axis.title.x = element_blank(),
        legend.box.background= element_rect(fill="white"),
        panel.background = element_rect(fill="white"), 
        strip.text.x = element_text(family = "Myriad Pro", size = 9, colour="black"),
        
        panel.grid.major.y =  element_line(color="white", size=0.2),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.placement = "outside") 

Fig_S5


tiff("MS_Fig_S5.tiff", width = 7.3, height = 4, units = "in", res = 600)
print(Fig_S5)  # Use print() to ensure the plot is rendered
dev.off() 

library(scales)


Fig_S6 = 
  ggplot(Lipid, aes(x=Names, y=MeanValue, fill=Names)) + 
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  geom_errorbar(mapping=aes(ymin=MeanValue-SD, ymax=MeanValue+SD), size=0.25, width=0.5, position=position_dodge(.9)) +
  scale_y_continuous("Relative Abundance", labels = scales::scientific)+ 
    scale_fill_manual(values=alpha(c("#E9EAEC", "#C0C0C0", "#8D918D"), c(0.8, 0.8, 1))) +
    guides(fill=guide_legend(title = "Knockdown")) +
    
    facet_wrap(~AminoAcid, ncol=3, scales="free", labeller = labeller(cluster = facet_labels)) +
    PlotTHEME +
    theme(legend.background= element_rect(fill="white"),  
          axis.text.x = element_text(family = "Myriad Pro", size = 7, colour="black"),
          axis.title.x = element_blank(),
          
          legend.box.background= element_rect(fill="white"),
          panel.background = element_rect(fill="white"), 
          strip.text.x = element_text(family = "Myriad Pro", size = 8, colour="black"),
          
          panel.grid.major.y =  element_line(color="white", size=0.2),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_blank(), strip.placement = "outside") 

          
  
  
Fig_S6


tiff("MS_Fig_S6.tiff", width = 7.3, height = 10, units = "in", res = 600)
print(Fig_S6)  # Use print() to ensure the plot is rendered
dev.off() 

getwd()
