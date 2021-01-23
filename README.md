# Roseate-Tern-Age-dependent-Survival-Analysis - Honors project
# Writer:  Adam McGregor
# Date Created:  15/12/2020  

# Using R2Ucare &amp; RMark, I perform an age-dependent analysis of survival 
# Best viewed "raw"

# ########################################################################### #
# #-----------------------------Where to start------------------------------# #

# 1) Set working Directory to file with raw data set

     setwd("/Users/adamm/Desktop")
     setwd("Honors Project/R")
     getwd()

# 2) Retrieve required packages for data handling & analysis

     library(RMark)
     library(reshape)
     library(reshape2)
     library(R2ucare)
     library(unmarked)
     library(tidyr)

# 3.1) Before the raw data file was retrieved, it first had to be manipulated.
#      The raw data resided in an excel file. The required data columns, age, 
#      year, month, day & id were selected and dealt with as described below:

#        Id, was given the header - animal_id.

#        Age, was the age determined at capture & thusly given the 
#        header - capture_age.

#        Year, month & day were organised into a date format required for 
#        our used programmes, for example: 29-Aug-15. 
#        It was given the header - date.

#        Only encounters between 2009-2016 were used.

#      Note: Entries for days were sometimes not formatted consistantly, 
#      for example: 17-20. 
#      So, these would be averaged and rounded up, 
#      for example: 17+20=37,  37/2=18.5,   rounded up =19.
#      The specific day was not of importance, as we organised individuals
#      by the year and month of capture.

# 3.2) These columns were copy and pasted & then saved as a CVS file, 
#      using "Notepad++"

# ########################################################################### #
# #----------Creating CMR Encounter Histories For Analysis in R-------------# #

# 1) First, the raw data, saved as a CVS file is read into R.
#    The date was turned into an object called date.code, that helped if R
#    needed to organise via date.

     raw.data<-read.csv("Raw_Data_Tern.cvs",header = TRUE, sep = "")
     raw.data$date.code<-strptime(raw.data$date, format = "%d-%b-%y")
     
# 2) A subset of the raw data frame was created to group encounter histories
#    by the  age at first capture.

     ini.age.df <- as.data.frame(with(raw.data,table(animal_id,capture_age)))
     ini.age.df <- subset(ini.age.df,Freq>0,select=c(capture_age))

# 3) Using the package,Reshape2, I melted down the raw data frame for 
#    manipulation.    

     junk <- melt.data.frame(raw.data,id.var=c("animal_id",
                                               "date.code",
                                               "capture_age",
                                               "year"),
                                               measure.var="detect")

#  4) The melted data frame was cast to create encounter histories, ordered 
#     to include their capture age (age at first capture), and their Id.

      y = dcast(junk, capture_age + animal_id ~ year,fun.aggregate=mean)

#  5) Clean up of encounters containing NA -> 0, and counts >2 to 1.

      y [is.na(y)]=0
      k <- dim(y)[2]
      y[,2:k]<-(y[,2:k]>0)*1

# 6) This creates a data frame of encounters histories, 
#    where only individuals that the age at first capture code = 1.
#    (Safring ring guide, 1 = 0-0.5 yrs old)

     true_age.df <- y[y$capture_age == 1, ]

# 7) This pasty function allows the above data frame to turn into a CH for
#    in the next step.

     pasty<-function(x) 
     {
  
     k<-ncol(x)
     n<-nrow(x)
     out<-array(dim=n)
     for (i in 1:n)
     {
     y<-(x[i,]>0)*1
     out[i]<-paste(y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],sep="")
     }
     return(out)
     }
     
# 8) I now turn the data frame into a CH format readable by RMark.
#    This code selects the columns of encounters of the years 2009-2016, 
#    ignoring the first two "capture age" and "animal id" columns.

     cap_age.df = data.frame(ch=pasty(y[,3:10]),
                             capture_age=ini.age.df$capture_age)

# 9) Now, I created a data frame including CH, where it only has 
#    individuals that entered the data set at birth. This means 
#    analysis of age in RMark, can be treated as true age, rather 
#    than time since seen (TSM).

     true_age.CH <- cap_age.df[cap_age.df$capture_age == 1, ]


# ########################################################################### #
# #---------------------Readying Age Data For GOF Tests---------------------# #

# 1) This removes the columns from the true_age.df 
#    not needed to produce a matrix usable in R2Ucare

     drop <- c("capture_age","animal_id")
     enc.his.age.df = true_age.df[,!(names(true_age.df) %in% drop)]

# 2) Next, I turn the data frame (enc.his.age.df) 
#    into a matrix (enc.his.age.mtx)

     enc.his.age.mtx = data.matrix(enc.his.age.df, rownames.force = NA)

# 3) I then remove the column/row names, to simplify
     colnames(enc.his.age.mtx) <- NULL
     rownames(enc.his.age.mtx) <- NULL

# 4) This is a vector of sample_size, for CH of individuals first captured 
#    at birth.
#    It uses the observed number of rows with capture histories = 330.

     sample_size <- matrix(c(1),nrow = 330,ncol = 1,)

# 5) Make a list of the required matrices used for GOF Tests in R2Ucare,
#    sample size & encounter histories.

     tern_age <- list(enc.his.age.mtx, sample_size)


# ########################################################################### #
# #---------------------------------R2UCare---------------------------------# #

##   CJS GOF Tests:

# 1) Get encounter histories & frequencies required to perform GOF tests

     age.hist = tern_age[[1]]
     age.freq = tern_age[[2]]
  
# 2) CJS GOF overall test results for inspection
     overall_CJS(age.hist, age.freq)
  
# 3) Further inspection of test components
     test3sr_age = test3sr(age.hist, age.freq)
     test3sm_age = test3sm(age.hist, age.freq)
     test2ct_age = test2ct(age.hist, age.freq)
     test2cl_age = test2cl(age.hist, age.freq)  
     

# ########################################################################### #
# #----------------------------------RMark----------------------------------# #
     
     
# 1) Creates needed constructs from the release-recapture history. 

     true_age.proc <- process.data(true_age.CH, 
                              begin.time  = 2009  , 
                              model  = "CJS" ,  
                              initial.ages  = (0)    )
     
# 2) Making the design of parameters used in the model.

     ddl = make.design.data(true_age.proc)
  
  
# 3) These functions can be called for their specific data design in age bins.
          
#    Age structure - a -  age binned as 0-1, 1+ 
     ddl = add.design.data(true_age.proc,ddl, parameter = "Phi", type  = "age",
                                  bins=c(0,1,10),  name = "a",  right = FALSE)
     
#    Age structure - b -  age binned as 0-3, 3+ 
     ddl = add.design.data(true_age.proc,ddl,parameter="Phi",type="age",
                                  bins=c(0,3,10),name="b", right=FALSE)
     
#    Age structure - c -  age binned as 0-1, 1-3, 3+ 
     ddl = add.design.data(true_age.proc,ddl, parameter = "Phi",type="age",
                           bins=c(0,1,3,10),name="c", right=FALSE)
                           
#                  ---List of Model Parameter Designs---
     
# 4.1) Survival
     
     Phi.time.c          = list(formula = ~ -1 +     c  * time          )
     Phi.time.b          = list(formula = ~ -1 +     b  * time          )
     Phi.time.a          = list(formula = ~ -1 +     a  * time          )
     
     Phi.time.a.c        = list(formula = ~ -1 +     c  + time          )
     Phi.time.a.b        = list(formula = ~ -1 +     b  + time          )
     Phi.time.a.a        = list(formula = ~ -1 +     a  + time          )
     
     Phi.c               = list(formula = ~          c                  )
     Phi.b               = list(formula = ~          b                  )
     Phi.a               = list(formula = ~          a                  )
     Phi.d               = list(formula = ~          d                  )
  
# 4.2) Recapture
     
     p.time              = list(formula = ~               time          )
     p.cst               = list(formula = ~  1                          ) 


#                            ---List of models---
#     
# x) Here are MARK model objects that contain MARK input files with PIMS and 
#    design matrix specific to the data and model structure formulas.
     
# ########################################################################### #
     
     Phi.c.t.p.t     = make.mark.model(true_age.proc, ddl, 
                                       parameters = list(Phi = Phi.time.c   , 
                                                         p   = p.time        ))
     
     Phi.c.t.p.cst   = make.mark.model(true_age.proc, ddl, 
                                       parameters = list(Phi = Phi.time.c   , 
                                                         p   = p.cst         ))
     
# ########################################################################### #
  
     Phi.c.a.t.p.t   = make.mark.model(true_age.proc, ddl, 
                                       parameters = list(Phi = Phi.time.a.c , 
                                                         p   = p.time        )
                                       
     Phi.c.p.t       = make.mark.model(true_age.proc, ddl, 
                                       parameters = list(Phi = Phi.c        , 
                                                         p   = p.time        ))
     
# ########################################################################### #
     
     Phi.b.a.t.p.t   = make.mark.model(true_age.proc, ddl, 
                                       parameters = list(Phi = Phi.time.a.b , 
                                                         p   = p.time        ))

     Phi.a.a.t.p.t   = make.mark.model(true_age.proc, ddl, 
                                       parameters = list(Phi = Phi.time.a.a , 
                                                         p   = p.time        ))
     
# ########################################################################### #     
     
     Phi.c.p.cst  = make.mark.model(true_age.proc, ddl, 
                                       parameters = list(Phi = Phi.c         , 
                                                         p   = p.cst         ))
     
     Phi.b.p.cst  = make.mark.model(true_age.proc, ddl, 
                                    parameters = list(Phi = Phi.b            , 
                                                      p   = p.cst            ))
     
     Phi.a.p.cst  = make.mark.model(true_age.proc, ddl, 
                                    parameters = list(Phi = Phi.a            , 
                                                      p   = p.cst            ))
                                                      
# ########################################################################### #


# x) Models Selection: Here, models were created sequentially, first testing 
#                      the best p, and then Phi. The c-hat = 1, 
#                      models were compared by AICc, the model list was called to
#                      assess the produced models.
     
     Phi.c.t.p.t.results    = run.mark.model(Phi.c.t.p.t   )
     Phi.c.t.p.cst.results  = run.mark.model(Phi.c.t.p.cst )
     
     
     Phi.c.a.t.p.t.results  = run.mark.model(Phi.c.a.t.p.t )
     Phi.c.a.p.cst.results  = run.mark.model(Phi.c.p.t     )
     
     Phi.a.a.t.p.t.results  = run.mark.model(Phi.a.a.t.p.t )
     Phi.b.a.t.p.t.results  = run.mark.model(Phi.b.a.t.p.t )
     
     
     Phi.c.p.cst.results = run.mark.model(Phi.c.p.cst )
     Phi.b.p.cst.results = run.mark.model(Phi.b.p.cst )
     Phi.a.p.cst.results = run.mark.model(Phi.a.p.cst )

#  Call all models that have ran

   collect.models(lx       = NULL, 
                  type     = NULL, 
                  table    = TRUE, 
                  adjust   = TRUE,
                  external = FALSE) 
     
# Results of estimates were retrieved by entering a model's results into the
# console, for example: "Phi.a.p.cst.results".
