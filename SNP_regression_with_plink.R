#Setup:
{
library(gdata)                                # load packages
library(tidyverse)
library(splitstackshape)                      # need to read .frq
library(data.table)
library(ggplot2)
library(shiny)
library(foreign)
library(shinydashboard)
library(DT)
}
#Functions:
{
  column_num <- function(name="",dataframe=mydata)
  {
    #This function returns the number of the named column or gives 1
    # e.g.: column_num("DEPR") ---> 27
    if (as.logical(length(which(names(dataframe)==name)))){
      return(which(names(dataframe)==name))}
    else{
      print(paste("Using 1st column. No column found with this name: ", name),quote=FALSE)
      return(1)}
  }
  del_here <- function(name="")
  {
    #This function deletes a named file in current wd, useful for temps
    # e.g.:del_here("temp.txt") --->Deletes it
    file.remove(paste(normalizePath(getwd()), "\\", name, sep=""))
  }
  get_snp_alleles <-function(freq = CLDN3_freq_maf, ped = CLDN3_ped)
  {
    #This function takes the plink created frq and ped data to make an snp allele table
    #example: get_snp_alleles() returns the data table
    myheaders <- freq[,2]   #Take the snp names, they will be the headers
    df<-ped[,7:length(ped)]   #Take the ped data and cut off the first 6 columns (allele data starts at 7)
    dt= data.table(df[,1:(ncol(df)/2)])   #Create a data table with as many rows as people we have and half as many columns as alleles.
    for(i in 1:ncol(dt))    #For loop going through columns
    {
      dt[,(i):=(paste(df[,(i*2-1):(i*2)][[1]],df[,(i*2-1):(i*2)][[2]],sep = ""))]   #Replace every data in i-th column with merged allele data. e.g.:take "A" and "G" from ped data and put "AG" here. 
    }
    names(dt) <- myheaders[[1]]    #Give the SNPs as headers for the data table
    return(dt)
  }
  get_mean_data <- function( j=1,df=myframe, snpd=snp_data, PHENO="SYDEPADD_4a")
  {
    #This function takes a number, ped and snp data and the sought pheno name.
    #It returns the means and sample standard deviation of the pheno category scores,
    #grouped by the alleles of the SNP, specified by the number.
    #It also returns the absolute sum of differences between means and the sum of SSD-s.
    #get_mean_data() returns a list with the data, using the first SNP in default order.
    if(j<=(ncol(snpd)) &j>0)
    {
      mean_data <- drop_na(data.frame(df,snpd),column_num(name=PHENO,dataframe = df),(ncol(df)+j)) %>%
        group_by(drop_na(data.frame(df,snpd),column_num(name=PHENO,dataframe = df),(ncol(df)+j))[(ncol(df)+j)]) %>%
        summarise(mean=mean(.data[[PHENO]]),sd = sd(.data[[PHENO]]),.groups='drop')
      
      absmean <- sum(abs(diff(mean_data[[2]])))+abs(mean_data[[1,2]]-mean_data[[3,2]])
      abssd <- sum(abs(mean_data[[3]]))
      
      return(list("mean_diff"=absmean,"sd_sum"=abssd,"mean_data"=mean_data))
    }
    else if(j==0)
    {
      print(paste("Please refer to the 1st SNP with a '1'."))
    }
    else
    {
      print(paste(j," Is out of bounds for the snp data. The number of SNPs in the data is: ",(ncol(snpd))))
    }
  }
  compare_snp_means_sd<-function(from=1,to=(as.integer(ncol(snpd))),snpd=snp_data,PHENO="SYDEPADD_4a",returndata=FALSE,df=myframe)
  {
    #This function gives back the difference of means and sum of sample standard deviations of the phenotype scores
    #in relation to the alleles of the chosen SNPs. It displays the most promising snp candidates based on these scores.
    if(from<=(as.integer(ncol(snpd))) & from<=to & to<=(as.integer(ncol(snpd))))
    {
      mns <- lapply(1:(to-from+1), c) 
      sds <- lapply(1:(to-from+1), c)
      maxmns <- lapply(1:(to-from+1), c) 
      for(j in 1:(to-from+1))
      {
        mns[(j)]<-(get_mean_data(j=(from+j-1),snpd=snpd,df=df,PHENO = PHENO))$mean_diff
        sds[(j)]<-(get_mean_data(j=(from+j-1),snpd=snpd,df=df,PHENO = PHENO))$sd_sum
        maxmns[(j)]<-((get_mean_data(j=(from+j-1),snpd=snpd,df=df,PHENO = PHENO))$mean_data[,2])
      }
      print(paste("The highest sum of differences between means was: ",max(unlist(mns))," at: ",((which(mns==max(unlist(mns))))[1]+from-1)
                  ,", The lowest sum of SSDs was at: ", (which(sds==min(unlist(sds))))[1]+from-1 ) )
      
      if(returndata)
      {
        return(list(
          "means_index"=which(mns==max(unlist(mns))),
          "SSD_index"=which(sds==min(unlist(sds))),
          "delta_means"=mns,
          "means"=maxmns,
          "SDs"=sds))
      }
    }
    else{print(paste("Snps between ",from," and ",to," is out of bounds for given snp data: ",(as.integer(ncol(snpd)))))}
  }
  citeall<-function()
  {
    #This function prints citations for all the libraries used in the creation of the application
    #citeall() will print to the console.
    print(citation())
    print(citation("gdata"))                              
    print(citation("tidyverse"))
    print(citation("splitstackshape"))                        
    print(citation("data.table"))
    print(citation("ggplot2"))
    print(citation("shiny"))
    print(citation("foreign"))
    print(citation("shinydashboard"))
    print(citation("DT"))
  }
  install-packages()<-function()
  {
    install.packages("gdata")
    install.packages("tidyverse")
    install.packages("splitstackshape")
    install.packages("data.table")
    install.packages("ggplot2")
    install.packages("shiny")
    install.packages("foreign")
    install.packages("shinydashboard")
    install.packages("DT")
  }
}
### APP ############################# 
{
  ui <- dashboardPage(
      # Sidebar ---------------------------------
    dashboardHeader(title = "Basic dashboard"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Data", tabName = "data_tab", icon = icon("database")),
        menuItem("PLINK", tabName = "plink_tab", icon = icon("product-hunt")),
        menuItem("Basic plots", tabName = "basicplots_tab", icon = icon("chart-line")),
        menuItem("SNP plots", tabName = "snpplot_tab", icon = icon("chart-bar"))
      )#<- end of sidebarMenu
    ),#<- end of dashboardSidebar
    dashboardBody(
      tabItems(
        # Data tab ---------------------------------
        tabItem(tabName = "data_tab",h2("Data selection"),
                box(
                  textInput("WD", "Working directory:", getwd(),placeholder = "Set global working directory."),
                  actionButton("SET_WD", "Set"),
                  verbatimTextOutput("LOADED_WD")),
                box(
                  textInput("PHENO_DATA", "Phenotype file:","PHENO_BPMAN.CSV",placeholder = "eg.: Phenotype.CSV"),
                  actionButton("LOAD_PHENO", "Load"),
                  actionButton("MOD_DATA", "Tidy"),
                  verbatimTextOutput("LOADED_PHENO")),
                box(
                  textInput("FAMBIMBED", "Fam/Bim/Bed file:","NM_bpman_1",placeholder = "Please use fam/bim/bed files."),
                  actionButton("LOAD_FAMBIMBED", "Load"),
                  verbatimTextOutput("LOADED_FAMBIMBED")),
                box(
                  textInput("PEDMAP", "Ped/Map file:","CLDN3_10kb",placeholder = "Please use ped/map files."),
                  actionButton("LOAD_PEDMAP", "Load"),
                  verbatimTextOutput("LOADED_PEDMAP")),
                box(
                  textInput("SNPLIST", "Snplist file:","CLDN3_10kb",placeholder = "Please use snplist files."),
                  actionButton("LOAD_SNPLIST", "Load"),
                  verbatimTextOutput("LOADED_SNPLIST")),
                box(
                  textInput("FREQMAF", "Frq (MAF) file:","CLDN3_10kb_freq_MAF",placeholder = "Please use frq files."),
                  actionButton("LOAD_FREQMAF", "Load"),
                  verbatimTextOutput("LOADED_FREQMAF")),
                box(
                  textInput("HARDY", "Hwe (Hardy Weinberg) file:","CLDN3_10kb_hardy",placeholder = "Please use hwe files."),
                  actionButton("LOAD_HARDY", "Load"),
                  verbatimTextOutput("LOADED_HARDY")),
                box(
                  title = "Data viewer", status = "warning",
                  actionButton("SHOW_PHENO", "Phenotype"),
                  actionButton("SHOW_FAMBIMBED", "Fam/bim/bed"),
                  actionButton("SHOW_PEDMAP", "Ped/map"),
                  actionButton("SHOW_SNPLIST", "Snplist"),
                  actionButton("SHOW_FREQMAF", "Freq"),
                  actionButton("SHOW_HARDY", "Hardy")),
                fluidRow(DT::dataTableOutput("DATA_SHOW"), width=600)
                
                
        ),#<- end of data tab 
        # Plink tab ---------------------------------
        tabItem(
          tabName = "plink_tab",h2("PLINK commands"),
          box(title = getwd(),width = 500),
          box(width = 500,
            textInput("PLINK1", "","plink --bfile NM_bpman_1 --chr 7 --from-kb 73173 --to-kb 73194 --recode --out mypedmap",placeholder = "Type in a command."),
            verbatimTextOutput("PLINK1_DEFCOMMAND"),
            actionButton("PLINK1_RUN", "Run")),
          box(width = 500,
              textInput("PLINK2", "","plink --file mypedmap --write-snplist --out mypedmap_snplist",placeholder = "Type in a command."),
              verbatimTextOutput("PLINK2_DEFCOMMAND"),
              actionButton("PLINK2_RUN", "Run")),
          box(width = 500,
              textInput("PLINK3", "","plink --file mypedmap --freq --out mypedmap_freq_MAF",placeholder = "Type in a command."),
              verbatimTextOutput("PLINK3_DEFCOMMAND"),
              actionButton("PLINK3_RUN", "Run")),
          box(width = 500,
              textInput("PLINK4", "","plink --file mypedmap --hardy --out mypedmap_hardy",placeholder = "Type in a command."),
              verbatimTextOutput("PLINK4_DEFCOMMAND"),
              actionButton("PLINK4_RUN", "Run")),
          box(width = 500,status = "warning",
              textInput("PLINK5", "","plink --file mypedmap --pheno PHENO_BPMAN.txt --pheno-name SYDEPADD_4a --covar PHENO_BPMAN.txt --covar-name POPULATION,AGE,GENDER --linear --out mypedmap_depr_linreg",placeholder = "Type in a command."),
              verbatimTextOutput("PLINK5_DEFCOMMAND"),
              selectInput("PLINK5_SELECT", "Regression:",choices=c("linear","logistic")),
              selectInput("PLINK5_SELECT2", "Model:",choices=c("additive","dominant","recessive")),
              actionButton("PLINK5_RUN", "Run"),
              actionButton("PLINK52_RUN", "Add/Dom/Rec")),
        ),#<- end of plink tab 
        # basicplots tab ---------------------------------
        tabItem(
          tabName = "basicplots_tab",h2("Basic plots"),
          fluidRow(
            box(title = "Gender ratio:",
              plotOutput("basicplot1")
            ),
            box(title = "Age groups:",
                plotOutput("basicplot2")
            ),
            box(title = "Age density map:",
                plotOutput("basicplot3"),
                textInput("BASICPLOT3_PHENO", "Phenotype:", "SYDEPADD_4a",placeholder = "e.g.:SYDEPADD_4a"),
                actionButton("BASICPLOT3_SUBMIT","Submit")
            ),
            box(title = "Logistic regression:",
                plotOutput("basicplot4"),
                textInput("BASICPLOT4_PHENO", "Phenotype Y:", "CHS_LP",placeholder = "e.g.:CHS_LP"),
                textInput("BASICPLOT42_PHENO", "Phenotype X:", "DEPR",placeholder = "e.g.:DEPR"),
                actionButton("BASICPLOT4_SUBMIT","Submit")
            ),
            box(title = "Linear regression:",
                plotOutput("basicplot5"),
                textInput("BASICPLOT5_PHENO", "Phenotype Y:", "Recent_LESUM",placeholder = "e.g.:Recent_LESUM"),
                textInput("BASICPLOT52_PHENO", "Phenotype X:", "SYDEPADD_4a",placeholder = "e.g.:SYDEPADD_4a"),
                actionButton("BASICPLOT5_SUBMIT","Submit")
            ),
            box(title = "Gender split plot:",
                plotOutput("basicplot6"),
                textInput("BASICPLOT6_PHENO", "Phenotype Y:", "AGE",placeholder = "e.g.:AGE"),
                textInput("BASICPLOT62_PHENO", "Phenotype X:", "SYDEPADD_4a",placeholder = "e.g.:SYDEPADD_4a"),
                checkboxInput("BASICPLOT6_CHECKBOX", "Density shading", TRUE),
                selectInput("BASICPLOT6_SELECT", "Regression:",choices=c("None","Linear","Logistic")),
                actionButton("BASICPLOT6_SUBMIT","Submit")
            )
          )#<- end of fluidrow
        ),#<- end of basicplots tab 
        # snp plots tab ---------------------------------
        tabItem(
          tabName = "snpplot_tab",h2("SNP plots"),
          fluidRow(
            box(title = "SNP means plot",
                width = 400,
                plotOutput("snpplot1"),
                textInput("SNPPLOT1_PHENO", "Phenotype:", "SYDEPADD_4a",placeholder = "e.g.:SYDEPADD_4a"),
                actionButton("SNPPLOT1_SUBMIT","Submit")
                ),
            box(title="Single SNP plot",
                plotOutput("snpplot2", height = 350),
                textInput("PHENO1", "Phenotype:", "SYDEPADD_4a",placeholder = "write something dummy"),
                numericInput("SNP1", "Index of SNP:", 1, min = 1, max = 21),
                actionButton("submit","Submit"),
                verbatimTextOutput("value")
                )
          )
        )#<- end of snp plot tab
      )#<- end of tabItemS
    )#<- end of dashboard body
  )#<- end of dashboard page
### APP SERVER ############################# 
  server <- function(input, output)
  {
    ##values =================================
    values <- reactiveValues(
      mywd =NULL,
      myval1 = "SYDEPADD_4a",
      mynum1 = 1,
      mypheno = "PHENO_BPMAN.txt",
      myfambimbed = NULL,
      mypedmap = NULL,
      mysnplist =NULL,
      myfreqmaf =NULL,
      myhardy =NULL,
      mydatashow =NULL,
      myi= 0,
      mydata = NULL,
      myframe = NULL,
      basicplot3pheno = "SYDEPADD_4a",
      basicplot4pheno = "CHS_LP",
      basicplot42pheno= "DEPR",
      basicplot5pheno = "Recent_LESUM",
      basicplot52pheno= "SYDEPADD_4a",
      basicplot6pheno = "AGE",
      basicplot62pheno= "SYDEPADD_4a",
      basicplot6checkbox =TRUE,
      basicplot6select= "NULL",
      mysnp_data=NULL,
      myped_data=NULL,
      myfreq_data=NULL,
      mysnplist_data=NULL,
      snpplot1_pheno="SYDEPADD_4a",
      myname="mypedmap"
      )#<- end of values
    ##DATA Tab: =================================
      #Load buttons: ---------------------------------
      observeEvent(input$SET_WD, {
        if(dir.exists(isolate(input$WD)))
        {
          values$mywd <- isolate(input$WD)
          setwd(values$mywd)
          print(paste("WD set to:",values$mywd))
        }
        else
        {
          print("Directory does not exist.")
        }
      })
      observeEvent(input$LOAD_PHENO, {
        if(file.exists(isolate(input$PHENO_DATA)))
        {
          values$mypheno <- isolate(input$PHENO_DATA)
          values$mydata <- read.csv(isolate(input$PHENO_DATA),sep=";",na="-9", stringsAsFactors = F, encoding = "UTF-8")
          print(paste(isolate(input$PHENO_DATA), "loaded."))
          values$myframe <- values$mydata
        }
        else
        {
          print(paste("Could not load:",isolate(input$PHENO_DATA)))
        }
  
      })
      observeEvent(input$MOD_DATA, { 
        values$mydata[ ,column_num("GENDER",dataframe = values$mydata)][values$mydata[ ,column_num("GENDER",dataframe = values$mydata)] == 1] <-"Male";
        values$mydata[ ,column_num("GENDER",dataframe = values$mydata)][values$mydata[ ,column_num("GENDER",dataframe = values$mydata)] == 2] <-"Female";
        values$mydata[ ,column_num("DEPR",dataframe = values$mydata)] <- as.logical(values$mydata[ ,column_num("DEPR",dataframe = values$mydata)]-1);
        values$mydata[ ,column_num("Age_range",dataframe = values$mydata)] <-case_when(values$mydata[ ,column_num("Age_range",dataframe = values$mydata)]==1 ~ "<20",values$mydata[ ,column_num("Age_range",dataframe = values$mydata)]==2 ~ "21-30",values$mydata[ ,column_num("Age_range",dataframe = values$mydata)]==3 ~ "31-40",values$mydata[ ,column_num("Age_range",dataframe = values$mydata)]==4 ~ "41-50",values$mydata[ ,column_num("Age_range",dataframe = values$mydata)]==5 ~ "51-60",values$mydata[ ,column_num("Age_range",dataframe = values$mydata)]==6 ~ "61<");
        values$myframe <- values$mydata
        if(!is.null(values$myped_data) && !is.null(values$myfreq_data)){
          values$mysnp_data <-na_if(   get_snp_alleles(freq = values$myfreq_data,ped = values$myped_data),   "00")
          print("Made snp data from ped and freq data.")
        }
        else{print("Load ped and freq data to create snp data file")}
      })
      observeEvent(input$LOAD_FAMBIMBED, {
        if (any(file.exists(paste(isolate(input$FAMBIMBED),".fam",sep = ""),paste(isolate(input$FAMBIMBED),".bim",sep = ""),paste(isolate(input$FAMBIMBED),".bed",sep = ""),isolate(input$FAMBIMBED) )))
        {
          values$myfambimbed <- isolate(input$FAMBIMBED)
          if(any(str_detect(isolate(input$FAMBIMBED),"\\.fam"), str_detect(isolate(input$FAMBIMBED),"\\.bim"),str_detect(isolate(input$FAMBIMBED),"\\.bed")) )
          {
            values$myfambimbed <- gsub('.{4}$','',values$myfambimbed)
          }
          print(paste(values$mywd,"/",values$myfambimbed," loaded.",sep = ""))
          
        }
        else
        {
          print(paste(isolate(input$FAMBIMBED),"fam/bim/bed file not found in working directory."))
        }
       
        
      })
      observeEvent(input$LOAD_PEDMAP, {
        if (any(file.exists(paste(isolate(input$PEDMAP),".ped",sep = ""),paste(isolate(input$PEDMAP),".map",sep = ""),isolate(input$PEDMAP) )))
        {
          values$mypedmap <- isolate(input$PEDMAP)
          if(any(str_detect(isolate(input$PEDMAP),"\\.ped"), str_detect(isolate(input$PEDMAP),"\\.map")) )
          {
            values$mypedmap <- gsub('.{4}$','',values$mypedmap)
          }
          values$myped_data = read.table(paste(values$mypedmap,".ped",sep=""), sep ="\n", header = FALSE)
          values$myped_data = cSplit(values$myped_data, 1:ncol(values$myped_data), sep=" ", stripWhite=TRUE, type.convert=FALSE)
          print(paste(values$mywd,"/",values$mypedmap," loaded.",sep = ""))
        }
        else
        {
          print(paste(isolate(input$PEDMAP),"ped/map file not found in working directory."))
        }
        
      })
      observeEvent(input$LOAD_SNPLIST, {
        if (any(file.exists(paste(isolate(input$SNPLIST),".snplist",sep = ""),isolate(input$SNPLIST) )))
        {
          values$mysnplist <- isolate(input$SNPLIST)
          if(str_detect(isolate(input$SNPLIST),"\\.snplist") ) 
          {
            values$mysnplist <- gsub('.{8}$','',values$mysnplist)
          }
          values$mysnplist_data <- read.table(paste(values$mysnplist,".snplist",sep=""))
          print(paste(values$mywd,"/",values$mysnplist," loaded.",sep = ""))
          
        }
        else
        {
          print(paste(isolate(input$SNPLIST),"snplist file not found in working directory."))
        }
        
        
      })
      observeEvent(input$LOAD_FREQMAF, {
        if (any(file.exists(paste(isolate(input$FREQMAF),".frq",sep = ""),isolate(input$FREQMAF) )))
        {
          values$myfreqmaf <- isolate(input$FREQMAF)
          if(str_detect(isolate(input$FREQMAF),"\\.frq")  )
          {
            values$myfreqmaf <- gsub('.{4}$','',values$myfreqmaf)
          }

          values$myfreq_data = read.table(paste(values$myfreqmaf,".frq",sep=""), sep ="\n",skip=1, header = FALSE)
          values$myfreq_data = cSplit(values$myfreq_data, 1:ncol(values$myfreq_data), sep=" ", stripWhite=TRUE, type.convert=FALSE)
          
          print(paste(values$mywd,"/",values$myfreqmaf," loaded.",sep = ""))
          
        }
        else
        {
          print(paste(isolate(input$FREQMAF),"frq file not found in working directory."))
        }
        
        
      })
      observeEvent(input$LOAD_HARDY, {
        if (any(file.exists(paste(isolate(input$HARDY),".hwe",sep = ""),isolate(input$HARDY) )))
        {
          values$myhardy <- isolate(input$HARDY)
          if(str_detect(isolate(input$HARDY),"\\.hwe") )
          {
            values$myhardy <- gsub('.{4}$','',values$myhardy)
          }
          print(paste(values$mywd,"/",values$myhardy," loaded.",sep = ""))
          
        }
        else
        {
          print(paste(isolate(input$HARDY),"hwe file not found in working directory."))
        }
        
        
      })
      #Data viewer buttons: ---------------------------------
      observeEvent(input$SHOW_PHENO, {
        if( is.null(values$mydata) )
        {
          print("Phenotype data not loaded.")
        }
        else
        {
          values$mydatashow = datatable(values$mydata, options=list(scrollX=TRUE,  ordering=FALSE), rownames = FALSE)
        }
        values$myi=0
      })
      observeEvent(input$SHOW_FAMBIMBED, {
        if( is.null(values$myfambimbed) )
        {
          print("Fam/bim/bed data not loaded.")
        }
        else
        {
          if(values$myi == 3)
          {
            print("Reading .bim file. Please wait...")
            values$mydatashow = read.table(paste(values$myfambimbed,".bim",sep = ""),header = FALSE)
            print("Displaying .bim file.")
          }
          else
          {
            print("Displaying .fam file.")
            values$mydatashow = read.table(paste(values$myfambimbed,".fam",sep = ""),header = FALSE)
            if(values$myi == 2)
            {
              print("Displaying large .bim files can take a long time. Clicking again will try displaying the .bim file.")
            }
          }
          values$myi=values$myi +1
        }
      }) # Reading bim takes a long time. Bed is not to be read (incomprehensible).
      observeEvent(input$SHOW_PEDMAP, {
        if( is.null(values$mypedmap) )
        {
          print("Ped/map data not loaded.")
        }
        else
        {
          if(values$myi %% 2)
          {
          print("Displaying .ped file.")
          values$mydatashow = datatable(read.table(paste(values$mypedmap,".ped",sep = ""),header = FALSE), options=list(scrollX=TRUE,  ordering=FALSE), rownames = FALSE)
          }
          else
          {
            print("Displaying .map file.")
            values$mydatashow = read.table(paste(values$mypedmap,".map",sep = ""),header = FALSE)
          }
          if(values$myi == 0){values$myi = 1}
          else{values$myi = 0}
        }
      })
      observeEvent(input$SHOW_SNPLIST, {
        if( is.null(values$mysnplist) )
        {
          print("Snplist data not loaded.")
        }
        else
        {
          values$mydatashow = read.table(paste(values$mysnplist,".snplist",sep = ""),header = FALSE)
        }
        values$myi=0
      })
      observeEvent(input$SHOW_FREQMAF, {
        if( is.null(values$myfreqmaf) )
        {
          print("Freq data not loaded.")
        }
        else
        {
          values$mydatashow = read.table(paste(values$myfreqmaf,".frq",sep = ""),header = TRUE)
        }
        values$myi=0
      })
      observeEvent(input$SHOW_HARDY, {
        if( is.null(values$myhardy) )
        {
          print("Hardy data not loaded.")
        }
        else
        {
          values$mydatashow = read.table(paste(values$myhardy,".hwe",sep = ""),header = TRUE)
        }
        values$myi=0
      })
      #Texts: ---------------------------------
      output$LOADED_WD <- renderText({paste("Current:",values$mywd)})
      output$LOADED_PHENO <- renderText({paste("Current:",values$mypheno)})
      output$LOADED_FAMBIMBED <- renderText({paste("Current:",values$myfambimbed)})
      output$LOADED_PEDMAP <- renderText({paste("Current:",values$mypedmap)})
      output$LOADED_FREQMAF<- renderText({paste("Current:",values$myfreqmaf)})
      output$LOADED_SNPLIST <- renderText({paste("Current:",values$mysnplist)})
      output$LOADED_HARDY <- renderText({paste("Current:",values$myhardy)})
      #Data viewer ---------------------------------
      output$DATA_SHOW <- renderDataTable({values$mydatashow})
      
    ##PLINK tab =================================
      #Buttons ---------------------------------
      observeEvent(input$PLINK1_RUN, {
      shell(input$PLINK1)
      if (gregexpr("--out",input$PLINK1,fixed = TRUE)>0)
      {
        values$myname<-substring(input$PLINK1, unlist(gregexpr("--out",input$PLINK1))+6 )
      }
      })
      observeEvent(input$PLINK2_RUN, {
        shell(input$PLINK2)
      })
      observeEvent(input$PLINK3_RUN, {
        shell(input$PLINK3)
      })
      observeEvent(input$PLINK4_RUN, {
        shell(input$PLINK4)
      })
      observeEvent(input$PLINK5_RUN, {
        shell(input$PLINK5)
      })
      observeEvent(input$PLINK52_RUN, {
        myinput=input$PLINK5
        if(any(str_detect(myinput," --dominant"),str_detect(myinput," --recessive") ))
        {
          myinput=gsub(" --dominant","",myinput,fixed = TRUE)
          myinput=gsub(" --recessive","",myinput,fixed = TRUE)
        }
        add_command=gsub(substr(myinput,regexpr("--out",myinput,fixed=TRUE)[1],nchar(myinput)),paste("--out ",substr(myinput,regexpr("--out",myinput,fixed=TRUE)[1]+6,nchar(myinput)),"_additive",sep=""),myinput)
        dom_command=gsub(substr(myinput,regexpr("--out",myinput,fixed=TRUE)[1],nchar(myinput)),paste("--out ",substr(myinput,regexpr("--out",myinput,fixed=TRUE)[1]+6,nchar(myinput)),"_dominant",sep=""),myinput)
        dom_command=paste(substr(dom_command, 1, regexpr("--out",myinput,fixed=TRUE)[1]-1), "--dominant ", substr(dom_command, regexpr("--out",myinput,fixed=TRUE)[1], nchar(dom_command)), sep = "")
        rec_command=gsub(substr(myinput,regexpr("--out",myinput,fixed=TRUE)[1],nchar(myinput)),paste("--out ",substr(myinput,regexpr("--out",myinput,fixed=TRUE)[1]+6,nchar(myinput)),"_recessive",sep=""),myinput)
        rec_command=paste(substr(rec_command, 1, regexpr("--out",myinput,fixed=TRUE)[1]-1), "--recessive ", substr(rec_command, regexpr("--out",myinput,fixed=TRUE)[1], nchar(rec_command)), sep = "")
        shell(add_command)
        shell(dom_command)
        shell(rec_command)
      })
      #Texts ---------------------------------
      output$PLINK1_DEFCOMMAND <- renderText({paste("plink --bfile NM_bpman_1 --chr 7 --from-kb 73173 --to-kb 73194 --recode --out",values$myname)})
      output$PLINK2_DEFCOMMAND <- renderText({paste("plink --file ",values$myname," --write-snplist --out ",values$myname,"_snplist",sep = "")})
      output$PLINK3_DEFCOMMAND <- renderText({paste("plink --file ",values$myname," --freq --out ",values$myname,"_freq_MAF",sep = "")})
      output$PLINK4_DEFCOMMAND <- renderText({paste("plink --file ",values$myname," --hardy --out ",values$myname,"_hardy",sep = "")})
      output$PLINK5_DEFCOMMAND <- renderText({if (input$PLINK5_SELECT2=="additive"){a=""}else{a=paste(" --",input$PLINK5_SELECT2,sep="")}
                                              paste("plink --file ",values$myname," --pheno ",gsub('.{4}$','.txt',values$mypheno)," --pheno-name SYDEPADD_4a --covar ",gsub('.{4}$','.txt',values$mypheno)," --covar-name POPULATION,AGE,GENDER --",input$PLINK5_SELECT,a," --out ",values$myname,"_",substr(input$PLINK5_SELECT2,1,3),"_reg",sep = "")})
      
    ##Basic plot tab =================================
      #Genders plot: ---------------------------------
      output$basicplot1 <- renderPlot( { 
        ggplot(drop_na(values$myframe,GENDER), aes(x=as.factor(GENDER)),t)+ geom_bar() +labs(x="Genders")
      } )
      #Age groups plot: ---------------------------------
      output$basicplot2 <- renderPlot( { 
        ggplot(drop_na(data=values$myframe,Age_range), aes(x=as.factor(Age_range)),t)+ geom_bar() +labs(x="Age groups")
      } )
      #Density plot: ---------------------------------
      observeEvent(input$BASICPLOT3_SUBMIT, {
          values$basicplot3pheno<-input$BASICPLOT3_PHENO
      })
      output$basicplot3 <- renderPlot( { 
        ggplot(drop_na(values$myframe,AGE,.data[[values$basicplot3pheno]]), aes(x=drop_na(values$myframe,AGE,.data[[values$basicplot3pheno]])[ ,column_num(values$basicplot3pheno,dataframe = values$mydata)],y=drop_na(values$myframe,AGE,.data[[values$basicplot3pheno]])[ ,column_num("AGE",dataframe = values$mydata)])) + 
          stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
          scale_fill_distiller(palette=8, direction=1) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          theme(legend.position='none') +
          geom_point(aes(shape = "cross"), size = .5,position = "jitter") + 
          labs(x=values$basicplot3pheno, y="Age")
      } )
      #logistic regression plot: ---------------------------------
      observeEvent(input$BASICPLOT4_SUBMIT, {
        values$basicplot4pheno<-input$BASICPLOT4_PHENO
        values$basicplot42pheno<-input$BASICPLOT42_PHENO
      })
      output$basicplot4 <- renderPlot( { 
        ggplot(drop_na(values$myframe,values$basicplot4pheno,values$basicplot42pheno), aes(y=(as.integer(drop_na(values$myframe,values$basicplot4pheno,values$basicplot42pheno)[ ,column_num(values$basicplot42pheno,dataframe = values$mydata)])),x=drop_na(values$myframe,values$basicplot4pheno,values$basicplot42pheno)[ ,column_num(values$basicplot4pheno,dataframe = values$mydata)])) + 
          geom_point(aes(shape = "cross"), size = .5) +
          geom_jitter(width = 0,height = 0.1,size=0.5) +
          geom_smooth(method = glm,method.args=list(family="binomial"),formula = 'y~x') +
          labs(y=values$basicplot42pheno,x=values$basicplot4pheno) +
          theme(legend.position='none') 
      } )
      #linear regression plot: ---------------------------------
      observeEvent(input$BASICPLOT5_SUBMIT, {
        values$basicplot5pheno<-input$BASICPLOT5_PHENO
        values$basicplot52pheno<-input$BASICPLOT52_PHENO
      })
      output$basicplot5 <- renderPlot( { 
        ggplot(drop_na(values$myframe,values$basicplot5pheno,values$basicplot52pheno),aes(x=(drop_na(values$myframe,values$basicplot5pheno,values$basicplot52pheno)[ ,column_num(values$basicplot5pheno,dataframe = values$mydata)]),y=drop_na(values$myframe,values$basicplot5pheno,values$basicplot52pheno)[ ,column_num(values$basicplot52pheno,dataframe = values$mydata)]))+
            labs(x=values$basicplot52pheno,y=values$basicplot5pheno)+
            geom_point(aes(shape = "cross"), size = .5, position = "jitter") +
            geom_smooth(method = lm,formula='y ~ x')+           
            theme(legend.position='none')
      } )
      #gender split: ---------------------------------
      observeEvent(input$BASICPLOT6_SUBMIT, {
        values$basicplot6pheno<-input$BASICPLOT6_PHENO
        values$basicplot62pheno<-input$BASICPLOT62_PHENO
        values$basicplot6checkbox<-input$BASICPLOT6_CHECKBOX
        values$basicplot6select<-input$BASICPLOT6_SELECT
        
      })
      output$basicplot6 <- renderPlot({
        if(values$basicplot6checkbox)
        {
          densityshade <- stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)
        }
        else
        {
          densityshade <- NULL
        }
        if(values$basicplot6select=="Linear")
        {
          regressionmen <- geom_smooth(data = drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Male"),],AGE,SYDEPADD_4a),aes(x=(drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Male"),],AGE,SYDEPADD_4a)[ ,column_num("AGE",dataframe = values$myframe)]),y=drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Male"),],AGE,SYDEPADD_4a)[ ,column_num("SYDEPADD_4a",dataframe = values$myframe)]),method = lm,formula = 'y~x', color="red")
          regressionwomen <- geom_smooth(data = drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Female"),],AGE,SYDEPADD_4a),aes(x=(drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Female"),],AGE,SYDEPADD_4a)[ ,column_num("AGE",dataframe = values$myframe)]),y=drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Female"),],AGE,SYDEPADD_4a)[ ,column_num("SYDEPADD_4a",dataframe = values$myframe)]),method = lm,formula = 'y~x', color="blue")
        }
        else if(values$basicplot6select=="Logistic")
        {
          regressionmen <- geom_smooth(data = drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Male"),],AGE,SYDEPADD_4a),aes(x=(drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Male"),],AGE,SYDEPADD_4a)[ ,column_num("AGE",dataframe = values$myframe)]),y=drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Male"),],AGE,SYDEPADD_4a)[ ,column_num("SYDEPADD_4a",dataframe = values$myframe)]),method = glm,formula = 'y~log(x)', color="red")
          regressionwomen <- geom_smooth(data = drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Female"),],AGE,SYDEPADD_4a),aes(x=(drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Female"),],AGE,SYDEPADD_4a)[ ,column_num("AGE",dataframe = values$myframe)]),y=drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Female"),],AGE,SYDEPADD_4a)[ ,column_num("SYDEPADD_4a",dataframe = values$myframe)]),method = glm,formula = 'y~log(x)', color="blue")
        }
        else
        {
          regressionmen <- NULL
          regressionwomen <- NULL
        }
        ggplot(drop_na(values$myframe,values$basicplot6pheno,values$basicplot62pheno),aes(x=(drop_na(values$myframe,values$basicplot6pheno,values$basicplot62pheno)[ ,column_num(values$basicplot6pheno,dataframe = values$myframe)]),y=drop_na(values$myframe,values$basicplot6pheno,values$basicplot62pheno)[ ,column_num(values$basicplot62pheno,dataframe = values$myframe)]))+
          labs(x=values$BASICPLOT6_PHENO,y=values$BASICPLOT62_PHENO) +
          densityshade +
          scale_fill_distiller(palette=6, direction=1) +
          scale_shape_identity() +
          regressionmen +
          regressionwomen + 
          geom_point(data=drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Male"),],values$basicplot6pheno,values$basicplot62pheno),aes(x=(drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Male"),],values$basicplot6pheno,values$basicplot62pheno)[ ,column_num(values$basicplot6pheno,dataframe = values$myframe)]),y=drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Male"),],values$basicplot6pheno,values$basicplot62pheno)[ ,column_num(values$basicplot62pheno,dataframe = values$myframe)]),shape = "triangle filled", position="jitter",size = .7, fill="red") +
          geom_point(data=drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Female"),],values$basicplot6pheno,values$basicplot62pheno),aes(x=(drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Female"),],values$basicplot6pheno,values$basicplot62pheno)[ ,column_num(values$basicplot6pheno,dataframe = values$myframe)]),y=drop_na(values$myframe[which(values$myframe[,column_num("GENDER",dataframe = values$myframe)]=="Female"),],values$basicplot6pheno,values$basicplot62pheno)[ ,column_num(values$basicplot62pheno,dataframe = values$myframe)]),position = "jitter",shape="circle filled", fill="blue", size=.7) +
          theme(legend.position='none')
          
      })
    ##SNP plot tab =================================
      #All SNPs plot: ---------------------------------
      observeEvent(input$SNPPLOT1_SUBMIT, {
        values$snpplot1_pheno<-input$SNPPLOT1_PHENO
      })
      output$snpplot1 <- renderPlot( { 
        if(is.null(values$mysnplist_data)){print("No snp list loaded.")}
        else
        {
          
          mns<-compare_snp_means_sd(returndata = TRUE,PHENO = values$snpplot1_pheno,snpd = values$mysnp_data,df = values$myframe)$means
          
          a<-data.frame(unlist(mns))
          b<-rep(sample(x=colors(),size=length(unlist(values$mysnplist_data)),replace = FALSE),each=3)
          c<-cbind(a,b)
         
          ggplot(c,aes(x=1:nrow(a),y=a[,1],fill=b) )+
            geom_bar(stat = "identity", width=0.5)+
            geom_hline(yintercept = mean(drop_na(values$myframe,values$snpplot1_pheno)[,column_num(values$snpplot1_pheno,dataframe = values$myframe)]), color="blue")+
            coord_cartesian(ylim=c(min(a)-0.1,max(a)+0.1))+
            labs(y=paste("Mean ",values$snpplot1_pheno), x="SNP index",fill="SNP")+
            scale_x_continuous(breaks = seq(1,length(b),by=3),labels=paste(seq(1,length(unlist(values$mysnplist_data))),unlist(values$mysnplist_data),sep='; '))+
            theme(legend.position = "none",axis.text.x = element_text(size = 10,angle = 90),legend.key.size = unit(0.2,"cm"))
        }
      } )
      #One SNP plot: ---------------------------------
      observe({
        if(input$submit >0)
        {values$myval1<- isolate(input$PHENO1)
        values$mynum1<-isolate(input$SNP1)}
      })
      output$snpplot2 <- renderPlot( { 
          #if myval>0 and <maxsnps then show
        if(is.null(values$mysnplist_data)){print("No snp list loaded.")}
        else
        if(values$mynum1>0 && values$mynum1<=length(unlist(values$mysnplist_data)) )
        {
          PHENO<- values$myval1
          ggplot(drop_na(data.frame(values$myframe,values$mysnp_data),PHENO,(ncol(values$myframe)+values$mynum1)), aes(x=drop_na(data.frame(values$myframe,values$mysnp_data),PHENO,(ncol(values$myframe)+values$mynum1))[ ,  (ncol(values$myframe)+values$mynum1)], y=drop_na(data.frame(values$myframe,values$mysnp_data),PHENO,(ncol(values$myframe)+values$mynum1))[ ,column_num(PHENO,dataframe = values$mydata)]))+
          geom_point(aes(shape = "cross"), size = .5,position = "jitter")+
          stat_summary(fun.data=mean_sdl,fun.args = list(mult=1), geom="pointrange",shape=18,size=1,color="red")+
          theme(legend.position='none') +
          geom_hline(yintercept = mean(drop_na(values$myframe,PHENO)[ ,column_num(PHENO,dataframe = values$mydata)]), color="blue")+
          labs(y=paste(PHENO," scores"), x=colnames(values$mysnp_data)[(values$mynum1)])
        }
        else {print(paste("There are only",length(unlist(values$mysnplist_data)),"SNPs loaded.",values$mynum1,"is out of bounds."))}
        } )
  }#<- end of the server function
  shinyApp(ui, server)
}#<- end of APP ----------------------------------------------------------------------------------------------------
