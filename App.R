library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(shinythemes)
library(reshape2)
library(lubridate)
library(shinydashboard)
library(shinyWidgets)
library(shinydashboardPlus)
library(shinyBS)
library(shinyjs)


CSS <- "
.qtip-big { 
  font-size: 15px;
  line-height: 18px;
  white-space: wrap;
  word-spacing: 1px;
  background: #F0FFFF;
  border-color: #F0FFFF;
  background-color:#F0FFFF;
}
"
createDiv <- function(n, prefixID){
  sprintf(paste(
    "for(var i = 1; i <= %d; i++){",
    "  var div;",
    sprintf("  var id = '%s-' + i.toString();", prefixID),
    "  if(document.getElementById(id) === null){",
    "    div = document.createElement('div');",
    "    div.setAttribute('id', id);",
    "    div.setAttribute('class', 'qtip-big');",
    "    div.style.display = 'none';",
    "    document.body.appendChild(div);",
    "  }",
    "}",
    sep = "\n"
  ), n)
}




fillDiv <- function(dat, i,infodat, prefixID){
  x <- dat[[i]]
  y<-infodat[1,i+1]
 
  sprintf(paste(
    "var div = document.getElementById('%s-%d');",
    "var html = '<b>  </b> %s</br>';",
    "div.innerHTML = html;",
    sep = "\n"
  ), 
  prefixID,
  i,
  y
  
)
  
}


tooltips <- function(n, prefixID){
  settings <- sprintf(paste(
    "{",
    "  overwrite: true,",
    "  content: {",
    sprintf("    text: $('#%s-%%s').clone()", prefixID),
    "  },",
    "  show: {",
    "    ready: false",
    "  },",
    "  position: {",
    "    my: 'bottom %%s',",
    "    at: 'top center'",
    "  },",
    "  style: {",
    "    classes: 'qtip-youtube'",
    "  }",
    "}",
    sep = "\n"
  ), 1:n)
  settings <- sprintf(settings, ifelse(1:n > n/2, "right", "left"))
  sprintf("var tooltips = [%s];", paste0(settings, collapse=","))
}

##### App User Interface #####

###### Header Content ######

ui<-fluidPage(title ="ASF app",theme = shinytheme("spacelab"),
              # setBackgroundColor(color = "ghostwhite"),
              useShinydashboard(),
              
              tags$style(HTML("
    .tabbable > .nav > li > a                  {background-color: sage;  color:black}
    .tabbable > .nav > li > a[data-value='t1'] {background-color: red;   color:white}
    .tabbable > .nav > li > a[data-value='t2'] {background-color: blue;  color:white}
    .tabbable > .nav > li > a[data-value='t3'] {background-color: green; color:white}
    .tabbable > .nav > li[class=active]    > a {background-color: black; color:white}
  ")),
              
              tags$head(tags$style(".shiny-notification {position: fixed; top: 60% ;left: 50%}")),            
    
 fluidRow(
   titlePanel( div(column(width = 8,h4(div(HTML("<b>ASF Shiny App (ver. L2)</b>"),style = "color:#556B2F")),
                        h4(div(HTML("<b>What is the best pre-movement testing approach to use during a disease outbreak?</b>"),style = "color:maroon"),br(),
                                       div(HTML("<b>Use this interactive application to explore how ASF transmission dynamics impact detection via premovement surveillance</b>"),style = "color:darkblue"))), 
                 
                 column(width = 4,br(),
                        tags$img(src = "ulog.jpg",height="120%", width="60%", align="right"))),
             windowTitle="MyPage"
 )),
 
 
 ###### Side bar panel UI ######
 
 
  fluidRow (
    column(4,offset=0,#entire side bar panel
           br(),
           br(),
          tags$style(HTML(".box.box-solid.box-primary>.box-header {
                                color:#FFFFFF;
                                background-color:#778899;}

                                .box.box-solid.box-primary{
                                border-bottom-color:#778899;
                                border-left-color:#778899;
                                border-right-color:#778899;
                                border-top-color:#778899;
                                }")),
          tags$style(HTML("


.box.box-solid.box-danger>.box-header {
  color:#8D6E63;
  background:#8D6E63
                    }

.box.box-solid.box-danger{
border-bottom-color:#8D6E63;
border-left-color:#8D6E63;
border-right-color:#8D6E63;
border-top-color:#8D6E63;
}

.box.box-danger>.box-header {
  color:#8D6E63;
  background:#8D6E63
                    }

.box.box-danger{
border-bottom-color:#8D6E63;
border-left-color:#8D6E63;
border-right-color:#8D6E63;
border-top-color:#8D6E63;
}

                                    ")),
 
          tags$style(HTML("


.box.box-solid.box-info>.box-header {
  color:#4c5b85;
  background:#4c5b85;
                    }

.box.box-solid.box-info{
border-bottom-color:#4c5b85;
border-left-color:#4c5b85;
border-right-color:#4c5b85;
border-top-color:#4c5b85;
}

.box.box-info>.box-header {
  color:#4c5b85;
  background:#4c5b85;
                    }

.box.box-info{
border-bottom-color:#4c5b85;
border-left-color:#4c5b85;
border-right-color:#4c5b85;
border-top-color:#4c5b85;
}

                                    ")),         
          
          tags$head(
            tags$link(rel = "stylesheet", href = "jquery.qtip.min.css"),
            tags$script(src = "jquery.qtip.min.js"),
            tags$style(CSS)
          ),
          useShinyjs(),
          
           
 #######Input transmission dynamics parameters side bar panel UI######
           box(solidHeader = T, collapsible = T, title = "Input transmission dynamics parameters", status = "primary",width=12,collapsed = TRUE,     
           
             fluidRow (
            
               fluidRow(
                 
                 column(width=1),
               column(width=8,offset = 0 ,
                      tipify(selectInput("inp_contact_scenario", (div(HTML("How fast is disease spreading within and between pens(contact rates)?"),style = "color:black")),
                                  list("Slow"=1,"Fast"=2,"Medium"=3),selected =3 
                      ), placement="right",title = "The within and between pen adequate contact rates determine the speed of ASF spread within the barn", trigger = "hover")       
               ),
               br(),
               br(),
               column(width=3,offset =0 ,
                      tipify(checkboxInput("chbox_det", label=div(HTML("View <br/> details")), value = FALSE, width = NULL),
                             placement="bottom",title = "Check to view contact rate distribution details", trigger = "hover")
                   
               )
               
                ),#end column,
               
               
               conditionalPanel(
                 condition = "input.inp_contact_scenario == '4'||input.chbox_det",
               fluidRow(
               column(width=1),
               column(width=8,offset = 0,
                      h5(div(HTML("Input the within pen contact rate <br> distribution values"))),
                      fluidRow(
                        
                       column(width=4,offset =0 , 
                              
                              tipify(textInput("txtinp_betamin", "Lower", value = "1", width = NULL,
                                               placeholder = NULL),placement="bottom",title = "Enter the lower value of<br> the adequate contact rate", trigger = "hover"),
                      # tipify(textInput("txtinp_betamin", "Lower", value = "2", width = 3,
                    ),
                      column(width=4,offset = 0 ,
                      tipify(textInput("txtinp_betamode", "Central", value = "1.64", width = NULL,
                                       placeholder = NULL),placement="bottom",title = "Enter the most likely value of<br> the adequate contact rate", trigger = "hover"),
                      ),
                      column(width=4,offset = 0 ,
                      tipify(textInput("txtinp_betahigh", "High", value = "2.74", width = NULL,
                                       placeholder = NULL),placement="bottom",title = "Enter the upper value of<br> the adequate contact rate", trigger = "hover")
                       
                       # tipify(textInput("txtinp_betamode", "On click"), "Hello again! This is a click-able pop-up", 
                      )),
                    
                    #between pen
                   
                   
                    h5(div(HTML("Input the between pen contact rate <br> distribution values"))),
                    fluidRow(
                      
                      column(width=4,offset =0 , 
                             
                             tipify(textInput("txtinp_bpbetamin", "Lower", value = "1", width = NULL,
                                              placeholder = NULL),placement="bottom",title = "Enter the lower value of<br> the adequate contact rate", trigger = "hover"),
                             # tipify(textInput("txtinp_betamin", "Lower", value = "2", width = 3,
                      ),
                      column(width=4,offset = 0 ,
                             tipify(textInput("txtinp_bpbetahigh", "High", value = "2.74", width = NULL,
                                              placeholder = NULL),placement="bottom",title = "Enter the upper value of<br> the adequate contact rate", trigger = "hover")
                             
                             # tipify(textInput("txtinp_betamode", "On click"), "Hello again! This is a click-able pop-up", 
                      ),
                    
                      column(width=4,offset = 0 ,
                             shinyjs::hidden(textInput("txtinp_bpbetamode", "Central", value = "1.64", width = NULL,
                                                              placeholder = NULL)),
                      ))
                    
                    
               ),
              
               )),
              
                
            
              fluidRow(
                column(width=1),
                column(width=5,
                       
                       tipify(selectInput("inp_ASFstrain", div(HTML("What is the ASFV Strain?"),style = "color:black"),
                                   list("Moderately virulent"=1,"Highly virulent Georgia hu"=3,"Highly virulent Olesen long"=4),selected ="Moderately_virulent" 
                       ),placement="right",title = "This is the scenario for input parameters related to strain charecteristics such as disease state durations and disease mortality", trigger = "hover")
                       
                ),
                
                
              column(width=5, 
               # tags$head(
               #   tags$style(type="text/css", "#inline label{ margin-left: 10px;display: table-cell; text-align: center; vertical-align: middle; } 
               #  #inline .form-group { display: table-row;}")
               # ),
               selectInput("Inp_pensize", div(HTML("How many pigs per pen?"),style="color:black"),
                           list("Small 40 pigs "=40,"Large 120 pigs"=120),selected = 40)
               
              ) 
              
              
              
              ),
              
              fluidRow(
                column(width=1),
                column(width=5,
                       selectInput("Inp_barnsize", div(HTML("How many pigs <br> per barn?"),style="color:black"),
                                   list("Small 1200 pigs"=1200,"Medium 2400 pigs"=2400),selected = 1200)
               ) ,
               
               column(width=5,
                     conditionalPanel(condition = "2>1",
                      tags$div(id = "inline",tipify(numericInput("txtinp_numiter", div(HTML("Number of <br>iterations(runs)  "),style="color:black"), value = 500,max=2000,min=1),placement="bottom",title = "Enter the number of simulation iterations<br>", trigger = "hover")))
               ) ),
              
               
               
             
            
              
              #end fluid row
              fluidRow(
                column(width=1),
              column(width=3, 
                     
                     
                     
                     actionButton("do", HTML("Simulate <br> Transmission"))
                     
                     
              ),
              
              column(width=4, 
                     
                     
                     
                     
                     actionButton("oth_inp", HTML("Change advanced <br> options")),
                     bsModal(id = "clinmodal", title = "Normal Morbidity and Mortality",trigger =  "oth_inp", size = "large",
                          #testing radio buttons
                          
                          
                          
                            
                             fluidRow(
                               
                               column(width = 3,
                                      br(),
                                      radioButtons("clinscenario", 
                                                   HTML("<strong>What is the normal morbidity?</strong> "), 
                                                   choices = list("High number of sick pigs routine causes" = "High", 
                                                                  "Low number of sick pigs routine causes" = "Low", 
                                                                  "Custom" = "Custom"),
                                                   # "Micro" = "micro"),
                                                   selected = "High"),
                                      actionButton(inputId = "Updatemorb",label = HTML("Update mortality<br> & morbidity"))
                                      
                               ),
                             column(width = 5, 
                                    br(),
                                    
                                    
                                     htmlOutput(outputId = "clintextoutput"),
                                     column(width=3,offset =0 , 
                                      
                                      tipify(textInput("txtinp_sickmin", "Lower", value = "1", width = NULL,
                                                       placeholder = NULL),placement="bottom",title = "Enter the lower value of<br> the sick pigs in routine production", trigger = "hover"),
                
                               ),
                               column(width=3,offset = 0 ,
                                      tipify(textInput("txtinp_sickmode", "Central", value = "2", width = NULL,
                                                       placeholder = NULL),placement="bottom",title = "Enter the most likely value of<br> the sick pigs in routine production", trigger = "hover"),
                               ),
                               column(width=3,offset = 0 ,
                                      tipify(textInput("txtinp_sickhigh", "High", value = "4.5", width = NULL,
                                                       placeholder = NULL),placement="bottom",title = "Enter the upper value of<br> the sick pigs in routine production", trigger = "hover")
                                      
                                     
                               
                                      
                                      
                                      ),
                               
                               br(),
                               br(),
                               br(),
                               tipify(checkboxInput("chkinp_usehighmort", "Use high normal mortality parameters", value = TRUE,
                                               ),placement="bottom",title = "Select to use data from<br> high normal mortality herds in simulation", trigger = "hover")
                               
                               
                               
                               
                               
                               ), 
                             
                            
                             
                           ####modal for sick pigs UI####  
                             column(width = 4,
                                    br(),
                             conditionalPanel(
                              condition = "input.inp_ASFstrain != '5'",
                             sliderInput(inputId = "sickidprob",label = HTML("What percent of pigs with mild clinical signs are identified for sampling by farm personnel?"),min = 0,max = 100,value = 85, step = 1,ticks = FALSE)
                             ),
                            
                             actionButton(inputId = "deflin_but",label = "Use default")
                             )
                           ),
                             
                          conditionalPanel( 
                            condition = "1>2",
                             hr( color="black" ,size=1.2),
                             
                             fluidRow(
                               column(width = 5, 
                                      
                                      
                              tipify(el = sliderInput(inputId = "Between_pen_people_percent",label = HTML("Enter the percent of between pen spread via distance independent pathways (fomites)"),min = 0,max = 100,value = 5,step = 1,ticks = FALSE),
                                     title = "This is the percent of between pen spread via distance independent pathways  such as people or fomites. A higher value results in a greater disease jumping vs clustering",placement = "bottom")      
                                      
                                      
                                      )) )
                             
                             ),
              ),
              
              column(width=3, 
                     
                     
                     checkboxInput("chbox_comp", label="Compare vs previous simulation ", value = FALSE, width = NULL)
              )),
              
              column(width = 10,offset = 1,
              htmlOutput(outputId = 'text12') 
              )
             )
           ),
           
          
####Premovement surveillance input side bar panel UI####

           box(solidHeader = T, collapsible = T, title = "Pre-movement surveillance inputs", status = "danger",width=12,collapsed = TRUE,
             
           
             fluidRow (
               fluidRow(
                 column(width = 1),
                 column(width=10,offset = 0,
                        conditionalPanel(condition = "1>2",         
               div(h5(HTML("<b>   Enter potential days post herd exposure when movement may occur</b>"))),
                ) ),
               ),
            
               fluidRow(
               column(width = 1),
               column(width=5, 
                      
                      
                      conditionalPanel(condition = "1>2",                  
             textInput("txtinp_startpredays",  h5("Start day"), value = "1", 
                       placeholder = NULL)),
                   
             sliderInput("txtinp_pmip_len",  h5("PMIP duration"), value = "0",min = 0,max = 20,step = 1,ticks = FALSE 
                       ),
           
             actionButton(width = "100%","runsurv", div(HTML("Run premovement<br/> surveillance")))
               ) ,      
            

             
             column(width=5,
                    conditionalPanel(condition = "1>2", 
                      textInput("txtinp_endpredays",  h5("End day"), value = "30", 
                                placeholder = NULL))
                 
                  ,
                    sliderInput("txtinp_pmip_eff",  h5("PMIP effectiveness"),min = 0,max = 100, value = "100", step = 1,ticks = FALSE
                              ),
                  
             
                    actionButton(width = "100%","Edit_protocol", div(HTML("Edit premovement<br/>surveillance  protocol"))),
                    br(),
                    br(),
                    downloadButton(width = "100%",outputId = "downsurvres",label = HTML("Download <br>simulation results"),
                                   style = "color: #fff; background-color: #27ae60; border-color: #fff;width: 100%")
                    
           
                    )
             )),
            fluidRow( 
              column(width = 1),  
           
            
  
            column(width = 5,offset = 0,
                                   br(),
            column(width = 5,
                   br(),

            )
            )
            
            )
           ),
         
####Detection curve inputs side bar panel UI####
conditionalPanel(condition = "2>1", 
           box(id = "detbox",solidHeader = T,status='info',collapsible = T, title = "Detection curve inputs",width=12,headerBorder = FALSE,collapsed = TRUE, 
               fluidRow (
                 fluidRow(
                   column(width = 1),
                   column(width=10,offset = 0,
                          div(h5(HTML("<b>   Enter potential movement days post exposure <br> to generate detection curve</b>"))),
                   ),
                 ),
                
                 fluidRow(
                   column(width = 1),
                   column(width=3, 
                          textInput("txtinp_startdetdays",  h5("Start day"), value = "1", 
                                    placeholder = NULL)
                       
                   ),
                   
                   
                   column(width=3,
                          
                          textInput("txtinp_enddetdays",  h5("End day"), value = "30", 
                                    placeholder = NULL)
                    
                   )
                   ,
                   
                   column(width=3,
                          
                          textInput("txtinp_enddetby",  h5("by"), value = "1", 
                                    placeholder = NULL)
                         
                   )
                   
                 ),
                 fluidRow( 
                   column(width = 1),  
                   
                   
                   
                   column(width = 5,offset = 0,
                          br(),
                          actionButton(width = "100%","rundet", div(HTML("Generate detection<br/> curve")))
                   ),
                   column(width = 5,
                          br(),
                          
                          downloadButton(outputId = "downdetres",label = HTML("Download <br>simulation results"),
                                         style = "color: #fff; background-color: #27ae60; border-color: #fff;width: 100%")
                          
                          
                   )
                 )
          
           )))

    ),
    
    
    
 ##### Output Tabs UI #####   
    column(width=7,
           
           tabsetPanel(id = "inTabset",
                       
             tabPanel("Transmission plots",
                      plotOutput(outputId ='graph2'),
                      selectInput("state", "Choose disease state to plot:",
                                                 list("Susceptible","Latent recovering","Latent dying","Infectious recovering","Infectious dying","Recovered","Dead","Clinical","Blood positive"),selected ="Dead" 
                                    ),
                      plotOutput(outputId ='graph1')
             ),
           
           tabPanel("Transmission tables",
                    br(),
                    fluidRow(
                      tags$head(
                        tags$style(type="text/css", "#dumm1 label{ margin-left: 20px;width: 80px;display: table-cell; text-align: center; vertical-align: middle; }
                        #dumm1 input{ margin-left: 10px;width: 50px;}
                       #dumm1 .form-group { display: table-row;}")
                      ),
                      
                      tags$head(
                        tags$style(type="text/css", "#dumm2 label{ margin-left: 10px;width: 50px;display: table-cell; text-align: center; vertical-align: middle; }
                        #dumm2 input{width: 40px}
                       #dumm2 .form-group { display: table-row;}")
                      ),
                      column(width=9,
                      column(width = 3,
                      tags$div(id = "dumm1",tipify(textInput("txtinp_startday1", "Start day", value = "1", 
                                                              placeholder = NULL),placement="bottom",title = "Enter the start day post infection<br>", trigger = "hover")),
                      ),
                      column(width = 3,
                      tags$div(id = "dumm1",tipify(textInput("txtinp_stop1", "Stop day", value = "20", 
                                                             placeholder = NULL),placement="bottom",title = "Enter the start day post infection<br>", trigger = "hover")),
                      ),
                      column(width = 3,
                      tags$div(id = "dumm2",tipify(textInput("txtinp_by1", "By", value = "2", 
                                                             placeholder = NULL),placement="bottom",title = "Enter the start day post infection<br>", trigger = "hover"))
                      )
                      )
                      ),
                    br(),
                    DT::dataTableOutput('Dtbl_1inf'),
                    downloadButton('downloadinftbl', 'Download infectious table'),
                    hr(),
                    br(),
                 
                    selectInput("stateftbl", "Choose disease state",
                                list("Susceptible","Latent recovering","Latent dying","Infectious recovering","Infectious dying","Recovered","Dead","Clinical","Blood positive"),selected ="Dead" 
                    ),
                    DT::dataTableOutput('Dtbl_anystate'),
                    downloadButton('downloadstatetable', 'Download state table'),
                    
                    
           ),
           tabPanel(div(HTML("Morbidity & mortality")),
                  
                    br(),
                    fluidRow(
                      tags$head(
                        tags$style(type="text/css", "#dumm1 label{ margin-left: 20px;width: 80px;display: table-cell; text-align: center; vertical-align: middle; }
                        #dumm1 input{ margin-left: 10px;width: 50px;}
                       #dumm1 .form-group { display: table-row;}")
                      ),
                      
                      tags$head(
                        tags$style(type="text/css", "#dumm2 label{ margin-left: 10px;width: 50px;display: table-cell; text-align: center; vertical-align: middle; }
                        #dumm2 input{width: 40px}
                       #dumm2 .form-group { display: table-row;}")
                      ),
                      column(width=9,
                             column(width = 3,
                                    tags$div(id = "dumm1",tipify(textInput("txtinp_startday1mm", "Start day", value = "1", 
                                                                           placeholder = NULL),placement="bottom",title = "Enter the start day post infection<br>", trigger = "hover")),
                             ),
                             column(width = 3,
                                    tags$div(id = "dumm1",tipify(textInput("txtinp_stop1mm", "Stop day", value = "20", 
                                                                           placeholder = NULL),placement="bottom",title = "Enter the start day post infection<br>", trigger = "hover")),
                             ),
                             column(width = 3,
                                    tags$div(id = "dumm2",tipify(textInput("txtinp_by1mm", "By", value = "2", 
                                                                           placeholder = NULL),placement="bottom",title = "Enter the start day post infection<br>", trigger = "hover"))
                             )
                      )
                    ),
                    br(),
                    
                  
                    
                      
                    DT::dataTableOutput('Dtbl_2clin'),
                  
                    downloadButton('downloadclin', 'Download clinical table'),
                    hr(),
                    br(),
                    DT::dataTableOutput('Dtbl_mort'),
                    downloadButton('downloaddailmort', 'Download mortality table'),
           ),
        
           tabPanel("Premove surveillance protocol",
                    fluidRow(
                      
                      column(width=4, 
                      tipify(checkboxInput(inputId = "chbox_editport1", label=(div(HTML("Edit surveillance protocol 1 <br> Only individual pig samples (e.g., blood)"),style="color:black;font-weight: bold;")), value = FALSE, width = NULL),placement="bottom",title = "select to edit protocol 1", trigger = "hover")),
                      column(width=8,tipify(checkboxInput(inputId = "chbox_editport2", label=(div(HTML("Edit surveillance protocol 2 <br> Only individual pig samples (e.g., blood)"),style="color:black;font-weight: bold;")), value = FALSE, width = NULL),placement="bottom",title = "select to edit protocol 2", trigger = "hover")
                      ),
                      
                      
                     ),
            
                    
                    div(HTML("For this educational application, the diagnostic sensitivity of the RRT-PCR test for individual pig testing was set to 95%"),style="color:black;font-weight: bold;"),
                  br(),
                    
                    fluidRow(
                      conditionalPanel(condition = "input.chbox_editport1",
                      tags$head(
                        tags$style(type="text/css", "#dumm3 label{ margin-left: 20px;width: 200px;display: table-cell; text-align: center; vertical-align: middle; }
                        #dumm3 input{ margin-left: 10px;width: 100px;}
                       #dumm3 .form-group { display: table-row;}")
                      ),
                      
                      tags$head(
                        tags$style(type="text/css", "#dumm2 label{ margin-left: 10px;width: 50px;display: table-cell; text-align: center; vertical-align: middle; }
                        #dumm2 input{width: 40px}
                       #dumm2 .form-group { display: table-row;}")
                      ),
                      
                     
                      div(HTML("<b>Surveillance protocol 1 inputs:</b>"),style = "color:darkblue"),
                      column(width=4,
                             br(),
                             tags$div(id = "inline",tipify(numericInput("ninp_prot1samplesize", div(HTML("Number of individual <br> pig samples  "),style="color:black"), value = 20,max=100,min=1),placement="bottom",title = "Enter the sample size for individual pig testing<br>", trigger = "hover")),
                      
                             br(),
                             
                             tipify(sliderInput(inputId = "slinp_prot1_morttrig",min = 1,max=20,label = div(HTML("Enter the daily mortality trigger per 1000 pigs"),style="color:black"),value = 5,step = 1,ticks = FALSE),
                                    placement="right",title = "Enter the daily mortality trigger threshold", trigger = "hover")
                             
                             ),
                      column(width=4,
                             br(),
                             tipify(
                               selectInput(
                                 inputId = "selinp_indprot1days", 
                                 label = div(HTML("Days prior movement <br> for individual pig samples"),style="color:black"),
                                 selected = 1,
                                 choices = c("None",c(1:30)),multiple = TRUE
                                 
                               )
                               
                           ,placement="right",title = "Enter the days before movement when samples are collected<br>", trigger = "hover"),
                           
                           
                           
                           
                      ),
                      column(width=4,
                             br(),
                             tipify(selectInput("selinp_prot1_priority", (div(HTML("Choose priority scheme for individual pig samples"),style = "color:black")),
                                                list("Dead, sick, live"=1,"Sick, dead, live"=2,"Dead,live"=3,"Dead only"=4,"Dead,sick"=5,"Sick,apparently healthy"=6,"random live"=7),selected =2
                             ), placement="right",title = "Choose sample type priority for individual pig samples", trigger = "hover") 
                             
                             
                      ),
                      hr()
                    )
                    
                
                    ),
                    
                    conditionalPanel(condition = "input.chbox_editport2",
                    fluidRow(
                      tags$head(
                        tags$style(HTML("hr {border-top: 1px solid #000000;}"))
                      ),
                      
                      div(HTML("<b>Surveillance protocol 2 inputs:</b>"),style = "color:darkblue"),
                 
                   
                      column(width=4,
                             br(),
                             tags$div(id = "inline",tipify(numericInput("ninp_prot2samplesize", div(HTML("Number of individual<br> pig samples  "),style="color:black"), value = 20,max=100,min=1),placement="bottom",title = "Enter the sample size for individual pig testing<br>", trigger = "hover")),
                             
                             br(),
                             
                             tipify(sliderInput(inputId = "slinp_prot2_morttrig",min = 1,max=20,label = div(HTML("Enter the daily mortality trigger per 1000 pigs"),style="color:black"),value = 5,step = 1,ticks = FALSE),
                                    placement="right",title = "Enter the daily mortality trigger threshold", trigger = "hover"),
                             conditionalPanel(condition = "1>2", 
                             tipify(sliderInput(inputId = "slinp_percentpospensprot2",min = 0,max=100,label = div(HTML("Enter the percent of pens sampled via ropes"),style="color:black"),value = 0,step = 1,ticks = FALSE),
                                    placement="right",title = "Enter the percent of pens sampled via ropes", trigger = "hover"))
                             
                      ),
                      column(width=4,
                             br(),
                             tipify(
                               selectInput(
                                 inputId = "selinp_indprot2days", 
                                 label = div(HTML("Days prior movement <br> for individual pig samples"),style="color:black"),
                                 selected = 1,
                                 choices = c("None",c(1:30)),multiple = TRUE
                                 
                               )
                               
                               ,placement="right",title = "Enter the days before movement when samples are collected<br>", trigger = "hover"),
                             br(),
                             conditionalPanel(condition = "1>2", 
                             tipify(
                               
                               selectInput(
                                 inputId = "selinp_oralprot2days", 
                                 label = div(HTML("Days prior movement <br> for oral fluid samples"),style="color:black"),
                                 selected = "None",
                                 choices = c("None",c(1:30)),multiple = TRUE
                                 
                               )
                               
                               ,placement="right",title = "Enter the days before movement when samples are collected<br>", trigger = "hover")),
                            
                             br(),
                             conditionalPanel(condition = "1>2", 
                              actionButton(width = "100%","view_oral", div(HTML("Edit oral fluids <br/> sensitivity"))),
                             bsModal(id = "survmodal", title = "Oral fluids sensitvity",trigger =  "view_oral", size = "large",DTOutput(outputId = "oral_sens_frame2"),
                                     fluidRow(column(width = 3,
                                                     br(),
                                                     
                                                     downloadButton('downloadoralsens', 'Download'),
                                                     br()    
                                                     
                                                     
                                     ),
                                     
                                     column(width=8,
                                            
                                            
                                            
                                            
                                            htmlOutput('error_msg'),
                                            fileInput(inputId ='FileInp_loadoralsens',label="",accept = ".csv",width = NULL,buttonLabel = "Load file" )
                                     )
                                   
                                     ), plotOutput('ropesenplot')
                                     
                                     
                             )
                      
                             
                             
                             )),
                      column(width=4,
                             br(),
                             tipify(selectInput("selinp_prot2_priority", (div(HTML("Choose priority scheme for individual pig samples"),style = "color:black")),
                                                list("Dead, sick, live"=1,"Sick, dead, live"=2,"Dead,live"=3,"Dead only"=4,"Dead,sick"=5,"Sick,apparently healthy"=6,"random live"=7),selected =2 
                             ), placement="right",title = "Choose sample type priority for individual pig samples", trigger = "hover"),
                             
                             br(),
                             conditionalPanel(condition = "1>2", 
                             tipify(selectInput("selinp_oral2_priority", (div(HTML("Choose priority scheme for oral fluids"),style = "color:black")),
                                                list("Random pen sampling"=1,"prioritizing pens with sick pigs"=2,"Prioritizing pens with dead pigs"=3,"Pens with dead and then sick pigs"=4),selected =4 
                             ), placement="right",title = "Choose sample type priority for individual pig samples", trigger = "hover") )
                             
                             
                             
                      )
                    )),
                    
                    shinyjs::hidden(box("Surveillance protocol",
                                        fluidRow(
                                          tags$head(
                                            tags$style(type="text/css", "#dumm3 label{ margin-left: 20px;width: 200px;display: table-cell; text-align: center; vertical-align: middle; }
                        #dumm3 input{ margin-left: 10px;width: 100px;}
                       #dumm3 .form-group { display: table-row;}")
                                          ),
                                          
                                          tags$head(
                                            tags$style(type="text/css", "#dumm2 label{ margin-left: 10px;width: 50px;display: table-cell; text-align: center; vertical-align: middle; }
                        #dumm2 input{width: 40px}
                       #dumm2 .form-group { display: table-row;}")
                                          ),
                                          column(width=11,
                                                 br(),
                                                 column(width = 6,
                                                        tags$div(id = "dumm3",tipify(textInput("txtinp_selected_protocol", "Selected protocol row", value = "1", 
                                                                                               placeholder = NULL),placement="bottom",title = "Enter the selected protocols to simulate<br> e.g., 1,2", trigger = "hover")),
                                                 ),
                                                 
                                                 column(width = 5,
                                                        tipify(checkboxInput(inputId = "chbox_selallprot", label=div(HTML("Simulate all rows")), value = TRUE, width = NULL),placement="bottom",title = "select to simulate all protocols", trigger = "hover")
                                                 )
                                          )
                                        ),
                                        
                                        fluidRow(column(align = "center",width = 12,
                                                        htmlOutput('error_msgprecheck'),                 
                                                        DTOutput(outputId = "pre_mov_frame",width = '100%'))),
                                        fluidRow(
                                          
                                          column(width=11,
                                                 htmlOutput('error_msgpre'),
                                                 fileInput(inputId ='FileInpLoadpremov',label="",accept = ".csv",width = NULL,buttonLabel = "Load file" )
                                          )),
                                        plotOutput('graph9',click = clickOpts("plot_click3"),dblclick = dblclickOpts("plot_dblclick3")),
                                        plotOutput('graph10',click = clickOpts("plot_click4"),dblclick = dblclickOpts("plot_dblclick4"))
                    )),
                    
                    
           ),
           
           
           tabPanel("Surveillance results",
                    fluidRow(column(align = "center",width = 12,
                                    DTOutput(outputId = "survresout",width = '100%')))
           ),
           
           
           
           tabPanel("Detection curve",
                    plotOutput('graph7',click = clickOpts("plot_click3"),dblclick = dblclickOpts("plot_dblclick3")),
                    fluidRow(column(align = "center",width = 12,
                                    DTOutput(outputId = "Deccurvoutput"))),
                    
                    plotOutput('graph8',click = clickOpts("plot_click4"),dblclick = dblclickOpts("plot_dblclick4"))
           ),
           
           tabPanel("Manual",
                    uiOutput("pdf_viewer")    
                   #htmlOutput("manual2")
           )
           
                    
                    
           )
           
           ) 
           
    )
  
 
  

  
)

####Server top####
server<-function(input,output,session)
{
  
  
  initplotreactval<-reactiveVal(0)
  initcontactselectinp<-reactiveVal(0)
  
  
  # Simulation parameters
  sim_parms <- list()
  
  # simulation parameters initialized below
  sim_parms$num_iterations = 500 # Number of transmission model iterations to run
  sim_parms$n_runs = 1 #do not change
  sim_parms$time_step = 0.01 # Transmission model time step in terms of days
  sim_parms$ndays = 35 # Number of days to run transmission model
  sim_parms$readdt <- 1 # Time step in days to store transmission model output
  initial_seed <- 40841751 # Sets seed for random number generation
  sim_parms$init_inf=1 # Initial number of infected pigs
  sim_parms$init_pen=NULL # Pen containing the initially infected pig(s). Randomly selected in transmission function 'get_single_sim'
  sim_parms$ignoredieouts <- TRUE # T/F to keep or exclude transmission model runs where infection does not take off in the barn
  # 1 = Guinat; 2 = royal vet med; 3 = hohle; 4 = ausvet; 5 = Perez and Klinkenberg; 6 = Nielson; 7 = adjacent + distance independent between pen spread
  # Note: delay in collecting dead pigs only implemented for transmission type 7
  sim_parms$trans_type<-7 # Allows for different formulations of between pen transmission
  
  #### Strain Scenario Parameters ####
  # These parameters are changed in 'sim_sate_parms_adjust' function defined in strain_beta_scenario_ftns.R
  sim_parms$is_georgia_strain=FALSE # modifies how the time to clinical signs appearing is simulated when True
  sim_parms$mort_lat_pars = c(13.299,0.3384482) #c(shape, scale) gamma for latent period of pigs that die following exposure
  sim_parms$mort_inf_pars = c(9.632,0.862) #c(shape, scale) gamma for infectious period of pigs that die following exposure
  sim_parms$recov_lat = c(13.299,0.3384482) #shape and scale parameter for latent period of pigs that recover following exposure
  sim_parms$recov_inf = c(55.42012,0.7950162) #c(shape, scale) gamma for infectious period of pigs that recover
  sim_parms$p_mort=0.1 # Probability that a pig dies due to ASF following exposure
  sim_parms$bpospar=c(-0.82, 0.82) #mean and standard deviation of norm distribution for when a pig is blood positive relative to the start of the infectious period
  sim_parms$clinpar=c(3.418,3.2)#gamma shape and scale duration of mild clinical signs
  sim_parms$sero_pars = c(26.257,0.214) #gamma shape and scale parameter for time to mild clinical signs appearing following exposure
  sim_parms$georgia_clinmean = c(4.985, 1.2463) #mean and standard deviation for normal distribution for mean time to clincal signs for Georgia strain
  sim_parms$georgia_clinmean_bounds = c(3, 7.5) #min/max bounds for normal distributed mean time to clinical signs Georgia strain
  sim_parms$p_sero=1 #Fraction of birds that develop mild clinical signs following exposure
  
  # ASF strain scenario list
  ttstrainnamelist=c("Moderately virulent","Highly virulent olsen","Highly virulent Georgia hu","Highly virulent Olensen long", "Modvir for severeclin")
  # Transmission scenario list
  ttcontactratelist=c("Slow (highly virulent strain estimate)","Fast (highly virulent strain estimate)","Moderately virulent strain estimate", "User defined contact rate")
  strain_namereac<-reactiveVal(ttstrainnamelist[1])
  contact_namereac<-reactiveVal(ttcontactratelist[1])
  
  state_parm_set_list=list("Moderately_virulent"=1,"olsen"=2,"hu"=3,"olsensens"=4, "Modvir_severeclin"=5)
  sim_parms$tsclin_par=c(26.257,0.214) #shape and scale of gamma distribution for time until onset of severe clinical signs
  sim_parms$lsclin_par=c(3.418,3.2) #shape and scale of gamma distribution of duration of severe clinical signs
  sim_parms$mort_par_bounds=c(-1,-1,-1, -1) #min/max bounds for latent and infectious periods for pigs that die
  sim_parms$recov_par_bounds=c(-1,-1,-1,-1) #min/max bounds for latent and infectious periods for pigs that recover
  sim_parms$bpos_par_bounds=c(-4,4) #min/max bounds for time of blood positive relative to onset of infectiousness
  sim_parms$mortdelaypar=c(0,0) #min/max in days of uniform distribution for how long a dead pig remains in the pen before removal
  need_run_normmorb_atleast<-reactiveVal(FALSE)
  ###########################################################################
  
  #### Contact Rate Scenario Parameters ####
  # These parameters are changed in 'sim_beta_pars' function defined in strain_beta_scenario_ftns.R
  sim_parms$betapars=c(1.14,2.01,3.55)# beta-PERT parameters (min, mode, max) for within-pen contact rate
  sim_parms$betabetween=c(0.16,0.46,1.06)# beta-PERT parameters (min, mode, max) for between pen contact rate
  sim_parms$betaroom=c(0.00001,0.01,0.1) # beta-PERT parameters (min, mode, max) for between room contact rate
  
  #Other contact rate parameters
  beta_set_list <- list("georgia_slow"=1, "georgia_fast"=2, "Moderately_virulent"=3)
  sim_parms$between_penpeoplefrac <- 0.05 # percent of between pen transmission that is distance independent (used when sim_parms$trans_type=7)
  sim_parms$morttrans_mult <- 1 #Scaling factor for increasing/decreasing disease transmission related to mortality left in pens
  
  ###########################################################################
  
  #### Premises Structure Parameters ####
  sim_parms$flock_size_pars=c(1200) #Total number of pigs on premises
  sim_parms$npens<-30 #Total number of pens
  sim_parms$nrooms=1 # Number of rooms in the barn
  sim_parms$prem_config<-2 # Options for pen layout within barn: 1-side alley two rooms; 2-central alleyway between two rows of pens
  
  ###########################################################################
  
  #### Surveillance Parameters ####
  sim_parms$move_day_bounds<-c(1,30) #Bounds for day post herd exposure that movement may occur (e.g. c(3,30) is exposure occurs between 3 and 30 days prior to movement day)
  sim_parms$PCRsens=0.95 #blood/tissue PCR test sensitivity
  sim_parms$pcr_targetswabs=5 #number of swabs per pool; divides up total number of samples given by premovsurvlist$Premov.sample.size
  sim_parms$dieout_mortthresh <- 5 #if max mortality less than or equal parameter the iteration is counted as a dieout
  sim_parms$normsickparms= c(0.01,0.02,0.045)# PERT distribution (min, mode, max) for number of healthy pigs with clinical signs. Scenarios include:
  # c(0.0025,0.005,0.04) for low proportion of pigs with mild clinical signs; c(0.01,0.02,0.045) for high proportion of pigs with mild clinical signs;
  # and c(0.01, 0.01, 0.01) for the proportion of pigs with severe clinical signs during routine production
  sim_parms$propsickobs = 0.85 #proportion of sick pigs that are successfully identified (e.g. if zero no pigs with clinical signs ever observed)
  sim_parms$boolretnumsamp = TRUE #input for perform_all_iter_premovsurv function; if TRUE returns number of sick, dead, and live pigs tested by PCR
  sim_parms$ret_infected_vs_infectious<-FALSE #if true returns number of infected pigs when undetected (infectious + latent), if false returns only number of infectious undetected pigs
  sim_parms$samp_ind_viablood_not_oral<-TRUE #if true does individual pig blood sampling, if false does oral sampling
  sim_parms$use_high_mortdata<-TRUE#if true uses high weekly mortality between 95th and 99.5th percentile else uses all herds in normal mortality simulation
  sim_parms$sampsize=30 #total number of swabs per premises; default value in PCR test functions (not actively used)
  
  infodatframe<-read.csv(file = paste("premoveinfodat.csv",sep=""),header = TRUE,stringsAsFactors = FALSE)
  
  sim_parms$nsurvdays = c(1,30) #min, max for days to perform surveillance; used in perform_alliter_test_timeseq
  sim_parms<<-sim_parms
  sim2parms<-reactiveVal(sim_parms)
  sim3parms<-reactiveVal(sim_parms)
  
  
  #####################################################################################################################
  
  # reads in functions for simulating virus transmission and formatting the simulation output
  source("run_trans_funcs.R")
  # reads in functions for simulating routine mortality and morbidity as well as performing PCR testing of dead tissue samples and/or blood samples from sick and apparently healthy pigs
  source('Blood_detect_exp.R')
  # reads in functions for simulating oral fluids sampling and testing
  source('Oral_fluids_detect_exp.R')
  # file contains helper functions for linearizing a 2-d array, a pert distribution simulator, and a wrapper for simulating from a hypergeometric distribution
  source("general_utility_function.R")
  # gets premises structure in array form (given by pen_distmult_linvec and room_distmult_linvec)
  source("Rooms_pens_distance_inp_nielson.R")
  # functions to simulate the active surveillance protocols
  # file calls Blood_detect_exp.R and Oral_fluids_detect_exp.R
  source('run_surv_scenario.R')
  # functions to make plots and tables of simulation results
  # file calls trans_post_process_func.R
  source('Transmission model visualization.R')
  
  
  state_parm_set_list<<-state_parm_set_list
  beta_set_list<<-beta_set_list
  #print(beta_set_list)
  # Reads in functions specifying parameters for strain and contact rate scenarios
  source('strain_beta_scenario_ftns.R')
  # Reads in functions for simulating and processing pre-movement active surveillance
  source('pre_movsurv_ftns.R')
  # Functions for making tables and plots of simulation results
  source('Plots for shiny.R')
  # Functions for making tables from transmission model output
  source('clin_mort_tableftns.R')
  
  
  pstate_list_names<<-c("Susceptible","Latent recovering","Latent dying","Infectious recovering","Infectious dying","Recovered","Dead","Clinical","Blood positive","Severe clinical")
  
  inp_strain_list_react<-reactiveVal(1)
  inp_beta_choice_react<-reactiveVal(3)
  #this is temperory to only update to the above value when they click inputdo
  temp_beta_choice_react<-reactiveVal(3)
  custom_withinpenB_reacvec<-reactiveVal(c(1.00,1.6,2.7))
  sim_inpcustom_withinpenB_reacvec<-reactiveVal(c(1.00,1.6,2.7))
  
  custom_betweenpenB_reacvec<-reactiveVal(c(0.1,0.3,0.5))
  sim_inpcustom_betweenpenB_reacvec<-reactiveVal(c(0.1,0.3,0.5))
  
  inp_graph1_disstate<-reactiveVal(7)
  #selected disease state for table
  inp_tbl1_disstate<-reactiveVal(7)
  #Number of days in data table
  tbl_outdays_reac<-reactiveVal(seq(1,20,by=2))
  tbl_outdays_reacmm<-reactiveVal(seq(1,20,by=2))
  inp_pensizereac<-reactiveVal(40)
  inp_total_herdsize<-reactiveVal(1200)
  barn_pen_arraylist<-reactiveVal()
  temp_barn_pen_arraylist<-reactiveVal()
  premovsurframe_init<-read.csv(file = "premovsurvlist.csv",header = TRUE,stringsAsFactors = FALSE)
  premovsurframe_init$Premove.survdays<-as.character(premovsurframe_init$Premove.survdays)
  premovsurframe_init$Premov.sample.size<-as.character(premovsurframe_init$Premov.sample.size)
  premovsurframe_init$Oral.surv.days<-as.character(premovsurframe_init$Oral.surv.days)
  premovsurframe_reac<-reactiveVal(premovsurframe_init)
  tableedit_premovsurframe_reac<-reactiveVal(premovsurframe_init)
  oral_fluid_sensframe_init<-read.csv(file = "oral_fluids_sens.csv",header = TRUE,stringsAsFactors = FALSE)
  oral_fluid_sensframe_reac<-reactiveVal(oral_fluid_sensframe_init)
  
  
  ####premovesurvlist####
  sim_premovsurframe_reac<-reactiveVal(premovsurframe_init)
  # functions to format active surveillance protocol details given in premovsurvlist.csv
  source('pre_movsurv_list_creator.R'  )
  pre_mov_sur_list<-reactiveVal()
  pre_mov_sur_list(get_premov_survlist(premovsurframe_init))
  premov_start_reac<-reactiveVal(1)
  premov_end_reac<-reactiveVal(30)
  premov_dur_reac<-reactiveVal(0)
  premov_eff_reactive<-reactiveVal(1)
  
  sel_premov_survreac<-reactiveVal(1)
  pensizearray_reac<-reactiveVal()
  input_protocol_vec_reac<-reactiveVal(1)

  
  
  ####clinical signs reactives####
  inp_normclin_reacvec<-reactiveVal(c(0.01,0.02,0.045))
  inp_probclinident<-reactive(0.85)
  testreactvevalues<-reactiveValues()
  clin_textout_reactive=reactiveVal(value = "<b>WCustom input: what is the percent of <i> sick pigs</i> in routine production 
(i.e., among ASF suceptible pigs)</b>")
  
  
  
  ###############################################
  # Random movement day premovement surveillance#
  ###############################################
 
  updateTabsetPanel(session, "inTabset",
                    selected ="Transmission plots")
  trans_shin_list<-readRDS(file = "initial_translist")
  
  prevtrans_listreactive<-reactiveVal(trans_shin_list)
  sim_outputreactive<-reactiveVal(trans_shin_list)
  is_detcurv_avil<-reactiveVal(FALSE)
  should_trans_runbeforesurv<-reactiveVal(FALSE)
  onlyupdatenormortmob<-reactiveVal(FALSE)
  first_t_selecinput<-reactiveVal(TRUE)
  first_t_selecinput2<-reactiveVal(TRUE)
  detc_outreac_frame<-reactiveVal()
  showprotocol1<-reactiveVal(FALSE)
  showproocol2<-reactiveVal(FALSE)
  
  
  observeEvent(input$state,{
    inp_graph1_disstate(which(pstate_list_names==input$state))
    
    
  })
  
  
  observeEvent(input$clinscenario,{
    
    if(input$clinscenario=="High"){
      inp_normclin_reacvec(c(.01, 0.02, 0.045))
      updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickmin',value =(inp_normclin_reacvec()[1]*100)  )
      updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickmode',value =(inp_normclin_reacvec()[2]*100)  )
      updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickhigh',value =(inp_normclin_reacvec()[3]*100)  )
      shinyjs::disable("txtinp_sickmin")
      shinyjs::disable("txtinp_sickmode")
      shinyjs::disable("txtinp_sickhigh")
    }else if (input$clinscenario=="Low"){
      inp_normclin_reacvec(c(0.0025, 0.005, 0.04))
      updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickmin',value =(inp_normclin_reacvec()[1]*100)  )
      updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickmode',value =(inp_normclin_reacvec()[2]*100)  )
      updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickhigh',value =(inp_normclin_reacvec()[3]*100)  )
      shinyjs::disable("txtinp_sickmin")
      shinyjs::disable("txtinp_sickmode")
      shinyjs::disable("txtinp_sickhigh")
    }else{
      shinyjs::enable("txtinp_sickmin")
      shinyjs::enable("txtinp_sickmode")
      shinyjs::enable("txtinp_sickhigh")
      
    }
    
    
  })
  
  
  observe({
    input$txtinp_betamin
    input$txtinp_betamode
    input$txtinp_betahigh
    
    custom_withinpenB_reacvec(c(as.numeric(input$txtinp_betamin),as.numeric(input$txtinp_betamode),as.numeric(input$txtinp_betahigh)))

    initcontactselectinp(1)
  })
 
  
  
  observe({
    input$txtinp_bpbetamin
    input$txtinp_bpbetamode
    input$txtinp_bpbetahigh
    custom_betweenpenB_reacvec(c(as.numeric(input$txtinp_bpbetamin),as.numeric(input$txtinp_bpbetamode),as.numeric(input$txtinp_bpbetahigh)))
    
    initcontactselectinp(1)
  })
  
  #clinical
  toListen <- reactive({
    list(input$txtinp_sickmin,
         input$txtinp_sickmode,
         input$txtinp_sickhigh,
         input$chkinp_usehighmort
         
         
         
         )
  })
  observeEvent(toListen() ,{
   
    if(as.numeric(initplotreactval())!=0){
     
     need_run_normmorb_atleast(TRUE)
    }
    isolate({inp_normclin_reacvec(c(as.numeric(input$txtinp_sickmin)/100,as.numeric(input$txtinp_sickmode)/100,as.numeric(input$txtinp_sickhigh)/100))
      
      
      })
    
    if(!first_t_selecinput()){
      shinyjs::disable("runsurv")
      shinyjs::disable("rundet")
      outtextreactive("Please process morbidity or run transmission before evaluating surveillance")
      
    }
    first_t_selecinput(FALSE) 

    
  })
  
  observeEvent(toListen2() ,{
   
    if(!first_t_selecinput2()){
      shinyjs::disable("runsurv")
      shinyjs::disable("rundet")
      outtextreactive("Please run transmission before evaluating surveillance")
      should_trans_runbeforesurv(TRUE)
      
    }
    first_t_selecinput2(FALSE) 
  })
  
  
  toListen2 <- reactive({
    list(
         input$Inp_barnsize,
         
        
         input$Between_pen_people_percent,
         input$inp_contact_scenario,
       
         input$inp_ASFstrain,
         input$Inp_pensize,
      
         input$txtinp_numiter
         
         )
  })
  
  
  
 

  observe({
    input$Inp_pensize
    input$Inp_barnsize
    
    
    inp_total_herdsize(as.numeric(input$Inp_barnsize))
    inp_pensizereac(as.numeric(input$Inp_pensize))
    
    tempsim_parms<-sim2parms()
    tempsim_parms$flock_size_pars=c(inp_total_herdsize()) 
    tempsim_parms$npens<-inp_total_herdsize()/inp_pensizereac()
  
    temp_barn_pen_arraylist(get_barn_pen_arraylist(sim_parms =tempsim_parms  ))
   
    
    }
  )
  
  observeEvent(input$deflin_but, {
    
    if(  input$inp_ASFstrain=='5'){
      inp_normclin_reacvec(c(0.0099, 0.01, 0.0101))
    }else{
      inp_normclin_reacvec(c(0.01,0.02,0.045))
    }
   
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickmin',value =(inp_normclin_reacvec()[1]*100)  )
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickmode',value =(inp_normclin_reacvec()[2]*100)  )
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickhigh',value =(inp_normclin_reacvec()[3]*100)  )
    shinyjs::disable("txtinp_sickmin")
    shinyjs::disable("txtinp_sickmode")
    shinyjs::disable("txtinp_sickhigh")
    updateRadioButtons(session = getDefaultReactiveDomain(),inputId ='clinscenario',selected = "High")
    updateSliderInput(session = getDefaultReactiveDomain(),inputId ='sickidprob',value ="85"  )
  })
  
  observeEvent(input$inp_contact_scenario,{
    temp_beta_choice_react(as.numeric(input$inp_contact_scenario))
    
    shinyjs::disable("txtinp_betamin")
    shinyjs::disable("txtinp_betamode")
    shinyjs::disable("txtinp_betahigh")
    shinyjs::disable("txtinp_bpbetamin")
    shinyjs::disable("txtinp_bpbetamode")
    shinyjs::disable("txtinp_bpbetahigh")
    
    
    if(  temp_beta_choice_react()==3){
      custom_withinpenB_reacvec(c(1.00,1.6,2.7))
  
    }else if(temp_beta_choice_react()==2){
      custom_withinpenB_reacvec(c(1,2.6,5.6)) 
    }else if(temp_beta_choice_react()==1){
      custom_withinpenB_reacvec(c(0.3,0.6,1.0)) 
    }else{
      shinyjs::enable("txtinp_betamin")
      shinyjs::enable("txtinp_betamode")
      shinyjs::enable("txtinp_betahigh")
      shinyjs::enable("txtinp_bpbetamin")
      shinyjs::enable("txtinp_bpbetamode")
      shinyjs::enable("txtinp_bpbetahigh")
      custom_withinpenB_reacvec(c(as.numeric(input$txtinp_betamin),as.numeric(input$txtinp_betamode),as.numeric(input$txtinp_betahigh)))
    }
    
    if(  temp_beta_choice_react()==3){
      custom_betweenpenB_reacvec(c(0.1,0.3,0.5)) 
    }else if(temp_beta_choice_react()==2){
      custom_betweenpenB_reacvec(c(0.31,0.99,1.98)) 
    }else if(temp_beta_choice_react()==1){
      custom_betweenpenB_reacvec(c(0.1,0.3,0.5)) 
    }else{
      custom_betweenpenB_reacvec(c(as.numeric(input$txtinp_bpbetamin),as.numeric(input$txtinp_bpbetamode),as.numeric(input$txtinp_bpbetahigh)))
    }
    
    isolate({updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_betamin',value =custom_withinpenB_reacvec()[1]  )
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_betamode',value =custom_withinpenB_reacvec()[2]  )
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_betahigh',value =custom_withinpenB_reacvec()[3]  )})
    
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_bpbetamin',value =custom_betweenpenB_reacvec()[1]  )
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_bpbetamode',value =custom_betweenpenB_reacvec()[2]  )
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_bpbetahigh',value =custom_betweenpenB_reacvec()[3]  )
    
    print("the contact rate scenario from select is")
    print(temp_beta_choice_react())
    
   })
  
  observeEvent(input$inp_ASFstrain,{

    if(  input$inp_ASFstrain=='5'){
      inp_normclin_reacvec(c(0.0099, 0.01, 0.0101))
    }else{
      inp_normclin_reacvec(c(0.01,0.02,0.045))
    }
    
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickmin',value =(inp_normclin_reacvec()[1]*100)  )
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickmode',value =(inp_normclin_reacvec()[2]*100)  )
    updateTextInput(session = getDefaultReactiveDomain(),inputId ='txtinp_sickhigh',value =(inp_normclin_reacvec()[3]*100)  )
    shinyjs::disable("txtinp_sickmin")
    shinyjs::disable("txtinp_sickmode")
    shinyjs::disable("txtinp_sickhigh")
    updateRadioButtons(session = getDefaultReactiveDomain(),inputId ='clinscenario',selected = "High")
    
    
   if( input$inp_ASFstrain =='5'){
    clin_textout_reactive( HTML("<b>Input percent of <i>severely sick pigs</i> in routine production 
(i.e., among ASF suceptible pigs)</b>"))
     
    
   }else{
     clin_textout_reactive(HTML( "<b>Custom input: what is the percent of <i> sick pigs</i> in routine production 
(i.e., among ASF suceptible pigs)</b>"))
     
     
     
   }
  })

  
  if_table_edited<-reactiveVal(FALSE)
  
  toListenpremove <- reactive({
    
    list(input$ninp_prot1samplesize,
         
         input$selinp_indprot1days,
         
         input$selinp_prot1_priority,
         input$slinp_prot1_morttrig,
         
         
         input$ninp_prot2samplesize,
         
         input$selinp_indprot2days,
         
         input$selinp_prot2_priority,
         input$slinp_prot2_morttrig,
         input$selinp_oralprot2days,
         input$selinp_oral2_priority,
         input$slinp_percentpospensprot2
         
         
    )
    
  })
  
  
  observeEvent(toListenpremove(),{
    temp_data=tableedit_premovsurframe_reac()
    temp_data[1,4]=input$ninp_prot1samplesize
  
    if(sum(grepl(pattern = "None",input$selinp_indprot1days))>0){
      temp_data[1,3]=NA
    }else{
      temp_data[1,3]=paste(input$selinp_indprot1days,collapse = ",")
      if(length(input$selinp_indprot1days)>1){
       
        temp_data[1,4]=paste(rep(input$ninp_prot1samplesize,length(input$selinp_indprot1days)),collapse = ",")
       
      }
    }
   temp_data[1,9]=as.numeric(input$selinp_prot1_priority) 
    
   temp_data[1,12]=as.numeric(input$slinp_prot1_morttrig) 
   
   
   ###prot 2
   temp_data[2,4]=input$ninp_prot2samplesize
   
   if(sum(grepl(pattern = "None",input$selinp_indprot2days))>0){
     temp_data[2,3]=NA
   }else{
     temp_data[2,3]=paste(input$selinp_indprot2days,collapse = ",")
     if(length(input$selinp_indprot2days)>1){
   
       temp_data[2,4]=paste(rep(input$ninp_prot2samplesize,length(input$selinp_indprot2days)),collapse = ",")
       
     }
   }
   
   
   temp_data[2,9]=as.numeric(input$selinp_prot2_priority) 
   
   temp_data[2,12]=as.numeric(input$slinp_prot2_morttrig) 
   
   #oral stuff
   temp_data[2,10]=as.numeric(input$selinp_oral2_priority) 
   temp_data[2,7]=as.numeric(input$slinp_percentpospensprot2)
   
   if(sum(grepl(pattern = "None",input$selinp_oralprot2days))>0){
     temp_data[2,6]=NA
   }else{
     temp_data[2,6]=paste(input$selinp_oralprot2days,collapse = ",")
    
   }
   
    tableedit_premovsurframe_reac(temp_data)
   
    if_table_edited(TRUE)
   
  })
  
  
  observeEvent(input$pre_mov_frame_cell_edit,{
    
    
    proxy = dataTableProxy("pre_mov_frame")
    temp_data=tableedit_premovsurframe_reac()
    
    info = input$pre_mov_frame_cell_edit
   
    if(is.numeric(temp_data[info$row,info$col+1])){
      temp_data[info$row,info$col+1]=as.numeric(info$value)
    }else{
      temp_data[info$row,info$col+1]=info$value
    }
    
   
    #storing edited values into a temporary reactive frame
    tableedit_premovsurframe_reac(temp_data)
    
    if_table_edited(TRUE)
    testselcol(info$col)
  })
  
  observeEvent(input$oral_sens_frame2_cell_edit,{
    
    proxy = dataTableProxy("oral_sens_frame2")
    temp_data=oral_fluid_sensframe_reac()
    
    info = input$oral_sens_frame2_cell_edit
    if(is.numeric(temp_data[info$row,info$col+1])){
    temp_data[info$row,info$col+1]=as.numeric(info$value)
    }else{
      temp_data[info$row,info$col+1]=info$value 
    }
    oral_fluid_sensframe_reac(temp_data)
   
  })
  
  
  observeEvent(input$stateftbl,{
    
    
    inp_tbl1_disstate(which(pstate_list_names==input$stateftbl))
   
    
  })
  
  observe({
    tto=as.numeric(input$txtinp_stop1)
    tfrom=as.numeric(input$txtinp_startday1)
    tby=as.numeric(input$txtinp_by1)
    
    if(!is.na(tfrom)&!is.na(tto)&!is.na(tby)&
       
       tto>=tby&tto>=tfrom ){
    tbl_outdays_reac(seq(from =  tfrom,
                         to = tto,
                         by = tby))
    
    }
  })
  
  observe({
    tto=as.numeric(input$txtinp_stop1mm)
    tfrom=as.numeric(input$txtinp_startday1mm)
    tby=as.numeric(input$txtinp_by1mm)
    
    if(!is.na(tfrom)&!is.na(tto)&!is.na(tby)&
       
       tto>=tby&tto>=tfrom ){
      tbl_outdays_reacmm(seq(from =  tfrom,
                           to = tto,
                           by = tby))
      
    }
  })
 
  
  oralfileloadreac<-reactiveVal(FALSE)
  
  observeEvent(input$FileInp_loadoralsens,
               {
                 
                 file <- input$FileInp_loadoralsens
                 ext <- tools::file_ext(file$datapath)
                 
                 req(file)
                 oralfileloadreac(TRUE)
                 output$error_msg <- renderText({ 
                 if(ext != "csv"){
                  HTML("<span style=color:red><b>Please load a CSV file with rope sampling sensitvity</b></span>")
                 }
                   
                          
                 })
                
                 
                 validate(need(ext == "csv", "Please upload a csv file"))
                
                 if(!is.null(file$datapath)){
                   tempdf<-read.csv(file =file$datapath,header = TRUE,stringsAsFactors = FALSE)
                  
                   
                   if(!is.null(tempdf)){
                     if(is.data.frame(tempdf)){
                       colnames(tempdf)<-c("Prevelence",  "Sensitivity")
                       oral_fluid_sensframe_reac(tempdf) 
                       
                     }}
                   
                  
                   
                 }
                 
               }          
  )
  
  
  observeEvent(input$txtinp_numiter,{
 
    tflag<-TRUE
    if(is.na(input$txtinp_numiter)){
      
      
    } else if(as.numeric(input$txtinp_numiter)>2000){
      
      
      tflag<-FALSE
      showModal(modalDialog(
        
        title = div(HTML("<b>Number of iterations has to be less than 2000</b>"),style = "color:maroon"),
        easyClose = TRUE,
        footer = modalButton("Close")
      )) 
      
      updateNumericInput(session = session,inputId = "txtinp_numiter",value = 2000 )
      
    }   
    
    
    
   
   
    validate(need(tflag, "The number of iterations should be between 500 and 2000"))
  }
  
  )
  
  
  
  
  #end loading an oral sensitvity file
  observeEvent(input$FileInpLoadpremov,
               {
                 
                
                 file <- input$FileInpLoadpremov
                 ext <- tools::file_ext(file$datapath)
                 
                 req(file)
                 output$error_msgpre <- renderText({ 
                   if(ext != "csv"){
                     HTML("<span style=color:red><b>Please load a CSV file with surveillance protocol</b></span>")
                   }
                   
                   
                 })
                 
                 
                 
                 
                 
                 
                 validate(need(ext == "csv", "Please upload a csv file"))
                 
               
                 if(!is.null(file$datapath)){
                   tempdf<-read.csv(file =file$datapath,header = TRUE,stringsAsFactors = FALSE)
                  
                   
                   
                   
                   if(!is.null(tempdf)){
                     if(is.data.frame(tempdf)){
                       
                       tempdf$Premove.survdays<-as.character(tempdf$Premove.survdays)
                       tempdf$Premov.sample.size<-as.character(tempdf$Premov.sample.size)
                       tempdf$Oral.surv.days<-as.character(tempdf$Oral.surv.days)
                       
                       
                       
                       templist<-get_premov_survlist(tempdf)
                       
                       
                       req(templist)
                       tflag=check_premovsurvlist(templist =templist )
                       
                       
                       
                       output$error_msgprecheck <- renderText({ 
                         if(!tflag){
                           HTML("<span style=color:red><b>Please check surveillance protocols in the file being loaded</b></span>")
                         }
                         
                         
                       })
                       validate(need(tflag, "Please check loaded protocols")) 
                       
                       premovsurframe_reac(tempdf)
                       tableedit_premovsurframe_reac(tempdf)
                       if_table_edited(FALSE)
                     }}
                   
                  
                 }
                 
               }          
  )
  
  

##detcurvebutton
  
  observeEvent(input$rundet, {
    
    progress <- shiny::Progress$new(min=0, max=100)
    progress$set(message = "Generating detection curve", value = 0)
    on.exit(progress$close())
    
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    
    tdetstart=as.numeric(input$txtinp_startdetdays)
    tdetend=as.numeric(input$txtinp_enddetdays)
    tdetby=as.numeric(input$txtinp_enddetby)
  
    predaysflag<-TRUE
    if(tdetend>35){
      predaysflag<-FALSE
      showModal(modalDialog(
        
        title = div(HTML("<b>End day has to be less than 35</b>"),style = "color:maroon"),
        easyClose = TRUE,
        footer = modalButton("Close")
      )) 
    }
    if(  tdetstart> tdetend){
      predaysflag<-FALSE
      showModal(modalDialog(
        
        title = div(HTML("<b>Please double check start and end days</b>"),style = "color:maroon"),
        easyClose = TRUE,
        footer = modalButton("Close")
      )) 
    }
    validate(need(predaysflag, "Please check days"))
    
    
    #reading in sim_parms
    temp_sim_parms=sim_outputreactive()$in_sim_parms
    
    #temp_sim_parms$move_day_bounds=c(premov_start_reac(),premov_end_reac())
    temp_sim_parms$normsickparms=inp_normclin_reacvec()
    temp_sim_parms$propsickobs<-as.numeric(input$sickidprob)/100
    temp_sim_parms$PMIPdur<- 0
    temp_sim_parms$PMIPeff<-premov_eff_reactive()
    isolate(pensizearray_reac(rep(round(temp_sim_parms$flock_size_pars/temp_sim_parms$npens),temp_sim_parms$npens)))
    
    if(if_table_edited()){
      premovsurframe_reac(tableedit_premovsurframe_reac())
    }
    pre_mov_sur_list(get_premov_survlist(premovsurframe =premovsurframe_reac() ))
    tflag=check_premovsurvlist(templist =  pre_mov_sur_list())
    if(!tflag)
      showModal(modalDialog(
        
        title = div(HTML("<b>Please double check pre-movement surveillance protocols</b>"),style = "color:maroon"),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    validate(need(tflag, "Please check loaded protocols")) 
    
    
    
    tstraintxt=c("Moderately virulent","Highly virulent olsen","Highly virulent Georgia hu","Highly virulent Olensen long", "Modvir for severeclin")[inp_strain_list_react()]
    txtcontactrate=c("Slow","Fast","Medium", "User defined contact rate")[inp_beta_choice_react()]
    isolate(
    if(input$chbox_selallprot==FALSE){
      input_protocol_vec_reac(as.numeric(strsplit(x = input$txtinp_selected_protocol,split = "," )[[1]]))
    }else{
      input_protocol_vec_reac(1:length(pre_mov_sur_list()))
    }
    )
   
    isolate(
    
      
      temp_survoutframe<-run_and_process_premovsurv_detcurve_shiny(updateprogress =updateProgress,min_day = tdetstart,max_day =  tdetend,by_day = tdetby, 
        premovsurvlist =pre_mov_sur_list(),
                                                                  trans_shin_list =sim_outputreactive(),inscenariosvec =input_protocol_vec_reac(),
                                                                  in_sim_parms = temp_sim_parms,pensizearray=pensizearray_reac(),sens_frame=oral_fluid_sensframe_reac(),contact_rate_string =txtcontactrate ,strain_scenario_string = tstraintxt
      ))
    
   
  
    detc_outreac_frame(   temp_survoutframe)
    is_detcurv_avil(TRUE)
  
    tempdf<-premovsurframe_reac()
    premovsurframe_reac(tempdf)
    
    updateTabsetPanel(session, "inTabset",
                      selected ="Detection curve")
    shinyjs::enable("runsurv")
    shinyjs::enable("rundet")
    })
  
  
  
####Surveillance action buttion####  
  observeEvent(input$runsurv, {
    progress <- shiny::Progress$new(min=0, max=100)
    progress$set(message = "Simulating surveillance", value = 0)
    on.exit(progress$close())
    
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    premov_start_reac(as.numeric(input$txtinp_startpredays))
    premov_end_reac(as.numeric(input$txtinp_endpredays))
    premov_dur_reac(as.numeric(input$txtinp_pmip_len))
    premov_eff_reactive(as.numeric(input$txtinp_pmip_eff)/100)
    
    predaysflag<-TRUE
    if(premov_end_reac()>35){
      predaysflag<-FALSE
      showModal(modalDialog(
        
        title = div(HTML("<b>End day has to be less than 35</b>"),style = "color:maroon"),
        easyClose = TRUE,
        footer = modalButton("Close")
      )) 
    }
    if(premov_start_reac()>premov_end_reac()){
      predaysflag<-FALSE
      showModal(modalDialog(
        
        title = div(HTML("<b>Please double check start and end days</b>"),style = "color:maroon"),
        easyClose = TRUE,
        footer = modalButton("Close")
      )) 
    }
    validate(need(predaysflag, "Please check days"))
    
    
    temp_sim_parms=sim_outputreactive()$in_sim_parms
    temp_sim_parms$move_day_bounds=c(premov_start_reac(),premov_end_reac())
    temp_sim_parms$normsickparms=inp_normclin_reacvec()
    temp_sim_parms$propsickobs<-as.numeric(input$sickidprob)/100
   
    temp_sim_parms$PMIPdur<- premov_dur_reac()
    temp_sim_parms$PMIPeff<-premov_eff_reactive()
    
    isolate(pensizearray_reac(rep(round(temp_sim_parms$flock_size_pars/temp_sim_parms$npens),temp_sim_parms$npens)))
    
    if(if_table_edited()){
      premovsurframe_reac(tableedit_premovsurframe_reac())
    }
    pre_mov_sur_list(get_premov_survlist(premovsurframe =premovsurframe_reac() ))
    
    
    tflag=check_premovsurvlist(templist =  pre_mov_sur_list())
    if(!tflag)
    showModal(modalDialog(

      title = div(HTML("<b>Please double check pre-movement surveillance protocols</b>"),style = "color:maroon"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
    validate(need(tflag, "Please check loaded protocols")) 
    
    tstraintxt=c("Moderately virulent","Highly virulent olsen","Highly virulent Georgia hu","Highly virulent Olensen long", "Modvir for severeclin")[inp_strain_list_react()]
    txtcontactrate=c("Slow","Fast","Medium", "User defined contact rate")[inp_beta_choice_react()]
   
    isolate(
    if(input$chbox_selallprot==FALSE){
      input_protocol_vec_reac(as.numeric(strsplit(x = input$txtinp_selected_protocol,split = "," )[[1]]))
    }else{
      input_protocol_vec_reac(1:length(pre_mov_sur_list()))
    })
   
      isolate(
        temp_survoutframe<-run_and_process_randday_premovsurv_shiny(updateprogress=updateProgress,premovsurvlist =pre_mov_sur_list(),
                                                                    trans_shin_list =sim_outputreactive(),inscenariosvec =input_protocol_vec_reac(),
                                                                    in_sim_parms = temp_sim_parms,pensizearray=pensizearray_reac(),sens_frame=oral_fluid_sensframe_reac(),contact_rate_string =txtcontactrate ,strain_scenario_string = tstraintxt
        ))
      
   

      temp_survoutframe<-temp_survoutframe[,c(5,1,2,4,6,13,20)]
    
  
      
        colnames(temp_survoutframe)<-c("Protocol",
            "ASF strain"  ,            "Transmission rate\n scenario" ,            "Individual sample \n prioritization" ,
      "Detect fraction"   ,   "Infectious pigs at movement if not detected", "Sample size"   
        )
      
 
        
        temp_survoutframe$`Test days prior to movement`<-premovsurframe_reac()$Premove.survdays
        temp_survoutframe$`Mortality trigger per 1000`<-premovsurframe_reac()$Mortrigper1000
        
    surv_out_reactive(temp_survoutframe)
    updateTabsetPanel(session, "inTabset",
                      selected ="Surveillance results"
    )
    isolate(sim3parms(temp_sim_parms))
    
    shinyjs::enable("runsurv")
    shinyjs::enable("rundet")
  })
  observeEvent(input$Edit_protocol, {
    updateTabsetPanel(session, "inTabset",
                      selected ="Premove surveillance protocol"
    )
  })
  
 
  
  observeEvent(input$Updatemorb, {
    
    tcounter<-as.numeric(initplotreactval())+1
    initplotreactval(tcounter)
    prevtrans_listreactive(sim_outputreactive())
   
    tempsim_parms<-sim2parms()
    tempsim_parms$normsickparms=inp_normclin_reacvec()
    
    nflag<-TRUE
    if(tempsim_parms$normsickparms[1]>=tempsim_parms$normsickparms[2]|
       tempsim_parms$normsickparms[2]>=tempsim_parms$normsickparms[3]){
      showModal(modalDialog(
        
        title = div(HTML("<b>Check lower central and high values for normal morbidity</b>"),style = "color:maroon"),
        easyClose = TRUE,
        footer = modalButton("Close")
      )) 
      
      nflag<-FALSE
    }
    validate(need(nflag, "maximum iterations is 2000"))
    
    tempsim_parms$propsickobs<-as.numeric(input$sickidprob)/100
    tempsim_parms$use_high_mortdata=input$chkinp_usehighmort
    if(input$inp_ASFstrain=='5'){
      tempsim_parms$propsickobs=1
    }
 
    progress <- shiny::Progress$new(min=0, max=100)
    progress$set(message = "Simulating normal mortality\n morbidity", value = 50)
    on.exit(progress$close())
    
    tpensizearray<-rep(round(tempsim_parms$flock_size_pars/tempsim_parms$npens),tempsim_parms$npens)
    
    sim_normmort <- sim_normmort_func(sim_parms  =tempsim_parms,pensizearray = tpensizearray )
    sim_normsick <- sim_normsick_func(transmission_results = sim_outputreactive()$transmission_results, in_sim_parms = tempsim_parms,pensizearray = tpensizearray)
    

    trans_outlist<-sim_outputreactive()
    trans_outlist$sim_normmort=sim_normmort
    trans_outlist$sim_normsick=sim_normsick
    trans_outlist$in_sim_parms=tempsim_parms
    sim2parms(tempsim_parms)
    sim_outputreactive(trans_outlist)
    
    if(should_trans_runbeforesurv()==FALSE){
      
    shinyjs::enable("runsurv")
    shinyjs::enable("rundet")
    outtextreactive("")
    need_run_normmorb_atleast(FALSE)
    }else{
      need_run_normmorb_atleast(FALSE)
      outtextreactive("Please run transmission before running surveillance") 
      
    }
  
  })
  
  
  ####Run transmission####
  observeEvent(input$do, {
    onlyupdatenormortmob(FALSE)
    
    
    updateTabsetPanel(session, "inTabset",
                      selected ="Transmission plots"
    )
    prevtrans_listreactive(sim_outputreactive())

    tempsim_parms<-sim2parms()
    ttflag<-TRUE
    if (as.numeric(input$txtinp_numiter)<500) {
      
      ttflag<-FALSE
      showModal(modalDialog(
        
        title = div(HTML("<b>Number of iterations has to be greater than or equal to 500</b>"),style = "color:maroon"),
        easyClose = TRUE,
        footer = modalButton("Close")
      )) 
      
      updateNumericInput(session = session,inputId = "txtinp_numiter",value = 500 )
      
      
      
    }
    
    validate(need(ttflag, "The number of iterations should be between 500 and 2000"))

    
    
    
    tempsim_parms$num_iterations<-as.numeric(input$txtinp_numiter)
    
    
    
    
    
    tempsim_parms$flock_size_pars=c(inp_total_herdsize()) #1200 4800 Total number of pigs on premises
    tempsim_parms$npens<-inp_total_herdsize()/inp_pensizereac()
    tempsim_parms$normsickparms=inp_normclin_reacvec()
    
    nflag<-TRUE
    if(tempsim_parms$normsickparms[1]>=tempsim_parms$normsickparms[2]|
       tempsim_parms$normsickparms[2]>=tempsim_parms$normsickparms[3]){
      showModal(modalDialog(
        
        title = div(HTML("<b>Check lower central and high values for normal morbidity</b>"),style = "color:maroon"),
        easyClose = TRUE,
        footer = modalButton("Close")
      )) 
      
      nflag<-FALSE
    }
    validate(need(nflag, "maximum iterations is 2000"))
    tempsim_parms$propsickobs<-as.numeric(input$sickidprob)/100
    tempsim_parms$use_high_mortdata=input$chkinp_usehighmort
    tempsim_parms$between_penpeoplefrac=as.numeric(input$Between_pen_people_percent)/100
   
    if(input$inp_ASFstrain=='5'){
      tempsim_parms$propsickobs=1
    }
    
  
    strain_namereac(ttstrainnamelist[as.numeric(input$inp_ASFstrain)])
    contact_namereac(ttcontactratelist[(temp_beta_choice_react())])
  
    
    sim2parms(tempsim_parms)
    sim3parms(tempsim_parms)
   
    
    barn_pen_arraylist(temp_barn_pen_arraylist())
    
    inp_beta_choice_react(temp_beta_choice_react())
    sim_inpcustom_withinpenB_reacvec(custom_withinpenB_reacvec())
    sim_inpcustom_betweenpenB_reacvec(custom_betweenpenB_reacvec())
  
    inp_strain_list_react(as.numeric(input$inp_ASFstrain))
    tcounter<-as.numeric(initplotreactval())+1
    initplotreactval(tcounter)
    should_trans_runbeforesurv(FALSE)
    
    
      progress <- shiny::Progress$new(min=0, max=100)
      progress$set(message = "Simulating transmission", value = 0)
      on.exit(progress$close())

      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }


      if(as.numeric(initplotreactval())==0){
        
        trans_shin_list<-prevtrans_listreactive()
        
        

      } else{
      
      trans_shin_list=get_single_trans_results_progressbar(updateprogress =updateProgress, in_sim_parms = sim2parms(),in_strain = inp_strain_list_react(),in_beta_choice = inp_beta_choice_react(),inbetavec = sim_inpcustom_withinpenB_reacvec(),in_barn_pen_arraylist = barn_pen_arraylist(),
                                                           inbetaBbvec = sim_inpcustom_betweenpenB_reacvec()
                                                           )
      
      
      sim_outputreactive(trans_shin_list)
      }
    
    shinyjs::enable("runsurv")
    shinyjs::enable("rundet")
    need_run_normmorb_atleast(FALSE)
    should_trans_runbeforesurv(FALSE)
    outtextreactive("")
      })
  
  outinftablereactive<-reactive({
    
    if(input$chbox_comp==FALSE){
      outtable=getdatatable_infectious_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate = 1,transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
    } else {
      outtable= Comb_getdatatable_infectious_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate = 4,
                                             transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel,
                                             past_transmission_results_premlevel = prevtrans_listreactive()$transmission_results_premlevel)
    }
    outtable
  }
    
  )
  
  
  outclintablereactive<-reactive({
    
    if(input$chbox_comp==FALSE){
    outtable=getdatatable_clin_prem(outdays = tbl_outdays_reacmm(),sim_parms = sim2parms(),sim_normsick = sim_outputreactive()$sim_normsick,
                                    transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
    } else {
    
      outtable=getdatatable_clin_prem(outdays = tbl_outdays_reacmm(),sim_parms = sim2parms(),sim_normsick = sim_outputreactive()$sim_normsick,
                                      transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
      
       }
    outtable
  }
  
  )
  
  outdailmorttablereactive<-reactive({
    
    if(input$chbox_comp==FALSE){
      outtable=getdatatable_dead_prem(outdays = tbl_outdays_reacmm(),sim_parms = sim2parms(),sim_normdead = sim_outputreactive()$sim_normmort,
                                      transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
    } else {
      
    
      outtable=getdatatable_dead_prem(outdays = tbl_outdays_reacmm(),sim_parms = sim2parms(),sim_normdead = sim_outputreactive()$sim_normmort,
                                      transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
      }
    outtable
  }
  
  )
  
  outtextreactive<-reactiveVal("")
  
 
   outdisstatebasicreactive<-reactive({
   
    if(input$chbox_comp==FALSE){
     
      outtable=get_datatbl_state_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate =inp_tbl1_disstate() ,transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
    } else {
      outtable= Comb_tbl_state_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate = inp_tbl1_disstate(),
                                                  transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel,
                                                  past_transmission_results_premlevel = prevtrans_listreactive()$transmission_results_premlevel)
    }
    outtable
  }
  
  )
  
  genericwritefile<-function(inobject,tfilename){
    write(x = inobject,file=tfilename)
    
  }
 
  output$downloadstatetable <- downloadHandler(
    
    filename = function(){
      
      tempname=paste(pstate_list_names[inp_tbl1_disstate()], "tbl.csv",sep = " ")
      tempname}, 
    content = function(fname){
      write.csv(outdisstatebasicreactive(), fname)
    }
  )
  
  output$downloadinftbl <- downloadHandler(
    filename = function(){"Infectious table.csv"}, 
    content = function(fname){
      write.csv(outinftablereactive(), fname)
    }
  )
 
  
  output$downloadoralsens <- downloadHandler(
    filename = function(){"oral_sens_frameapp.csv"}, 
    content = function(fname){
      write.csv(oral_fluid_sensframe_reac(), fname)
    }
  )
  
  
  output$downloadclin <- downloadHandler(
    filename = function(){"Clinical signs.csv"}, 
    content = function(fname){
      write.csv(outclintablereactive(), fname)
    }
  )
  
  output$downloaddailmort <- downloadHandler(
    filename = function(){"Daily mort.csv"}, 
    content = function(fname){
      write.csv(outdailmorttablereactive(), fname)
    }
  )
  
  output$download <- downloadHandler(
    filename = function(){"Infectious_pigs.csv"}, 
    content = function(fname){
      write.csv(outinftablereactive(), fname)
    }
  )
  
 ####Download tough#### 
  output$downdetres <- downloadHandler(
  
    
    filename = function() {
      paste("detoutput", "zip", sep=".")
    },
    
    
    content = function(fname) {
   
      olddir<-getwd()
      tmpdir <- tempdir()
      setwd(tempdir())
      
      fs <- c()
      path <- "sim_parms.RData"
      fs <- c(fs, path)
      tempobj=sim2parms()
      save(tempobj,file =  path)
      
      path <- paste0("simparms.txt")
      fs <- c(fs, path)
      sink( path)
      print(as.data.frame.complex(sim3parms()))
      print(paste0("Contact rate scenario: ",contact_namereac()))
      print(paste0("Strain scenario: ",strain_namereac()))
      sink()
      
      
      
      
      path <- paste0("infectiousovertime.csv")
      fs <- c(fs, path)
      write.csv(x =  outinftablereactive(), file = path)
      
      
      path <- paste0("sickpigs.csv")
      fs <- c(fs, path)
      write.csv(x =  outclintablereactive(), file = path)
      
      path <- paste0("pred_mortality.csv")
      fs <- c(fs, path)
      write.csv(x =  outdailmorttablereactive(), file = path)
      
      
      if(is_detcurv_avil()){
      path <- paste0("det_results_table.csv")
      fs <- c(fs, path)
      write.csv(x =  detc_outreac_frame(), file = path)
       
      }
      
      zip(zipfile=fname, files=fs)
    
      setwd(olddir)
      
    },
    contentType = "application/zip"
  )
  
  
  output$downsurvres <- downloadHandler(

    
    filename = function() {
      paste("output", "zip", sep=".")
    },
    
    
    content = function(fname) {
     
      olddir<-getwd()
      tmpdir <- tempdir()
      setwd(tempdir())
     
      fs <- c()
        path <- "sim_parms.RData"
        fs <- c(fs, path)
        tempobj=sim2parms()
        save(tempobj,file =  path)
        
        path <- paste0("simparms.txt")
        fs <- c(fs, path)
        sink( path)
        print(as.data.frame.complex(sim3parms()))
        print(paste0("Contact rate scenario: ",contact_namereac()))
        print(paste0("Strain scenario: ",strain_namereac()))
        sink()
        
        
        
        
        path <- paste0("infectiousovertime.csv")
        fs <- c(fs, path)
        write.csv(x =  outinftablereactive(), file = path)
        
          
          
        tempname=paste(pstate_list_names[inp_tbl1_disstate()], "tbl.csv",sep = " ")
        path <- paste0(tempname)
        fs <- c(fs, path)
        write.csv(x =  outdisstatebasicreactive(), file = path)
        
       
        if(input$chbox_comp==FALSE){
         
          touttable=get_datatbl_state_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate =1 ,transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
        } else {
          touttable= Comb_tbl_state_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate = 1,
                                        transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel,
                                        past_transmission_results_premlevel = prevtrans_listreactive()$transmission_results_premlevel)
        }
        
        path <- paste0("Susceptible pigs.csv")
        fs <- c(fs, path)
        write.csv(x =  touttable, file = path)
        
        ###recovered
        
        if(input$chbox_comp==FALSE){
          
          touttable=get_datatbl_state_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate =6 ,transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
        } else {
          touttable= Comb_tbl_state_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate = 6,
                                         transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel,
                                         past_transmission_results_premlevel = prevtrans_listreactive()$transmission_results_premlevel)
        }
        
        path <- paste0("recovered.csv")
        fs <- c(fs, path)
        write.csv(x =  touttable, file = path)
        
        ##recovered
        
        
        if(input$chbox_comp==FALSE){
          
          touttable=get_datatbl_state_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate =7 ,transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
        } else {
          touttable= Comb_tbl_state_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate = 7,
                                         transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel,
                                         past_transmission_results_premlevel = prevtrans_listreactive()$transmission_results_premlevel)
        }
        
        path <- paste0("cumulativedead.csv")
        fs <- c(fs, path)
        write.csv(x =  touttable, file = path)
        
        ### dead cumulative
        ####latent
        ##recovered
        ###dead cumulative
        
        if(input$chbox_comp==FALSE){
          
          touttable=getdatatable_latent_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate = 2,transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
          
           } else {
          touttable= Comb_getdatatable_latent_prem(outdays = tbl_outdays_reac(),sim_parms = sim2parms(),instate = 2,
                                         transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel,
                                         past_transmission_results_premlevel = prevtrans_listreactive()$transmission_results_premlevel)
        }
        
        path <- paste0("latent.csv")
        fs <- c(fs, path)
        write.csv(x =  touttable, file = path)
        
        
       
        path <- paste0("sickpigs.csv")
        fs <- c(fs, path)
        write.csv(x =  outclintablereactive(), file = path)
        
        path <- paste0("pred_mortality.csv")
        fs <- c(fs, path)
          write.csv(x =  outdailmorttablereactive(), file = path)
      
          
        
        path <- paste0("Surv_results_table.csv")
        fs <- c(fs, path)
        write.csv(x =  surv_out_reactive(), file = path)
            
            
        
        zip(zipfile=fname, files=fs)
        
        setwd(olddir)
       
    },
    contentType = "application/zip"
  )
  
  
  
  savefilename = "/SDA_request_7192021/severeclinsigns.csv"
  
  output$graph2<-renderPlot({
    if(input$chbox_comp==FALSE){
    Visualize_geom_ribbon_infectious_prem(sim_parms = sim2parms(),instate = 1,transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
    } else {
      Comb_Visualize_geom_ribbon_infectious_prem(sim_parms = sim2parms(),instate = 4,
                                            transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel,
                                            past_transmission_results_premlevel = prevtrans_listreactive()$transmission_results_premlevel)
  }})
  
  
  output$graph7<-renderPlot({
    if(is_detcurv_avil()){
      plot_geom_detcurve(survscenstring = "1",in_detcurveoutframe = detc_outreac_frame())
    }
    
    })
  
  output$ropesenplot<-renderPlot({
    get_oral_sens_plot(sens_frame = oral_fluid_sensframe_reac())
  })
  
  output$graph1<-renderPlot({
    if(input$chbox_comp==FALSE){
      Visualize_geom_ribbon_state_prem(sim_parms = sim2parms(),instate = inp_graph1_disstate(),transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel)
    } else {
     
      Comb_Visualize_geom_ribbon_state_prem(sim_parms = sim2parms(),instate = inp_graph1_disstate(),
                                            transmission_results_premlevel = sim_outputreactive()$transmission_results_premlevel,
                                              past_transmission_results_premlevel = prevtrans_listreactive()$transmission_results_premlevel)
    }
    
    
    })
  
  
  output$Dtbl_1inf<-renderDataTable(server = FALSE, {
    
    DT::datatable(outinftablereactive(),rownames = NULL,
                  extensions = 'Buttons',  
    options = list(  searchable = FALSE,
                    scrollY = TRUE,   ## enable scrolling on Y axis
                    autoWidth = TRUE,
                    columnDefs = list(list(className = 'dt-center', targets="_all")),
                    dom = 'rtip',
                 
                    exportOptions = list(
                      modifier = list(page = "all"))
                   
                                        
                    )
                    
    
    
    )})
  
  
  output$clintextoutput<-renderText({
   
      HTML(paste0(clin_textout_reactive()))
   
    
  })
  
  output$Dtbl_anystate<-renderDataTable({
    
    DT::datatable(outdisstatebasicreactive(),rownames = NULL, 
                  extensions = 'Buttons',  
                  options = list(
                    searchable = FALSE,
                    scrollY = TRUE,   ## enable scrolling on Y axis
                    autoWidth = TRUE,
                    columnDefs = list(list(className = 'dt-center', targets="_all")),
                    dom = 'rtip'
                    
                  )
                  
                  
                  
    )})
  
  output$Dtbl_2clin<-renderDataTable({
   
    DT::datatable(outclintablereactive(),rownames = NULL,
                  extensions = 'Buttons',  
                  options = list(  searchable = FALSE,
                                   scrollY = TRUE,   
                                   autoWidth = TRUE,
                                   columnDefs = list(list(className = 'dt-center', targets="_all")),
                                   dom = 'rtip'
                                 
                  ))
    
                  
                  
                  
    })
  #premov
  testselcol<-reactiveVal(NULL)
  
 output$pre_mov_frame<-renderDataTable(server = FALSE,{
   
   
   
   tdf<-premovsurframe_reac()
   colnames(tdf)<-as.character(infodatframe[2,])
 

   for(i in 1:ncol(tdf)){
     runjs(createDiv(i, "TOOLTIP"))
     runjs(fillDiv(tdf, i,infodatframe, "TOOLTIP"))
   }
   
   
   
   
   headerCallback <- c(
     "function(thead, data, start, end, display){",
     "  var ncols = data[0].length;",
     tooltips(ncol(tdf), "TOOLTIP"),
     "  for(var i = 1; i < ncols; i++){",
     "    $('th:eq(' + i + ')', thead).qtip(tooltips[i-1]);",
     "  }",
     "}"
   )

   DT::datatable(tdf,rownames = NULL,
                 extensions = 'Buttons',
                 editable = 'cell',callback = JS('table.columns.adjust().draw();'),
                  
                 options = list(  searchable = FALSE,
                                  scrollY = TRUE,   ## enable scrolling on Y axis
                                  scrollX=TRUE,
                                  autoWidth = TRUE,
                                  columnDefs = list(list(className = 'dt-center', targets="_all")),
                                  dom = 'Bfrtip',
                                 headerCallback = JS(headerCallback),
                                  buttons = c('csv', 'excel')
                                  
                 ))
  
    }
   
   )
 
 
 output$oral_sens_frame<-renderDataTable(server = FALSE,{
  
   tdf<-oral_fluid_sensframe_reac()
   colnames(tdf)<-c("Within Pen Prevalence","Oral fluids sensitivity")
   DT::datatable(tdf,rownames = NULL,
                 
                 extensions = 'Buttons',  
                 editable = 'cell',
                 options = list(  searchable = FALSE,
                                  scrollY = TRUE,   ## enable scrolling on Y axis
                                  scrollX=TRUE,
                                  autoWidth = TRUE,
                                  columnDefs = list(list(className = 'dt-center', targets="_all")),
                                  dom = 'Bfrtip',
                                  
                                  buttons = c('csv', 'excel')
                                  
                                  
                                  
                 ))})
 
 
 output$oral_sens_frame2<-renderDataTable(server = FALSE,{
   tdf<-oral_fluid_sensframe_reac()
   
   colnames(tdf)<-c("Within pen prevalence","Oral fluids sensitivity")
   DT::datatable(tdf,rownames = NULL,
                 
                 extensions = 'Buttons',  
                 editable = 'cell',
                 options = list(  searchable = FALSE,
                                  scrollY = TRUE,   ## enable scrolling on Y axis
                                  scrollX=FALSE,
                                  autoWidth = TRUE,
                                  columnDefs = list(list(className = 'dt-center', targets="_all")),
                                  dom = 'Bfrtip',
                                  
                                  buttons = c('csv', 'excel')
                                  
                                  
                 ))})
 
 ####detoutrender####
 
 output$Deccurvoutput<-renderDataTable(server = FALSE,{
   
   if(is_detcurv_avil()){
   
    tdf<-detc_outreac_frame()[,c(8,1,2,4,5,7,23)]
   
    colnames(tdf)<-c("Protocol","Testday\npost exposure" ,"Detect.fraction","Strain"         ,     "Contact.rate"  ,  "Individual sample \npriority", "Sample size")
                     
  
   DT::datatable(tdf,rownames = NULL,
                 extensions = 'Buttons',  
                 editable = 'cell',
                 options = list(  searchable = FALSE,
                                  scrollY = TRUE,   ## enable scrolling on Y axis
                                  scrollX=TRUE,
                                  autoWidth = TRUE,
                                  columnDefs = list(list(className = 'dt-center', targets="_all")),
                                  dom = 'Bfrtip',
                                  
                                  buttons = c('csv', 'excel')
                                  
                                  
                                  
                 ))}
 
 }  
 )
 
 
 ####surv_output#### 
 surv_out_reactive<-reactiveVal()
 surv_out_reactive(data.frame())
 
 output$text12<-renderUI(
   {
    
     
     tt2<-outtextreactive() 
     
   HTML(paste("<h5 style = color:darkblue;><b>",tt2,"</b><h5>"))
     
     
   })
 output$pdf_viewer <- renderUI( tags$iframe(src = 'survmanual3.pdf',width=800, height = 800) ) 
 
 getPage<-function() {
   return(includeHTML("Instruction Manual surv.html"))
 }
 output$manual2 <- renderUI({

  getPage()
   
 })
 
 output$survresout<-renderDataTable(server = FALSE,{

   tdf<-surv_out_reactive()

   DT::datatable(tdf,rownames = NULL,
                 extensions = 'Buttons',  
                 editable = 'cell',
                 options = list(  searchable = FALSE,
                                  scrollY = TRUE,  
                                  scrollX=TRUE,
                                  autoWidth = TRUE,
                                  columnDefs = list(list(className = 'dt-center', targets="_all")),
                                  dom = 'Bfrtip',
                                  
                                  buttons = c('csv', 'excel')
                                  
                 ))}
   )  
   
 
  output$Dtbl_mort<-renderDataTable(server = FALSE,{
   
    DT::datatable(outdailmorttablereactive(),rownames = NULL,
                  extensions = 'Buttons',  
                  options = list(  searchable = FALSE,
                                   scrollY = TRUE,   
                                   autoWidth = TRUE,
                                   columnDefs = list(list(className = 'dt-center', targets="_all")),
                                   dom = 'rtip'
                                  
                                   
                                
                  ))
    
  })
  
  
}


tags$head(tags$script('document.title = "New Tab Title";'))

shinyApp(ui=ui,server=server)

