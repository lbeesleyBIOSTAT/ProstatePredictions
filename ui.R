
 shiny::fluidPage(
	  
	  
	  shiny::titlePanel("Outcome Prediction for Prostate Cancer Patients"),

	  shiny::sidebarLayout(
		 shiny::sidebarPanel(
		 	fluidRow(
		 	column(2,
				img(src="2color-bluebg.png",height=100,width=100)
			),
            column(10,
			print("This application provides outcome predictions for prostate cancer patients. Enter patient characteristics at the time 
			of initial treatment. 'CF' stands for metastatic clinical failure, defined as the finding of metastatic disease in 
			imaging. Baseline is defined as date of treatment (surgery or radiation).\n") ,
		 	a("\n Click here for more information about the model and how to interpret the plots. ",target="_blank",href="RShinyDoc.pdf")
             )
			),
		 	hr(),
		 	fluidRow(
		 	column(6,
			shiny::selectInput("plottype", "Plot Type:",
                c("State Occupancy Probabilities" = "S",
                  "Overall Survival Probability" = "OS",
                  "Metastatic-Free Survival Probability" = "EF")),
            
            shiny::sliderInput('maxtime', "Years from baseline to plot (up to 15): ",  min = 1, max = 15, value = 10),
            shiny::checkboxInput("plotNums", label = "Show Probabilities", value = FALSE) ,
            	shiny::conditionalPanel(condition="input.plottype !='HAZ'",
                shiny::sliderInput('CurTime', "Patient known to be alive at time (years): ", value = 0, min = 0, max = 15) ,
				shiny::checkboxInput("RecurEvent", label = "Patient had observed clinical failure", value = FALSE) ),
	   		shiny::conditionalPanel(condition="input.RecurEvent==true",
                shiny::sliderInput('RecurTime', "When did the clinical failure occur (years)?: ", min = 0, max=15, value = 0) )
               
            ),
            column(3,
            shiny::selectInput("tx", "Treatment:", c("Surgery", "Radiation")),
            
            shiny::selectInput("TStage", "Clinical T Stage:", c("T1","T2","T3"), selected = 'T2'),
            shiny::numericInput("basePSA", label = "Pre-Treatment* PSA (ng/mL)", value= 8, min = 0.2, max = 500, step = 1), 
           	shiny::numericInput("glandvol", label = "Gland Volume (mL)", value = 40), 
            shiny::selectInput("pni", "Perineural Invasion:", c("No","Yes"), selected = 'Yes')
             ),
            column(3,
            shiny::selectInput("gleason", "Gleason Score:", c("5-6", "7=3+4", "7=4+3", "8", "9-10"), selected = "8"),
           	shiny::numericInput("age", label = "Age",value = 65),
			shiny::selectInput("race", "Race:", c("White" = "W", "African American" = "NW", "Other" = "O"), selected = 'White'),
			shiny::selectInput("charlson", "Charlson Index:", c("0", "1", "2", "3+"), selected = '1')	
             )
			) ,
			hr(),
			print("*Value of PSA before the first treatment (neoadjuvant hormones if received, surgery, or radiation).\n") ,
			hr(),
    			print("Disclaimer: This model has not been validated and is intended to be used by clinicians and researchers.\n Contact: Lauren J Beesley, lbeesley@umich.edu") 
         , width = 7),
	
	    shiny::mainPanel(shiny::uiOutput("plots")	 , width = 5)	    
         
  
	  	)#end sidebarlayout
	)