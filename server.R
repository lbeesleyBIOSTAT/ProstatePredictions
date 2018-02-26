
 function(input, output ) {
	plotInput <- shiny::reactive({
		plotname <- 'Plot1'
		plot_output_list <-  list(shiny::plotOutput(plotname, inline=TRUE))
		do.call(tagList, plot_output_list)
	})			        
	 local({
    			plotname <- 'Plot1'
		output[[plotname]]   <- renderPlot({
			times = seq(from=0,to=as.numeric(input$maxtime),length.out = 15)
				times = sort(unique(c(times, c(0:input$maxtime))))
			Prob1Save = rep(0,length(times))
			Prob2Save = rep(0,length(times))
			Prob3Save = rep(0,length(times))
			Prob4Save = rep(0,length(times))	
			cuts = c(5,10,12.5)
			TRANS_TR = function(x){return((x-7)/15)}

			# if(input$RecurTime > input$CurTime & input$RecurEvent == TRUE){
					# stop('Clinical failure time must be at or before current time')
			# }
			validate(
				need(input$basePSA >= 0.01 & input$basePSA <= 500, 'PSA must be between 0.01 and 500'),
				need(input$age <= 85 & input$age >= 40, 'Age must be between 40 and 85'),
				need(input$glandvol >= 2 & input$glandvol <= 500, 'Gland Volume must be between 2 and 500'),
				need(input$RecurEvent == FALSE | (input$RecurTime <= input$CurTime & input$RecurEvent == TRUE), 'Clinical failure time must be at or before current time')
			)
			# if(input$basePSA < 0.01 | input$basePSA > 500){
					# stop('PSA must be between 0.01 and 500')
			# }
			# if(input$age > 80 | input$age < 40){
					# stop('Age must be between 40 and 85')
			# }
			# if(input$glandvol < 2 | input$glandvol > 500){
					# stop('Gland Volume must be between 2 and 500')
			# }			

			age = (input$age-60)/10
			tx = c(as.numeric(input$tx == 'Surgery'), as.numeric(input$tx == 'Radiation'))
			stage = c(as.numeric(input$TStage == 'T1'), as.numeric(input$TStage == 'T2'), as.numeric(input$TStage == 'T3'))
			comorbid = c(as.numeric(input$charlson == '1'), as.numeric(input$charlson == '0'), as.numeric(input$charlson == '2'), 
						as.numeric(input$charlson == '3+'))
			gleason = c(as.numeric(input$gleason == '7=4+3'), as.numeric(input$gleason == '5-6'), as.numeric(input$gleason == '7=3+4'), 
						as.numeric(input$gleason == '8'), as.numeric(input$gleason == '9-10'))
			pni = c(as.numeric(input$pni == 'No'), as.numeric(input$pni == 'Yes'))
			txyeargrp = c(0,0,1)#estimations for latest treatment year group
			logbasepsa = log(input$basePSA+0.0001) - 1.8 
			glandvol = log(input$glandvol+0.0001) - 3.6
			nonwhite = as.numeric(input$race == 'NW')

			
						theta=c(4.8709507311,2.4443406454,1.0917176159,1.0905351584,0.0051900386,0.0106892087,0.0257787826,0.0398417152,0.4166370924,
-0.1593926399,0.3875586805,-0.2236755854,0.1655972708,-0.5504278148,-0.6368917356,0.0970926747,0.3224947993,0.0157813935,
0.1439785086,0.7917331358,0.6409081723,0.5374207719,-1.9168444669,-0.0809351616,-1.0924832180,0.0004624098,0.8200872011,
0.0531704936,0.1374205618,1.5709207771,0.3032139002,0.5448889092,0.3845889406,0.4149689801,0.0584384732,0.2571316880,
0.2072396009,0.1120528523,0.0511385219,0.8669891248,0.1786899780,0.3947249170,-0.0458645464,0.1284344321,0.7767575958,
0.1310879270,-0.9506155263,-0.1557320700,-0.3571915055,-0.2429476526,-0.9127202087,-0.0992964812,0.1592985700,-0.5085043517,
0.4729513427,0.1630241167,0.1278362333,-0.0089745664,-0.3717502488,0.0449099100,-0.3292160293,-0.2864758526,0.0454940471,
0.4915619754,0.3645263499,-0.2917036919,0.1438459806,0.1536374716,-0.0822560599,0.3095071762,-0.0581970662,0.1590141543,
-0.3772074830,-0.2289168325,-0.0467752992,0.0247712755,-0.1124004846,-0.4835542220,-0.5380356358,0.0346875669,0.2281079297,
-0.1760832716,0.4492972070)
												
			
			integrate1_short = function(vr)
			{						
				exp23_temp<-exp(s3_23*(stage[2]+stage[3])+
			    				c2_23*comorbid[2]+c3_23*comorbid[3]+c4_23*comorbid[4]+
			    				g4_23*(gleason[1]+gleason[3]+gleason[4])+g5_23*gleason[5]+
			    				tx2_23*tx[2]+
			    				pni2_23*pni[2]+
			    				age_23*age+
			    				logpsa_23*logbasepsa+
			    				nonwhite_23*nonwhite + glandvol_23*glandvol+
			    				as.numeric(btr23)*TRANS_TR(vr)+
			    				stage_23Int*(stage[2]+ stage[3])*tx[2]+ #group stages 2 and 3 here
			   				comorbid_23Int*(comorbid[1] + 2*comorbid[3] + 3*comorbid[4])*tx[2]+
			   				gleason_23Int*(gleason[3] + gleason[1] + gleason[4] + 2*gleason[5] )*tx[2]+ 
			    				pni2_23Int*pni[2]*tx[2]+
			    				age_23Int*age*tx[2]+
			    				logpsa_23Int*logbasepsa*tx[2]+
			    				nonwhite_23Int*nonwhite*tx[2]+ glandvol_23Int*glandvol*tx[2])   					
			    hazard12_temp = (p12/(l12^p12))*exp12*vr^(p12-1)				
				Cumhazard12_temp = (1/(l12^p12))*exp12*vr^(p12)
				Cumhazard13_temp = base_Haz(vr)*exp13
				Cumhazard23_temp = (1/(l23^p23))*exp23_temp*(t-vr)^(p23)										
				Surv1_temp = exp(-Cumhazard12_temp-Cumhazard13_temp)
				Surv2_temp = exp(-Cumhazard23_temp)				
				return(hazard12_temp*Surv1_temp*(1-Surv2_temp))   
			}
			integrate2_short = function(vr)
			{
				exp23_temp<-exp(s3_23*(stage[2]+stage[3])+
			    				c2_23*comorbid[2]+c3_23*comorbid[3]+c4_23*comorbid[4]+
			    				g4_23*(gleason[1]+gleason[3]+gleason[4])+g5_23*gleason[5]+
			    				tx2_23*tx[2]+
			    				pni2_23*pni[2]+
			    				age_23*age+
			    				logpsa_23*logbasepsa+
			    				nonwhite_23*nonwhite + glandvol_23*glandvol+
			    				as.numeric(btr23)*TRANS_TR(vr)+
			    				stage_23Int*(stage[2]+ stage[3])*tx[2]+ #group stages 2 and 3 here
			   				comorbid_23Int*(comorbid[1] + 2*comorbid[3] + 3*comorbid[4])*tx[2]+
			   				gleason_23Int*(gleason[3] + gleason[1] + gleason[4] + 2*gleason[5] )*tx[2]+ 
			    				pni2_23Int*pni[2]*tx[2]+
			    				age_23Int*age*tx[2]+
			    				logpsa_23Int*logbasepsa*tx[2]+
			    				nonwhite_23Int*nonwhite*tx[2]+ glandvol_23Int*glandvol*tx[2])   					
			    hazard12_temp = (p12/(l12^p12))*exp12*vr^(p12-1)				
				Cumhazard12_temp = (1/(l12^p12))*exp12*vr^(p12)
				Cumhazard13_temp = base_Haz(vr)*exp13
				Cumhazard23_temp = (1/(l23^p23))*exp23_temp*(t-vr)^(p23)										
				Surv1_temp = exp(-Cumhazard12_temp-Cumhazard13_temp)
				Surv2_temp = exp(-Cumhazard23_temp)	
				return(hazard12_temp*Surv1_temp*Surv2_temp)   
			}
			integrate4_short = function(vd)
			{	
			    hazard13_temp = base_haz(vd)*exp13					
				Cumhazard12_temp = (1/(l12^p12))*exp12*vd^(p12)
				Cumhazard13_temp = base_Haz(vd)*exp13						
				Surv1_temp = exp(-Cumhazard12_temp-Cumhazard13_temp)						
				return(hazard13_temp*Surv1_temp)   
			}


      		l12<-exp(theta[1]) #weibull parameters
		    l23<-exp(theta[2])  
		    p12<-theta[3]  #weibull parameters.  
		    p23<-theta[4]     
		   	b1<-theta[5]#1->3 baseline parameters
		   	b2<-theta[6]
		    b3=theta[7]  
		    b4=theta[8] 
		    s2_12=theta[9]  # stage2
		    s2_13=theta[10]
		    s3_12=theta[11]  # stage3
		    s3_13=theta[12]
		    s3_23=theta[13]				    
		    c2_12=theta[14]  # comorbid2
		    c2_13=theta[15]
		    c2_23=theta[16]
		    c3_12=theta[17]  # comorbid3
		    c3_13=theta[18]
		    c3_23=theta[19]
		    c4_12=theta[20]  # comorbid4
		    c4_13=theta[21]
		    c4_23=theta[22]    				    
		    g2_12=theta[23]  # gleason2
		    g2_13=theta[24]
		    g3_12=theta[25]  # gleason3
		    g3_13=theta[26]
		    g4_12=theta[27]  # gleason4
		    g4_13=theta[28]
		    g4_23=theta[29]
		    g5_12=theta[30]  # gleason5
		    g5_13=theta[31]
		    g5_23=theta[32]    
		    tx2_12=theta[33]  # tx2
		    tx2_13=theta[34]
		    tx2_23=theta[35]				    
		    pni2_12=theta[36]  # pni2
		    pni2_13=theta[37]
		    pni2_23=theta[38]				    
		    age_12=theta[39]  # age
		    age_13=theta[40]
		    age_23=theta[41]				    
		    logpsa_12=theta[42]  # logbasepsa
		    logpsa_13=theta[43]
		    logpsa_23=theta[44]   				    
		    txyear2_12=theta[45]  # txyeargrp2
		    txyear2_13=theta[46]
		    txyear3_12=theta[47]  # txyeargrp3
		    txyear3_13=theta[48]			    
		    nonwhite_12=theta[49]  # logbasepsa
		    nonwhite_13=theta[50]
		    nonwhite_23=theta[51] 				
		    glandvol_12=theta[52]  # logbasepsa
		    glandvol_13 =theta[53]
		    glandvol_23 =theta[54] 				    
		    btr23=theta[55]  #trCenter
		    stage_12Int=theta[56]  # stage2
		   	stage_13Int=theta[57]
		    stage_23Int=theta[58] #2 and 3 grouped together for this transition
		    comorbid_12Int=theta[59]  # comorbid2
		    comorbid_13Int=theta[60]
		    comorbid_23Int=theta[61]
		    gleason_12Int=theta[62]  # gleason2
		    gleason_13Int=theta[63]
		    gleason_23Int=theta[64] #group gleason 1, 2, 3, 4, for 23 transition
		    pni2_12Int=theta[65]  # pni2
		    pni2_13Int=theta[66]
		    pni2_23Int=theta[67]
		    age_12Int=theta[68]  # age
		    age_13Int=theta[69]
		    age_23Int=theta[70]
		    logpsa_12Int=theta[71]  # logbasepsa
		    logpsa_13Int=theta[72]
		    logpsa_23Int=theta[73]   
		    txyear2_12Int=theta[74]  # txyeargrp2
		    txyear2_13Int=theta[75]
		    txyear3_12Int=theta[76]  # txyeargrp3
		    txyear3_13Int=theta[77]
		    nonwhite_12Int=theta[78]  # nonwhite
		    nonwhite_13Int =theta[79]
		    nonwhite_23Int =theta[80]   
		    glandvol_12Int=theta[81]  # glandvol
		    glandvol_13Int =theta[82]
		    glandvol_23Int =theta[83]   
    	
    	
	        exp12<-exp(s2_12*stage[2]+s3_12*stage[3]+
	    				c2_12*comorbid[2]+c3_12*comorbid[3]+c4_12*comorbid[4]+
	    				g2_12*gleason[2]+g3_12*gleason[3]+g4_12*gleason[4]+g5_12*gleason[5]+
	    				tx2_12*tx[2]+
	    				pni2_12*pni[2]+
	    				age_12*age+
	    				logpsa_12*logbasepsa+
	    				txyear2_12*txyeargrp[2]+txyear3_12*txyeargrp[3]+
	    				nonwhite_12*nonwhite + glandvol_12*glandvol+
	    				
	    				stage_12Int*(stage[2]+ 2*stage[3])*tx[2]+
	    				comorbid_12Int*(comorbid[1] + 2*comorbid[3] + 3*comorbid[4])*tx[2]+
	    				gleason_12Int*(gleason[3] + 2*gleason[1] + 3*gleason[4] + 4*gleason[5])*tx[2]+
	    				pni2_12Int*pni[2]*tx[2]+
	    				age_12Int*age*tx[2]+
	    				logpsa_12Int*logbasepsa*tx[2]+
	    				txyear2_12Int*txyeargrp[1]*tx[2]+txyear3_12Int*txyeargrp[3]*tx[2]+
	    				nonwhite_12Int*nonwhite*tx[2]+ glandvol_12Int*glandvol*tx[2])
	    		exp13<-exp(s2_13*stage[2]+s3_13*stage[3]+
	    				c2_13*comorbid[2]+c3_13*comorbid[3]+c4_13*comorbid[4]+
	    				g2_13*gleason[2]+g3_13*gleason[3]+g4_13*gleason[4]+g5_13*gleason[5]+
	    				tx2_13*tx[2]+
	    				pni2_13*pni[2]+
	    				age_13*age+
	    				logpsa_13*logbasepsa+
	    				txyear2_13*txyeargrp[2]+txyear3_13*txyeargrp[3]+
	    				nonwhite_13*nonwhite + glandvol_13*glandvol+
	    				
	    				stage_13Int*(stage[2]+ 2*stage[3])*tx[2]+
	    				comorbid_13Int*(comorbid[1] + 2*comorbid[3]+ 3*comorbid[4])*tx[2]+
	    				gleason_13Int*(gleason[3] + 2*gleason[1] + 3*gleason[4] + 4*gleason[5])*tx[2]+
	    				pni2_13Int*pni[2]*tx[2]+
	    				age_13Int*age*tx[2]+
	    				logpsa_13Int*logbasepsa*tx[2]+
	    				txyear2_13Int*txyeargrp[1]*tx[2]+txyear3_13Int*txyeargrp[3]*tx[2]+
	    				nonwhite_13Int*nonwhite*tx[2]+ glandvol_13Int*glandvol*tx[2]) 
    							
    				
			base_haz = stepfun(x=cuts, y = c(b1, b2, b3, b4), right = F)
			base_Haz = function(t)
			{
				result = rep(NA,length(t))
				result = ifelse(t<=cuts[1], t*b1, result)
				result = ifelse(t>cuts[1] & t<=cuts[2], cuts[1]*b1 + (t-cuts[1])*b2, result)
				result = ifelse(t>cuts[2] & t<=cuts[3], cuts[1]*b1 + (cuts[2]-cuts[1])*b2 + (t-cuts[2])*b3, result)
				result = ifelse(t>cuts[3], cuts[1]*b1 + (cuts[2]-cuts[1])*b2 + (cuts[3]-cuts[2])*b3 + (t-cuts[3])*b4, result)
			}	
		
		    				
		    								

			if(input$CurTime == 0 & input$RecurEvent==FALSE){
				for(k in 1:length(times))
				{
					t = times[k]		
					Group1 = stats::integrate( integrate1_short, lower = 0,upper = t )$value
					Group2 = stats::integrate( integrate2_short, lower = 0,upper = t )$value 
					Group3 = exp(-base_Haz(t)*exp13)*exp(- (1/(l12^p12))*exp12*t^(p12))
					Group4 = stats::integrate( integrate4_short, lower = 0,upper = t )$value
					
		#Group1=cubature::adaptIntegrate(Vectorize(integrate1_short), lowerLimit = 0, upperLimit = t, maxEval = 15)$integral  
		#Group2=cubature::adaptIntegrate(Vectorize(integrate2_short), lowerLimit = 0, upperLimit = t, maxEval = 15)$integral  
		#Group4=cubature::adaptIntegrate(Vectorize(integrate4_short), lowerLimit = 0, upperLimit = t, maxEval = 15)$integral  
					
					
					Prob1Save[k] = Group1
					Prob2Save[k] = Group2
					Prob3Save[k] = Group3
					Prob4Save[k] = Group4
				}

			}else{
				if(input$RecurEvent == FALSE){ #no observed recurrence
					for(k in which(times>=input$CurTime))
					{
						t = times[k]		
						Group1 = stats::integrate( integrate1_short, lower = 0,upper = t  )$value
						Group2 = stats::integrate( integrate2_short, lower = 0,upper = t )$value 
						Group3 = exp(-base_Haz(t)*exp13)*exp(- (1/(l12^p12))*exp12*t^(p12))
						Group4 = stats::integrate( integrate4_short, lower = 0,upper = t  )$value
						Group1STAR = stats::integrate( integrate1_short, lower = 0,upper = input$CurTime  )$value
						Group2STAR = stats::integrate( integrate2_short, lower = 0,upper = input$CurTime )$value 
						Group3STAR = exp(-base_Haz(input$CurTime)*exp13)*exp(- (1/(l12^p12))*exp12*input$CurTime^(p12))
						Group4STAR = stats::integrate( integrate4_short, lower = 0,upper = input$CurTime  )$value
						Prob1Save[k] = ((Group1-Group1STAR)/Group3STAR)
						Prob2Save[k] = ((Group2-Group2STAR)/Group3STAR)
						Prob3Save[k] = (Group3/Group3STAR)
						Prob4Save[k] = ((Group4-Group4STAR)/Group3STAR)
					}
				}else{ #observed recurrence
					exp23<-exp(s3_23*(stage[2]+stage[3])+
		    				c2_23*comorbid[2]+c3_23*comorbid[3]+c4_23*comorbid[4]+
		    				g4_23*(gleason[2]+gleason[3]+gleason[4])+g5_23*gleason[5]+
		    				tx2_23*tx[2]+
		    				pni2_23*pni[2]+
		    				age_23*age+
		    				logpsa_23*logbasepsa+
		    				nonwhite_23*nonwhite + glandvol_23*glandvol+
		    				btr23*TRANS_TR(input$RecurTime)+
		   
		    				stage_23Int*(stage[2]+ stage[3])*tx[2]+ #group stages 2 and 3 here
		    				comorbid_23Int*(comorbid[2] + 2*comorbid[3] + 3*comorbid[4])*tx[2]+
		    				gleason_23Int*(gleason[3] + 2*gleason[1] + 3*gleason[4] + 4*gleason[5])*tx[2]+ 
		    				pni2_23Int*pni[2]*tx[2]+
		    				age_23Int*age*tx[2]+
		    				logpsa_23Int*logbasepsa*tx[2]+
		    				nonwhite_23Int*nonwhite*tx[2] + glandvol_23Int*glandvol*tx[2])  
					for(k in which(times>=input$CurTime))
					{
						t = times[k]		
						Prob1Save[k] =  1-(exp(-(1/(l23^p23))*exp23*((t-input$RecurTime)^p23))/ 
											exp(-(1/(l23^p23))*exp23*((input$CurTime-input$RecurTime)^p23))	)
						Prob2Save[k] =  exp(-(1/(l23^p23))*exp23*((t-input$RecurTime)^p23))/ 
											exp(-(1/(l23^p23))*exp23*((input$CurTime-input$RecurTime)^p23))	
						Prob3Save[k] = 0
						Prob4Save[k] = 0
					}					
				}
										
			}
		#2, recurred and alive; 1, recurred and died; 4, died without recurrence; 3, alive without recurrence

  			cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
			if(input$plottype=='S'){
				SUBSET = which(times>=input$CurTime)
				#par(mar = c(5,5,5,6), xpd = TRUE)
				plot(c(),c(), main = paste0('State Occupancy Probabilities'), ylim = c(0,1), xlim = c(0,max(times)), 
					xlab = 'Years from Baseline', 
					ylab = 'Proportion of Subjects in Group', cex.main = 1.5, cex.lab = 1.2)
				
				lower = rep(0,length(SUBSET))
				upper = Prob1Save[SUBSET]
				lines(times[SUBSET], upper  , lwd = 3)		
				if(input$BW == TRUE){
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), density = 25, angle = 0)		
				}else{
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[5])
				}
				lower = upper
				upper = lower + Prob4Save[SUBSET]
				lines(times[SUBSET], upper  , lwd = 3)
				if(input$BW == TRUE){
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = 'gray')
				}else{
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[4])
				}
				lower = upper
				upper = lower + Prob2Save[SUBSET]
				lines(times[SUBSET], upper  , lwd = 3)
				if(input$BW == TRUE){
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), density = 25, angle = 45)
				}else{
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[8])
				}
				lower = upper
				#upper = lower + Prob3Save[SUBSET]
				upper = rep(1,length(SUBSET))
				lines(times[SUBSET], upper  , lwd = 3)
				if(input$BW == TRUE){
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = 'white')
				}else{
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[3])
				}
				if(input$CurTime>0){abline(v=min(times[SUBSET]), lwd = 3)}
				if(input$plotNums == FALSE){
					if(input$BW == FALSE){
						legend(x='topleft', legend = c('Died with prior CF', 'Died without prior CF', 'Alive with prior CF', 'Alive without prior CF'), 
						fill = cbPalette[c(5,4,8,3)], cex = 1.2)
					}else{
						legend(x='topleft', legend = c('Died with prior CF', 'Died without prior CF', 'Alive with prior CF', 'Alive without prior CF'), 
						fill = c('black', 'gray', 'black', 'white'), density = c(25,NA,25,NA), angle = c(0,NA,45,NA),cex = 1.2)
					}
				}else{
					LAST = length(times)
					VALS = formatC( round(c(Prob1Save[LAST],Prob4Save[LAST],Prob2Save[LAST], Prob3Save[LAST]),3), format='f', digits=3 )
					if(input$BW == FALSE){					
						legend(x='topleft', legend = c(paste('Died with prior CF,   Prop. = ',VALS[1]), 
						paste('Died without prior CF,   Prop. = ',VALS[2]), 
						paste('Alive with prior CF,   Prop. = ',VALS[3]), paste( 'Alive without prior CF,   Prop. = ', VALS[4])),
						fill = cbPalette[c(5,4,8,3)], cex = 1.2)
					}else{
						legend(x='topleft', legend = c(paste('Died with prior CF,   Prop. = ',VALS[1]), 
						paste('Died without prior CF,   Prop. = ',VALS[2]), 
						paste('Alive with prior CF,   Prop. = ',VALS[3]), paste( 'Alive without prior CF,   Prop. = ', VALS[4])),
						fill = c('black', 'gray', 'black', 'white'), density = c(25,NA,25,NA), angle = c(0,NA,45,NA), cex = 1.2)						
					}
					points(rep(input$maxtime,4), y = cumsum(VALS), pch = 16, col = 'black') 
				}
				axis(side = 4)
			}else if(input$plottype == 'OS'){
				SUBSET = which(times>=input$CurTime)
				#par(mar = c(5,5,5,6), xpd = TRUE)
				plot(c(),c(), main = paste0('Overall Survival Probability'), ylim = c(0,1), xlim = c(0,max(times)), 
					xlab = 'Years from Baseline', 
					ylab = 'Survival Probability', cex.main = 1.5, cex.lab = 1.2)
				
				lower = rep(0,length(SUBSET))
				upper = Prob2Save[SUBSET] + Prob3Save[SUBSET]
				lines(times[SUBSET], upper  , lwd = 3)
				if(input$BW == FALSE){
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[5])		
				}else{
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = 'gray')
				}
				lower = upper
				#upper = lower + Prob1Save[SUBSET] + Prob4Save[SUBSET]
				upper = rep(1,length(SUBSET))
				lines(times[SUBSET], upper  , lwd = 3)
				if(input$BW == FALSE){
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[3])
				}else{
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = 'white')
				}
				if(input$CurTime>0){abline(v=min(times[SUBSET]), lwd = 3)}
	
				legend(x='bottomleft', legend = c('Alive', 'Died'), 
						fill = cbPalette[c(5,3)])#,bty = 'n')
						
				if(input$plotNums == FALSE){
					if(input$BW == FALSE){
					legend(x='bottomleft', legend = c('Alive', 'Died'), 
						fill = cbPalette[c(5,3)], cex = 1.2)#,bty = 'n')
					}else{
						legend(x='bottomleft', legend = c('Alive', 'Died'), 
						fill = c('gray', 'white'), cex = 1.2)#,bty = 'n')
					}
				}else{
					LAST = length(times)
					VALS = formatC(round(c(Prob2Save[LAST]+Prob3Save[LAST], Prob1Save[LAST]+Prob4Save[LAST]),3), format='f', digits=3 )
					if(input$BW == FALSE){
					legend(x='bottomleft', legend = c(paste('Alive,   Prop. = ',VALS[1]), paste('Died,   Prop. = ', VALS[2])), 
						fill = cbPalette[c(5,3)], cex = 1.2)
					}else{
						legend(x='bottomleft', legend = c(paste('Alive,   Prop. = ',VALS[1]), paste('Died,   Prop. = ', VALS[2])), 
						fill = c('gray', 'white'), cex = 1.2)
					}
					points(input$maxtime, y = VALS[1], pch = 16, col = 'black') 	
				}
				axis(side = 4)
			}else if(input$plottype == 'EF'){
				SUBSET = which(times>=input$CurTime)
				#par(mar = c(5,5,5,6), xpd = TRUE)
				plot(c(),c(), main = paste0('Metastatic-Free Survival Probability'), ylim = c(0,1), xlim = c(0,max(times)), 
					xlab = 'Years from Baseline', 
					ylab = 'Metastatic-Free Survival Probability', cex.main = 1.5, cex.lab = 1.2)
				
				lower = rep(0,length(SUBSET))
				upper = Prob3Save[SUBSET]
				lines(times[SUBSET], upper  , lwd = 3)
				if(input$BW == FALSE){
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[5])		
				}else{
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = 'gray')
				}
				lower = upper
				#upper = lower + Prob1Save[SUBSET] + Prob4Save[SUBSET] + Prob2Save[SUBSET] 
				upper = rep(1,length(SUBSET))
				lines(times[SUBSET], upper  , lwd = 3)
				if(input$BW == FALSE){
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = cbPalette[3])
				}else{
					polygon(c(times[SUBSET],rev(times[SUBSET])), c(upper, rev(lower)), col = 'white')
				}
				if(input$CurTime>0){abline(v=min(times[SUBSET]), lwd = 3)}

				if(input$plotNums == FALSE){
					if(input$BW == FALSE){
					legend(x='bottomleft', legend = c('Alive without prior CF', 'Died and/or CF'), 
						fill = cbPalette[c(5,3)], cex = 1.2)#,bty = 'n')
					}else{
						legend(x='bottomleft', legend = c('Alive without prior CF', 'Died and/or CF'), 
						fill =  c('gray', 'white'), cex = 1.2)#,bty = 'n')
					}
				}else{
					LAST = length(times)
					VALS = formatC(round(c(Prob3Save[LAST], Prob1Save[LAST]+Prob4Save[LAST]+Prob2Save[LAST]),3), format='f', digits=3 )
					if(input$BW == FALSE){
					legend(x='bottomleft', legend = c(paste('Alive without prior CF,   Prop. = ',VALS[1]), 
						paste('Died and/or CF,   Prop. = ', VALS[2])), 
						fill = cbPalette[c(5,3)], cex = 1.2)#,bty = 'n')						
					}else{
					legend(x='bottomleft', legend = c(paste('Alive without prior CF,   Prop. = ',VALS[1]), 
						paste('Died and/or CF,   Prop. = ', VALS[2])), 
						fill = c('gray', 'white'), cex = 1.2)#,bty = 'n')						
					}
					points(input$maxtime, y = VALS[1], pch = 16, col = 'black') 
				}		
				axis(side = 4)	
			}else{
				SUBSET = which(times>=input$CurTime)
				
				x = data.frame(x1 = Prob2Save[times == input$maxtime], x2 =(Prob1Save[times == input$maxtime] + Prob4Save[times == input$maxtime]), 
					x3 = Prob3Save[times == input$maxtime])
				
				p = ggtern::ggtern(data=x,ggtern::aes(x=x1,y=x2,z=x3)) + 
			  ggplot2::geom_point(fill="black",shape=21,size=4) + 
			ggtern::theme_custom(base_size = 16, base_family = "",
			  tern.panel.background = 'white',
			  col.T = ifelse(input$BW == FALSE,"blue", 'black'), col.L = ifelse(input$BW == FALSE,"darkgreen", 'black'), col.R = ifelse(input$BW == FALSE,"purple", 'black'))+
			  ggplot2::labs(title="Ternary Plot for State Occupancy Probabilities")+
			 ggtern:: theme_showarrows()+  
			  ggtern::Llab(label = '',labelarrow = 'Alive with prior CF') + 
			  ggtern::Rlab(label = '', labelarrow = 'Alive without prior CF') + 
			  ggtern::Tlab(label = '', labelarrow = 'Died')+
			  ggplot2::theme(axis.text=ggplot2::element_text(size=16),
			        axis.title=ggplot2::element_text(size=18,face="bold"))
			  	print(p)	
			}		
			# if(input$SavePlot == TRUE){
				# grDevices::dev.print(pdf, paste(plotname, '.pdf', sep = ''))
			# }					
		}, height = 420, width = 420)#end of renderplot
	})#end of local			
 	output$plots <- shiny::renderUI({print(plotInput())})
}#end of server

