library(INLA)
library(rgdal)
library(tidyverse)
library(shiny)

fillmap2<-function(map, figtitle, y , leg.loc="beside", y.scl=NULL,
                   main.cex=1.5,main.line=0,map.lty=1,leg.rnd=0,
                   leg.cex=1){
  
  # 0: dark 1: light light Current shading ranges from darkest to light gray white (to distinguish with lakes).
  y.uq=sort(unique(c(y,y.scl)))
  cols<-viridis(length(y.uq),direction=-1)
  shading=y
  for (i in 1:length(y)){
    shading[i]<-cols[which(y.uq==y[i])]
  }
  
  par(mar=c(0,0,2,0))
  if (leg.loc=="beside"){
    layout(matrix(1:2,ncol=2),width=c(.8,.2))
  } else 
    if (leg.loc=="below"){
      layout(matrix(1:2,nrow=2),height=c(.6,.4))
    } else (print("leg.loc options are below or beside"))
  
  plot(map,col=shading,axes=F, lty=map.lty)
  title(main=figtitle,cex.main=main.cex,line=main.line) 
  
  par(mar=c(5, 4, 4, 2) + 0.1)
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
  cols.5=cols[seq(1,length(y.uq),length.out=5)]
  lab.5=cols.5
  for (i in 1:5){lab.5[i]=y.uq[which(cols==cols.5[i])[1]]}
  lab.5=round(as.numeric(lab.5),leg.rnd)
  par(mar=c(0,0,0,0))
  if (leg.loc=="beside"){
    legend_image <- as.raster(matrix(cols, ncol=1))
    text(x=1.6, 
         y = seq(0,length(y.uq),length.out=5)/length(y.uq),
         labels = rev(lab.5), cex=leg.cex)
    rasterImage(legend_image, 0, 0, 1,1)
  } else{
    legend_image <- as.raster(matrix(cols, nrow=1))
    text(y=-0.25, 
         x = seq(0,length(y.uq),length.out=5)/(length(y.uq)*.5),
         labels = lab.5, cex=leg.cex)
    rasterImage(legend_image, 0, 0, 2,1)
  }
}

NCtracts=readOGR("/Users/ryanskeete/UNCW/MSC DATA SCIENCE/DSC551 Spatial Temporal Analysis/CaseStudy4/tl_2016_37_tract/tl_2016_37_tract.shp")
dim(NCtracts)
NHtracts<-NCtracts[which(NCtracts$COUNTYFP=="129"),]
NHtracts<-NHtracts[which(NHtracts$TRACTCE!=990100),]
# sort NHtracts by census tract in ascending order
NHtracts$TRACTCE<-as.numeric(NHtracts$TRACTCE)
NHtractsR<-NHtracts[order(NHtracts$TRACTCE),]

data<-read_csv("/Users/ryanskeete/UNCW/MSC DATA SCIENCE/DSC551 Spatial Temporal Analysis/CaseStudy4/wpd_arrests_race_tract1018_clean.csv", show_col_types = FALSE)
EfixedEffects_Tot=read.csv("/Users/ryanskeete/UNCW/MSC DATA SCIENCE/DSC551 Spatial Temporal Analysis/CaseStudy4/ExpfixedEffects_Tot.csv")[-1]
EfixedEffects_B=read.csv("/Users/ryanskeete/UNCW/MSC DATA SCIENCE/DSC551 Spatial Temporal Analysis/CaseStudy4/ExpfixedEffects_Black.csv")[-1]
EfixedEffects_W=read.csv("/Users/ryanskeete/UNCW/MSC DATA SCIENCE/DSC551 Spatial Temporal Analysis/CaseStudy4/ExpfixedEffects_White.csv")[-1]
data$census_tract_code<-as.numeric(data$census_tract_code)
dataR=data[which(data$census_tract_code!=990100),]

# Total, black, and white arrests as a percent of total arrests
dataR$B_W_arrests_pct <- (dataR$arrests_B+dataR$arrests_W+0.1)/(dataR$arrests_total+0.1)

#Black arrests as a precent of total arrests
dataR$arrests_B_pct_tot <- (dataR$arrests_B+0.1)/(dataR$arrests_total+0.1)
dataR$arrests_W_pct_tot <- (dataR$arrests_W+0.1)/(dataR$arrests_total+0.1)

# Arrests as a percent of the appropriate population
dataR$Tot_arrests_pctTOTal_POP <-dataR$arrests_total/(dataR$ct_pop+0.1)
dataR$B_arrests_pctB_POP<-dataR$arrests_B/(dataR$ct_black+0.1)
dataR$W_arrests_pctW_POP<-dataR$arrests_W/(dataR$ct_white+0.1)

# SIR of arrests
#dataR$ct_pop[dataR$ct_pop==0] <- 0.1
dataR$SIR <- dataR$arrests_total/(sum(dataR$arrests_total)/sum(dataR$ct_pop)*dataR$ct_pop+0.1)
dataR$SIR_B <- dataR$arrests_B/(sum(dataR$arrests_B)/sum(dataR$ct_black)*(dataR$ct_black)+0.1)
dataR$SIR_W <- dataR$arrests_W/(sum(dataR$arrests_W)/sum(dataR$ct_white)*(dataR$ct_white)+0.1)

d.inla=read.csv("/Users/ryanskeete/UNCW/MSC DATA SCIENCE/DSC551 Spatial Temporal Analysis/CaseStudy4/INLAdata.csv")[,-1]
#Expected total arrests
EcountTotArrests<-(sum(dataR$arrests_total)/sum(dataR$ct_pop)*dataR$ct_pop+0.1)
#Create dummy variable called GovernmentEntity which 1 when either police depart or jail is 1
d.inla$GovernmentEntity <- ifelse(d.inla$police.dept==1|d.inla$jails==1,1,0)

resp <- "arrests_total" # Response variable

covar <- c("black", # Percent black population 
           "poverty", # Percent of the population living in poverty  
           "educBachPlus", # Percent of the population with a BSc or more
           "male", # Percent male population
           "secperctot",# Second homes as a percent of all homes
           "age1824.perc", # Percent of the population aged 18-24
           "GovernmentEntity")
ftot2 <- as.formula(paste0(resp, " ~ ",
                           paste(covar, collapse = " + ")," +  f(id2,model = 'iid',param=c(2,1))"))
restot2  <- inla(ftot2,data = d.inla,control.compute = list(dic=TRUE,waic=TRUE), verbose=FALSE, family="Poisson", E= EcountTotArrests )

#Expected total White arrests
EcountWArrests<-(sum(dataR$arrests_W)/sum(dataR$ct_white)*dataR$ct_white+0.1)
respW <- "arrests_W" # Response variable
ftotW2 <- as.formula(paste0(resp, " ~ ",
                            paste(covar, collapse = " + ")," +  f(id2,model = 'iid',param=c(2,1))"))
restotW2  <- inla(ftotW2,data = d.inla,control.compute = list(dic=TRUE,waic=TRUE), verbose=FALSE, family="Poisson", E= EcountWArrests )

#Expected total arrests
EcountBArrests<-(sum(dataR$arrests_B)/sum(dataR$ct_black)*dataR$ct_black+0.1)
respB <- "arrests_B" # Response variable
ftotB2 <- as.formula(paste0(resp, " ~ ",
                            paste(covar, collapse = " + ")," +  f(id2,model = 'iid',param=c(2,1))"))
restotB2  <- inla(ftotB2,data = d.inla,control.compute = list(dic=TRUE,waic=TRUE), verbose=FALSE, family="Poisson", E= EcountBArrests ) 


ui <- fluidPage(
  titlePanel("2010-2018 WPD Arrest Data"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = "year", 
                  label = "Choose a year", 
                  value=2011, min=2010, max= 2018, sep="",     animate=animationOptions(interval=1000,loop=TRUE)),
      radioButtons(inputId = "Data", label="Choose Data", choices=c("Total      Arrests"="TotArrests","White Arrests"="WArrests","Black Arrests"="BArrests")),
      radioButtons(inputId = "Analysis", label="Choose Analysis", choices=c("None"="none","Proportion of Population"="PerPop","Proportion of Total Arrests"="PerTotArrest", "Standardize Incidence Ratio"="SIR","Poisson Regression"="PReg"))),
    mainPanel(textOutput("text"),
              plotOutput("map"),
              tableOutput("table")
    )
  ))



server <- function(input, output) {
  output$map <- renderPlot({if (input$year==2010 & input$Data=="TotArrests"&input$Analysis=="none"){fillmap2(NHtractsR,"2010 Total Arrests",dataR$arrests_total[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_total,leg.loc="below")}
    if (input$year==2011 & input$Data=="TotArrests"&input$Analysis=="none"){fillmap2(NHtractsR,"2011 Total Arrests",dataR$arrests_total[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_total,leg.loc="below")}
    if (input$year==2012 & input$Data=="TotArrests"&input$Analysis=="none"){fillmap2(NHtractsR," 2012 Total Arrests",dataR$arrests_total[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_total,leg.loc="below")}
    if (input$year==2013 & input$Data=="TotArrests"&input$Analysis=="none"){fillmap2(NHtractsR,"2013 Total Arrests",dataR$arrests_total[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_total,leg.loc="below")}
    if (input$year==2014 & input$Data=="TotArrests"&input$Analysis=="none"){fillmap2(NHtractsR,"2014 Total Arrests",dataR$arrests_total[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_total,leg.loc="below")}
    if (input$year==2015 & input$Data=="TotArrests"&input$Analysis=="none"){fillmap2(NHtractsR,"2015 Total Arrests",dataR$arrests_total[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_total,leg.loc="below")}
    if (input$year==2016 & input$Data=="TotArrests"&input$Analysis=="none"){fillmap2(NHtractsR,"2016 Total Arrests",dataR$arrests_total[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_total,leg.loc="below")}
    if (input$year==2017 & input$Data=="TotArrests"&input$Analysis=="none"){fillmap2(NHtractsR,"2017 Total Arrests",dataR$arrests_total[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_total,leg.loc="below")}
    if (input$year==2018 & input$Data=="TotArrests"&input$Analysis=="none"){fillmap2(NHtractsR,"2018 Total Arrests",dataR$arrests_total[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_total,leg.loc="below")}
    if (input$year==2010 & input$Data=="WArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"White Arrests", sep=" "),dataR$arrests_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W,leg.loc="below")}
    if (input$year==2011 & input$Data=="WArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"White Arrests", sep=" "),dataR$arrests_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W,leg.loc="below")}
    if (input$year==2012 & input$Data=="WArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"White Arrests", sep=" "),dataR$arrests_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W,leg.loc="below")}
    if (input$year==2013 & input$Data=="WArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"White Arrests", sep=" "),dataR$arrests_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W,leg.loc="below")}
    if (input$year==2014 & input$Data=="WArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"White Arrests", sep=" "),dataR$arrests_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W,leg.loc="below")}
    if (input$year==2015 & input$Data=="WArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"White Arrests", sep=" "),dataR$arrests_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W,leg.loc="below")}
    if (input$year==2016 & input$Data=="WArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"White Arrests", sep=" "),dataR$arrests_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W,leg.loc="below")}
    if (input$year==2017 & input$Data=="WArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"White Arrests", sep=" "),dataR$arrests_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W,leg.loc="below")}
    if (input$year==2018 & input$Data=="WArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"White Arrests", sep=" "),dataR$arrests_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W,leg.loc="below")}
    if (input$year==2010 & input$Data=="BArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"Black Arrests", sep=" "),dataR$arrests_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B,leg.loc="below")}
    if (input$year==2011 & input$Data=="BArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"Black Arrests", sep=" "),dataR$arrests_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B,leg.loc="below")}
    if (input$year==2012 & input$Data=="BArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"Black Arrests", sep=" "),dataR$arrests_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B,leg.loc="below")}
    if (input$year==2013 & input$Data=="BArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"Black Arrests", sep=" "),dataR$arrests_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B,leg.loc="below")} 
    if (input$year==2014 & input$Data=="BArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"Black Arrests", sep=" "),dataR$arrests_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B,leg.loc="below")}
    if (input$year==2015 & input$Data=="BArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"Black Arrests", sep=" "),dataR$arrests_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B,leg.loc="below")}
    if (input$year==2016 & input$Data=="BArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"Black Arrests", sep=" "),dataR$arrests_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B,leg.loc="below")}
    if (input$year==2017 & input$Data=="BArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"Black Arrests", sep=" "),dataR$arrests_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B,leg.loc="below")} 
    if (input$year==2018 & input$Data=="BArrests"&input$Analysis=="none"){fillmap2(NHtractsR, paste(input$year,"Black Arrests", sep=" "),dataR$arrests_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B,leg.loc="below")} 
    if (input$year==2010 & input$Data=="TotArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black & White Arrests as Percentage of Total Arrests", sep=" "),dataR$B_W_arrests_pct[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_W_arrests_pct,leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="TotArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black & White Arrests as Percentage of Total Arrests", sep=" "),dataR$B_W_arrests_pct[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_W_arrests_pct,leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="TotArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black & White Arrests as Percentage of Total Arrests", sep=" "),dataR$B_W_arrests_pct[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_W_arrests_pct,leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="TotArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black & White Arrests as Percentage of Total Arrests", sep=" "),dataR$B_W_arrests_pct[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_W_arrests_pct,leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="TotArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black & White Arrests as Percentage of Total Arrests", sep=" "),dataR$B_W_arrests_pct[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_W_arrests_pct,leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="TotArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black & White Arrests as Percentage of Total Arrests", sep=" "),dataR$B_W_arrests_pct[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_W_arrests_pct,leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="TotArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black & White Arrests as Percentage of Total Arrests", sep=" "),dataR$B_W_arrests_pct[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_W_arrests_pct,leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="TotArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black & White Arrests as Percentage of Total Arrests", sep=" "),dataR$B_W_arrests_pct[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_W_arrests_pct,leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="TotArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black & White Arrests as Percentage of Total Arrests", sep=" "),dataR$B_W_arrests_pct[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_W_arrests_pct,leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="WArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_W_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="WArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_W_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="WArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_W_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="WArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_W_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="WArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_W_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W_pct_tot,leg.loc="below",leg.rnd=3)}  
    if (input$year==2015 & input$Data=="WArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_W_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="WArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_W_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="WArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_W_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="WArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_W_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_W_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="BArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_B_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="BArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_B_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="BArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_B_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="BArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_B_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="BArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_B_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="BArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_B_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="BArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_B_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="BArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_B_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="BArrests"&input$Analysis=="PerTotArrest"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Percentage of Total Arrests", sep=" "),dataR$arrests_B_pct_tot[dataR$year==input$year],leg.cex=1,y.scl=dataR$arrests_B_pct_tot,leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="TotArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Total Arrests as Proportion of Population", sep=" "),dataR$Tot_arrests_pctTOTal_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$Tot_arrests_pctTOTal_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="TotArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Total Arrests as Proportion of Population", sep=" "),dataR$Tot_arrests_pctTOTal_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$Tot_arrests_pctTOTal_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="TotArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Total Arrests as Proportion of Population", sep=" "),dataR$Tot_arrests_pctTOTal_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$Tot_arrests_pctTOTal_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="TotArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Total Arrests as Proportion of Population", sep=" "),dataR$Tot_arrests_pctTOTal_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$Tot_arrests_pctTOTal_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="TotArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Total Arrests as Proportion of Population", sep=" "),dataR$Tot_arrests_pctTOTal_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$Tot_arrests_pctTOTal_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="TotArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Total Arrests as Proportion of Population", sep=" "),dataR$Tot_arrests_pctTOTal_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$Tot_arrests_pctTOTal_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="TotArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Total Arrests as Proportion of Population", sep=" "),dataR$Tot_arrests_pctTOTal_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$Tot_arrests_pctTOTal_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="TotArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Total Arrests as Proportion of Population", sep=" "),dataR$Tot_arrests_pctTOTal_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$Tot_arrests_pctTOTal_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="TotArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Total Arrests as Proportion of Population", sep=" "),dataR$Tot_arrests_pctTOTal_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$Tot_arrests_pctTOTal_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="WArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Proportion of White Population", sep=" "),dataR$W_arrests_pctW_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$W_arrests_pctW_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="WArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Proportion of White Population", sep=" "),dataR$W_arrests_pctW_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$W_arrests_pctW_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="WArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Proportion of White Population", sep=" "),dataR$W_arrests_pctW_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$W_arrests_pctW_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="WArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Proportion of White Population", sep=" "),dataR$W_arrests_pctW_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$W_arrests_pctW_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="WArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Proportion of White Population", sep=" "),dataR$W_arrests_pctW_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$W_arrests_pctW_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="WArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Proportion of White Population", sep=" "),dataR$W_arrests_pctW_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$W_arrests_pctW_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="WArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Proportion of White Population", sep=" "),dataR$W_arrests_pctW_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$W_arrests_pctW_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="WArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Proportion of White Population", sep=" "),dataR$W_arrests_pctW_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$W_arrests_pctW_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="WArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"White Arrests as Proportion of White Population", sep=" "),dataR$W_arrests_pctW_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$W_arrests_pctW_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="BArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Proportion of Black Population", sep=" "),dataR$B_arrests_pctB_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_arrests_pctB_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="BArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Proportion of Black Population", sep=" "),dataR$B_arrests_pctB_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_arrests_pctB_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="BArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Proportion of Black Population", sep=" "),dataR$B_arrests_pctB_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_arrests_pctB_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="BArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Proportion of Black Population", sep=" "),dataR$B_arrests_pctB_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_arrests_pctB_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="BArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Proportion of Black Population", sep=" "),dataR$B_arrests_pctB_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_arrests_pctB_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="BArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Proportion of Black Population", sep=" "),dataR$B_arrests_pctB_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_arrests_pctB_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="BArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Proportion of Black Population", sep=" "),dataR$B_arrests_pctB_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_arrests_pctB_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="BArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Proportion of Black Population", sep=" "),dataR$B_arrests_pctB_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_arrests_pctB_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="BArrests"&input$Analysis=="PerPop"){fillmap2(NHtractsR, paste(input$year,"Black Arrests as Proportion of Black Population", sep=" "),dataR$B_arrests_pctB_POP[dataR$year==input$year],leg.cex=1,y.scl=dataR$B_arrests_pctB_POP,leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="TotArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Total Arrests", sep=" "),dataR$SIR[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR,leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="TotArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Total Arrests", sep=" "),dataR$SIR[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR,leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="TotArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Total Arrests", sep=" "),dataR$SIR[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR,leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="TotArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Total Arrests", sep=" "),dataR$SIR[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR,leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="TotArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Total Arrests", sep=" "),dataR$SIR[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR,leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="TotArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Total Arrests", sep=" "),dataR$SIR[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR,leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="TotArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Total Arrests", sep=" "),dataR$SIR[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR,leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="TotArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Total Arrests", sep=" "),dataR$SIR[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR,leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="TotArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Total Arrests", sep=" "),dataR$SIR[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR,leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="WArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of White Arrests", sep=" "),dataR$SIR_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_W,leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="WArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of White Arrests", sep=" "),dataR$SIR_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_W,leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="WArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of White Arrests", sep=" "),dataR$SIR_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_W,leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="WArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of White Arrests", sep=" "),dataR$SIR_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_W,leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="WArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of White Arrests", sep=" "),dataR$SIR_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_W,leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="WArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of White Arrests", sep=" "),dataR$SIR_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_W,leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="WArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of White Arrests", sep=" "),dataR$SIR_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_W,leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="WArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of White Arrests", sep=" "),dataR$SIR_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_W,leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="WArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of White Arrests", sep=" "),dataR$SIR_W[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_W,leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="BArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Black Arrests", sep=" "),dataR$SIR_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_B,leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="BArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Black Arrests", sep=" "),dataR$SIR_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_B,leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="BArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Black Arrests", sep=" "),dataR$SIR_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_B,leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="BArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Black Arrests", sep=" "),dataR$SIR_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_B,leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="BArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Black Arrests", sep=" "),dataR$SIR_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_B,leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="BArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Black Arrests", sep=" "),dataR$SIR_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_B,leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="BArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Black Arrests", sep=" "),dataR$SIR_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_B,leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="BArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Black Arrests", sep=" "),dataR$SIR_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_B,leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="BArrests"&input$Analysis=="SIR"){fillmap2(NHtractsR, paste(input$year,"Standardize Incidence Ratio of Black Arrests", sep=" "),dataR$SIR_B[dataR$year==input$year],leg.cex=1,y.scl=dataR$SIR_B,leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="TotArrests"&input$Analysis=="PReg"){i<-0
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Total Arrests", sep=" "),exp(restot2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restot2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="TotArrests"&input$Analysis=="PReg"){i<-1
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Total Arrests", sep=" "),exp(restot2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restot2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="TotArrests"&input$Analysis=="PReg"){i<-2
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Total Arrests", sep=" "),exp(restot2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restot2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="TotArrests"&input$Analysis=="PReg"){i<-3
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Total Arrests", sep=" "),exp(restot2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restot2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="TotArrests"&input$Analysis=="PReg"){i<-4
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Total Arrests", sep=" "),exp(restot2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restot2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="TotArrests"&input$Analysis=="PReg"){i<-5
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Total Arrests", sep=" "),exp(restot2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restot2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="TotArrests"&input$Analysis=="PReg"){i<-6
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Total Arrests", sep=" "),exp(restot2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restot2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="TotArrests"&input$Analysis=="PReg"){i<-7
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Total Arrests", sep=" "),exp(restot2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restot2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="TotArrests"&input$Analysis=="PReg"){i<-8
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Total Arrests", sep=" "),exp(restot2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restot2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="WArrests"&input$Analysis=="PReg"){i<-0
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of White Arrests", sep=" "),exp(restotW2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotW2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="WArrests"&input$Analysis=="PReg"){i<-1
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of White Arrests", sep=" "),exp(restotW2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotW2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="WArrests"&input$Analysis=="PReg"){i<-2
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of White Arrests", sep=" "),exp(restotW2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotW2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="WArrests"&input$Analysis=="PReg"){i<-3
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of White Arrests", sep=" "),exp(restotW2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotW2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="WArrests"&input$Analysis=="PReg"){i<-4
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of White Arrests", sep=" "),exp(restotW2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotW2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="WArrests"&input$Analysis=="PReg"){i<-5
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of White Arrests", sep=" "),exp(restotW2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotW2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="WArrests"&input$Analysis=="PReg"){i<-6
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of White Arrests", sep=" "),exp(restotW2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotW2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="WArrests"&input$Analysis=="PReg"){i<-7
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of White Arrests", sep=" "),exp(restotW2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotW2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="WArrests"&input$Analysis=="PReg"){i<-8
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of White Arrests", sep=" "),exp(restotW2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotW2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2010 & input$Data=="BArrests"&input$Analysis=="PReg"){i<-0
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Black Arrests", sep=" "),exp(restotB2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotB2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2011 & input$Data=="BArrests"&input$Analysis=="PReg"){i<-1
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Black Arrests", sep=" "),exp(restotB2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotB2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2012 & input$Data=="BArrests"&input$Analysis=="PReg"){i<-2
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Black Arrests", sep=" "),exp(restotB2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotB2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2013 & input$Data=="BArrests"&input$Analysis=="PReg"){i<-3
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Black Arrests", sep=" "),exp(restotB2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotB2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2014 & input$Data=="BArrests"&input$Analysis=="PReg"){i<-4
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Black Arrests", sep=" "),exp(restotB2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotB2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2015 & input$Data=="BArrests"&input$Analysis=="PReg"){i<-5
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Black Arrests", sep=" "),exp(restotB2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotB2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2016 & input$Data=="BArrests"&input$Analysis=="PReg"){i<-6
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Black Arrests", sep=" "),exp(restotB2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotB2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2017 & input$Data=="BArrests"&input$Analysis=="PReg"){i<-7
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Black Arrests", sep=" "),exp(restotB2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotB2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)}
    if (input$year==2018 & input$Data=="BArrests"&input$Analysis=="PReg"){i<-8
    fillmap2(NHtractsR, paste(input$year,"Random Effects Model of Black Arrests", sep=" "),exp(restotB2$summary.random$id2$mean[(1+44*i):(44+44*i)]),leg.cex=1,y.scl=exp(restotB2$summary.random$id2$mean),leg.loc="below",leg.rnd=3)} 
  })
  
  output$text <- renderText({if (input$Analysis=="none") {paste("No adjustments specified. These are counts of arrests for the selected data.")} else
    if (input$Analysis=="PerPop"){paste("The 'Proportion of the Population' adjustment displays the selected arrests counts divided by the appropriate population (e.g. Black population only when 'Black Only Arrests' is selected).")} else
      if (input$Analysis=="PerTotArrest"){paste("The 'Proportion of Total Arrests' adjustment displays the selected arrests counts divided by the total arrest counts .")} else
        if (input$Analysis=="SIR"){paste("The standardized incidence ratio (SIR) adjustment was applied here. SIR is a method of adjusting for tract population by calculating a ratio of the observed count of arrests to the expected counts of arrests. The expected count of arrests is calculated as a rate of arrests over all tracts and years times the tract population for a given year. The population used reflects the data being displayed (e.g. Black population only when 'Black Only Arrests' is selected). Values greater than 1 suggest more observed arrests than expected.")} else
          if (input$Analysis=="PReg"){paste("The Poisson regression option for adjustment was applied here. In this method of adjustment, a Poisson regression model with spatio-temporal covariate adjustment was applied (see table output). The mapped values display the residual spatial variation in arrests after adjustment where higher (darker) values indicate areas of increased risk. All estimates are transformed so that they can be interpreted as a multiplicative change in the relative risk of arrests. Tract population is indirectly adjusted for.")}
  })
  
  output$table <- renderTable({
    if (input$Analysis=="PReg"&input$Data=="TotArrests"){
      
      tabl=EfixedEffects_Tot[2:8,c(1,3,5)]
      colnames(tabl)<-c("Mean","95% CI Lower Bound","95% CI Upper Bound")
      rownames(tabl)<-c("% Black",
                        "% Living in Poverty",
                        "% Bachelors degree or more",
                        "% Male",
                        "% Secondary Homes",
                        "% Aged 18-24",
                        "Government Entity")
      tabl}
    else
      if (input$Analysis=="PReg"&input$Data=="WArrests"){
        tabl=EfixedEffects_W[2:8,c(1,3,5)]
        colnames(tabl)<-c("Mean","95% CI Lower Bound","95% CI Upper Bound")
        rownames(tabl)<-c("% Black",
                          "% Living in Poverty",
                          "% Bachelors degree or more",
                          "% Male",
                          "% Secondary Homes",
                          "% Aged 18-24",
                          "Government Entity")
        tabl}
    else
      if (input$Analysis=="PReg"&input$Data=="BArrests"){
        tabl=EfixedEffects_B[2:8,c(1,3,5)]
        colnames(tabl)<-c("Mean","95% CI Lower Bound","95% CI Upper Bound")
        rownames(tabl)<-c("% Black",
                          "% Living in Poverty",
                          "% Bachelors degree or more",
                          "% Male",
                          "% Secondary Homes",
                          "% Aged 18-24",
                          "Government Entity")
        tabl}},
    
    rownames=T,colnames=T,digits=3,width="auto")
  
}
shinyApp(ui = ui, server = server)