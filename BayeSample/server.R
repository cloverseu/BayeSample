library(shiny)
library(shinydashboard)
library(ggplot2)
library(DT)
library(ggvis)
library(shinyBS)
library(plotly)
library(rsconnect)
source("D:\\软件\\大论文\\sampleApp\\BayeSample\\classfication\\normal.R", encoding="utf-8")
source("D:\\软件\\大论文\\sampleApp\\BayeSample\\classfication\\proportion.R", encoding="utf-8")
source("D:\\软件\\大论文\\sampleApp\\BayeSample\\classfication\\assurance.R", encoding="utf-8")
source("D:\\软件\\大论文\\sampleApp\\BayeSample\\classfication\\classicalSize.R", encoding="utf-8")
source("D:\\软件\\大论文\\sampleApp\\BayeSample\\classfication\\bayesFactorSSD.R", encoding="utf-8")
source("D:\\软件\\大论文\\sampleApp\\BayeSample\\classfication\\ppSSR.R", encoding="utf-8")
#library(SampleSizeMeans)

wwwPath <- paste0(getwd(),"/www")
addResourcePath("localfiles", wwwPath)

# rsconnect::setAccountInfo(name='bayessd',
#                           token='5041A898C0793B05864F6CF83A38BC47',
#                           secret='ePq6SYrut9ILWcR8YpN5vze2wqY2qh/JNaTcw4oQ')
# rsconnect::deployApp("D:\\软件\\大论文\\sampleApp\\BayeSample")



rsconnect::deployApp("D:\\软件\\大论文\\sampleApp\\BayeSample")

shinyServer(function(input, output, session) {
  
  
  rv <- reactiveValues(stats = "mu.varunknow.seq",sidebarMenu="dashboard",datatype="normal")
  

  observe({
    if(input$trial == "two"){
      rv$stats  <- input$stats3
    }else{
      if(input$datatype == "binary"){
        rv$stats <- input$stats1
      }else if(input$datatype == "survival"){
        rv$stats <- input$stats2
      }else{
        rv$stats <- input$stats
      }
    }
    print("999")
    print(rv$stats)
  })
  
  observeEvent(input$stats, {
    # changes in this variable will make the UI "controls" to be rendered again
    #rv$stats <- input$stats
    print("888")
    print(rv$stats)
  })
  
  observeEvent(input$datatype, {
    # changes in this variable will make the UI "controls" to be rendered again
    rv$datatype <- input$datatype
    print("666")
    print(input$stats)
    print(rv$stats)
  })
  #output$menu <- renderMenu({
   # sidebarMenu(
    #  menuItem("Menu item", icon = icon("calendar"))
    #)
  #})
  
 
  observe({
    print("test")
    print(input)
    print(input$sbMenu)
  })
  
  #右侧参数控制面板
  output$param <- renderUI({
    print("777")
    print(rv$sidebarMenu)
    func <- switch( rv$stats,
        "mu.varunknow.seq" = mu.varunknow.seq,
        "mu.varknown.seq" = mu.varknown.seq,
        "mudiff.equalvar.seq" = mudiff.equalvar.seq,
        "mudiff.unequalvar.seq" = mudiff.unequalvar.seq,
        "propdiff.seq" = propdiff.seq,
        #"propmbl.seq" = propmbl.seq,
        "inv.ANDks" = inv.ANDks,
        "inv.ANDus" = inv.ANDus,
        "inv.AWBD" = inv.AWBD,
        "inv.ASDExp" = inv.ASDExp,
        "BFt" = BFt,
        "BFpp" = BFpp,
        "BFsur" = BFsur,
        "ppSSR" = ppSSR,
        "PSSD" = PSSD
        
      
    )
    funclist = formals(func)
    #list中的CNname为对应的中文名称
    #列出所有可能参数的中文名称或者在fun中定义
    #列出所有的参数
    paramlist <- names(funclist)
    CNname <- funclist$CNname
    Textparam <- funclist$Textparam
    tips <- funclist$tips
    defaultValue <- funclist$defaultValue
    ##print("555")
    ##print(Textparam$method[[2]])
    ##print(paramlist)
    ##print(CNname[1])
    print(defaultValue)
    print(tips)
    lapply(1:(length(paramlist)-4), function(k) {
      # get min and max value
      #minmax <- c(round(min(dat[[varNames[k]]])), round(max(dat[[varNames[k]]])) )
      print(attributes(defaultValue[k]))
      print(tips[[paramlist[k]]])
      
      fluidRow(
        column(12, 
               bsTooltip(paramlist[k], tips[[paramlist[k]]],
                         "top", options = list(container = "body")),
               if (paramlist[k] == "level") {
                 # a slider range will created only is the variable is selected
                 ##print(paramlist[k])
                 sliderInput(inputId = "level",
                             label = paste(CNname[k+1],":"),
                             min = 0.1,
                             max = 0.99,
                             value = c(0.5,0.95)
                          )
               } 
               else if(paramlist[k] == "assurance"){
                 sliderInput(inputId = "assurance",
                             label = paste(CNname[k+1],":"),
                             min = 0.01,
                             max = 0.99,
                             value = c(0.01,0.05)
                 )
               }else if(paramlist[k] == "N"){
                 sliderInput(inputId = "N",
                             label = paste(CNname[k+1],":"),
                             min = 10,
                             max = 1000,
                             value = c(50,200)
                 )
               }
               else{
                 if (paramlist[k]  %in% names(Textparam)){
                   for (i in 2:length(Textparam[[paramlist[k] ]])){
                      #print(Textparam[[paramlist[k] ]][[i]])  
                   }
                   #print(length(Textparam[[paramlist[k] ]]))
                   selectInput(paramlist[k], paste(CNname[k+1],":"),
                      Textparam[[paramlist[k] ]][2:length(Textparam[[paramlist[k] ]])]
                   )
                 }else{
                  # numericInput(inputId = paramlist[k], paste(CNname[k+1],":"), 0.5)
                 #}
                 # otherwise uses single value with a default value
                 numericInput(inputId = paramlist[k], paste(CNname[k+1],":"), defaultValue[[paramlist[k] ]])
                 }
               }
              
        )
      )
    })
  })
  
  #列出散点图
  #output$plot1 = renderPlot({
    # get the correct id name for the current slider
   
  # id <- paste0("slider_", rv$selected)
   # cat("id", id, "\n")
    # get the value from the input
    #val = input[[id]]
    # plot all points of the selected variable
    #plot(dat[,rv$selected])
    # fill out the points that are greater or equal to the value
    #points(dat[ dat[,rv$selected] >= val, rv$selected, drop = FALSE], pch = 19, cex = 2)
  #})
  
  ##print(val)
  
  menuItemSelected <- reactive({
    slt <- input$sidebarMenu
    
  })
  
  rvSelected <- reactive({
    slt <- rv$sidebarMenu
    
  })
  
  allChoices <- reactive({
    func <- switch( rv$stats,
                   "mu.varunknow.seq" = mu.varunknow.seq,
                   "mu.varknown.seq" = mu.varknown.seq,
                   "mudiff.equalvar.seq" = mudiff.equalvar.seq,
                   "mudiff.unequalvar.seq" = mudiff.unequalvar.seq,
                   "propdiff.seq" = propdiff.seq,
                   #"propmbl.seq" = propmbl.seq,
                   "inv.ANDks" = inv.ANDks,
                   "inv.ANDus" = inv.ANDus,
                   "inv.AWBD" = inv.AWBD,
                   "inv.ASDExp" = inv.ASDExp,
                   "BFt" = BFt,
                   "BFpp" = BFpp,
                   "BFsur" = BFsur,
                   "ppSSR" = ppSSR,
                   "PSSD" = PSSD
                   
                   
                   
    )
    funclist = formals(func)
    #list中的CNname为对应的中文名称
    #列出所有可能参数的中文名称或者在fun中定义
    #列出所有的参数
    paramlist <- names(funclist)
    val <- list()
    for (i in paramlist){
      id <- i
      cat("id", id, "\n")
      print(input[[id]])
    # get the value from the input
      val[id] = input[[id]]
    # plot all points of the selected variable
    }
    print(input)
    print("huhu")
    print(input$sidebarMenu)
    #print(input$stats)
    sampleData <- switch( rv$stats,
                         "mu.varunknow.seq" = mu.varunknow.seq(input$len, input$alpha, input$beta, input$n0, seq(input$level[1], input$level[2], length.out = 10), input$method),
                         "mu.varknown.seq" = mu.varknown.seq(input$len, input$lambda, input$n0, seq(input$level[1], input$level[2], length.out = 10)),
                         "mudiff.equalvar.seq" = mudiff.equalvar.seq(input$len, input$alpha, input$beta, input$n01, input$n02, seq(input$level[1], input$level[2], length.out = 10),method=input$method),
                         "mudiff.unequalvar.seq" = mudiff.unequalvar.seq(len=input$len, alpha1=input$alpha1, beta1=input$beta1, alpha2=input$alpha2, beta2=input$beta2, n01=input$n01, n02=input$n02,seq(input$level[1], input$level[2], length.out = 10),method = input$method),
                         "propdiff.seq" = propdiff.seq(len=input$len, c1=input$c1, d1=input$d1, c2=input$c2, d2=input$d2, level = seq(input$level[1], input$level[2], length.out = 10), method=input$method),
                         #"propmbl.seq" = propmbl.seq(len=input$len, c1=input$c1, d1=input$d1, c2=input$c2, d2=input$d2, level=seq(input$level[1], input$level[2], length.out = 10), method=input$method),
                         "inv.ANDks" = inv.ANDks(prior_mean=input$prior_mean,prior_sd=input$prior_sd,prior_size=input$prior_size,N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),alpha=input$alpha),
                         "inv.ANDus" = inv.ANDus(prior_mean=input$prior_mean,prior_sd=input$prior_sd,prior_size=input$prior_size,N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),alpha=input$alpha),
                         "inv.AWBD" = inv.AWBD(beta1param1=input$beta1param1, beta1param2=input$beta1param2,beta2param1=input$beta2param1, beta2param2=input$beta2param2, N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),alpha=input$alpha),
                         "inv.ASDExp" = inv.ASDExp(t0=input$t0,a=input$a,b=input$b,m=input$m,v=input$v,R=input$R,TT=input$TT,N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),alpha = input$alpha),
                         "BFt" = BFt(boundary=input$boundary, ES=input$ES, alternative=input$alternative,priortype =input$priortype, prior_location=input$prior_location,prior_scale = input$prior_scale,prior_df = input$prior_df,N=ceiling(seq(input$N[1], input$N[2], length.out = 5))),
                         "BFpp" =  BFpp(N=ceiling(seq(input$N[1], input$N[2], length.out = 5)),boundary=input$boundary, ES=input$ES, alternative=input$alternative,prior_mean = input$prior_mean, prior_sd = input$prior_sd,effecttype=input$effecttype,p1=input$p1),
                         "BFsur" = BFsur(N=ceiling(seq(input$N[1], input$N[2], length.out = 5)), boundary=input$boundary, beta0=input$beta0,prior.maxL=input$prior.maxL,prior.a0=input$prior.a0,prior.b0=input$prior.b0,Tsim=input$Tsim),
                         "ppSSR" = ppSSR(N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),y1=input$y1,n1=input$n1,y2=input$y2,n2=input$n2,delta=input$delta,eta=input$eta,a1=input$a1, b1=input$b1, a2=input$a2, b2=input$b2),
                         "PSSD" = PSSD(a_a=input$a_a, b_a=input$b_a, n_d=input$n_d, p_d=input$p_d, p0=input$p0, delta=input$delta, lambda=input$lambda, N=ceiling(seq(input$N[1], input$N[2], length.out = 10)))
    )
    print(sampleData)
    sampleData
    # Require that all input$checkboxgrp and 
    # the last generated numericInput are available.
    # (If the  last generated numericInput is available (is not NULL),
    # then all previous are available too)
    
    # "eval(parse(text = paste0("input$", input$checkboxgrp))))" yields
    # a value of the last generated numericInput. 
    
    # In this way we avoid multiple re-evaulation of allChoices() 
    # and errors
    #req(input$len)
      # get the correct id name for the current slider
      #id <- paste0("slider_", rv$selected)
      #cat("id", id, "\n")
      # get the value from the input
      #val = input[[id]]
      # plot all points of the selected variable
      #plot(dat[,rv$selected])
      # fill out the points that are greater or equal to the value
      #points(dat[ dat[,rv$selected] >= val, rv$selected, drop = FALSE], pch = 19, cex = 2)
  })
  
  
  reference <- reactive({
    a <- allChoices()$conclusion
    
  })
  
  output$summary <- renderText({
    reference()
  })
  
  output$pdfviewer <- renderUI({
    tags$iframe(id = "localFile",
                style = "height: 900px; width: 100%; scrolling = yes",
                # src = paste0("localfiles/", fileName)
                src = paste0("localfiles/",rv$stats,".pdf")
    )
  })  
  
  
  plotdata<- reactive({
    sampleData <- switch( rv$stats,
                         "mu.varunknow.seq" = mu.varunknow.seq(input$len, input$alpha, input$beta, input$n0, seq(input$level[1], input$level[2], length.out = 10), input$method),
                         "mu.varknown.seq" = mu.varknown.seq(input$len, input$lambda, input$n0, seq(input$level[1], input$level[2], length.out = 10)),
                         "mudiff.equalvar.seq" = mudiff.equalvar.seq(input$len, input$alpha, input$beta, input$n01, input$n02, seq(input$level[1], input$level[2], length.out = 10),method=input$method),
                         "mudiff.unequalvar.seq" = mudiff.unequalvar.seq(len=input$len, alpha1=input$alpha1, beta1=input$beta1, alpha2=input$alpha2, beta2=input$beta2, n01=input$n01, n02=input$n02,seq(input$level[1], input$level[2], length.out = 10),method = input$method),
                         "propdiff.seq" = propdiff.seq(len=input$len, c1=input$c1, d1=input$d1, c2=input$c2, d2=input$d2, level = seq(input$level[1], input$level[2], length.out = 10), method=input$method),
                         #"propmbl.seq" = propmbl.seq(len=input$len, c1=input$c1, d1=input$d1, c2=input$c2, d2=input$d2, level=seq(input$level[1], input$level[2], length.out = 10), method=input$method),
                         "inv.ANDks" = inv.ANDks(prior_mean=input$prior_mean,prior_sd=input$prior_sd,prior_size=input$prior_size,N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),alpha=input$alpha),
                         "inv.ANDus" = inv.ANDus(prior_mean=input$prior_mean,prior_sd=input$prior_sd,prior_size=input$prior_size,N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),alpha=input$alpha),
                         "inv.AWBD" = inv.AWBD(beta1param1=input$beta1param1, beta1param2=input$beta1param2,beta2param1=input$beta2param1, beta2param2=input$beta2param2, N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),alpha=input$alpha),
                         "inv.ASDExp" = inv.ASDExp(t0=input$t0,a=input$a,b=input$b,m=input$m,v=input$v,R=input$R,TT=input$TT,N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),alpha = input$alpha),
                         "BFt" = BFt(boundary=input$boundary, ES=input$ES, alternative=input$alternative,priortype =input$priortype, prior_location=input$prior_location,prior_scale = input$prior_scale,prior_df = input$prior_df,N=ceiling(seq(input$N[1], input$N[2], length.out = 5))),
                         "BFpp" =  BFpp(N=ceiling(seq(input$N[1], input$N[2], length.out = 5)),boundary=input$boundary, ES=input$ES, alternative=input$alternative,prior_mean = input$prior_mean, prior_sd = input$prior_sd,effecttype=input$effecttype,p1=input$p1),
                         "BFsur" = BFsur(N=ceiling(seq(input$N[1], input$N[2], length.out = 5)), boundary=input$boundary, beta0=input$beta0,prior.maxL=input$prior.maxL,prior.a0=input$prior.a0,prior.b0=input$prior.b0,Tsim=input$Tsim),
                         "ppSSR" = ppSSR(N=ceiling(seq(input$N[1], input$N[2], length.out = 10)),y1=input$y1,n1=input$n1,y2=input$y2,n2=input$n2,delta=input$delta,eta=input$eta,a1=input$a1, b1=input$b1, a2=input$a2, b2=input$b2),
                         "PSSD" = PSSD(a_a=input$a_a, b_a=input$b_a, n_d=input$n_d, p_d=input$p_d, p0=input$p0, delta=input$delta, lambda=input$lambda, N=ceiling(seq(input$N[1], input$N[2], length.out = 10)))
    )
  })
  
  
 
  vis <- reactive({      
    #yvar_name <- names(axis_vars_y)[axis_vars_y == input$yvar]     
    #yvar <- prop("y", as.symbol(input$yvar))      
    a = data.frame(allChoices())
    #print(a)
    p <- ggplot(a, aes(x = level, y = sample)) + geom_point(text = y)
    ggplotly(p,tooltip = "text")
    ##print(a$sample)
    #if (is.null(input$len)){
     # mtcars %>%       
      #  ggvis(x = ~mpg, y = ~wt) %>%           
       # layer_points() 
    #}else{
     # ggplot(a, aes(x = level, y = sample)) + geom_point()
        #p =ggvis(x = a$level, y = a$sample)         
        #layer_points() 
    #}
  })   
  
# output$plot1 <- renderPlot({
#   a = data.frame(allChoices())
#   format(a, digits = 3)
#   print(a)
#   if("sample_c"  %in% names(a)){
#     ggplot(a, aes(x=level)) + 
#     geom_point(aes(y=sample), color="red",size=5) + 
# geom_text( hjust = 0, nudge_x = 0.05,label = a$sample)+
# geom_line(aes(y=BROWN.P, , color="cyan")) +
 #     geom_point(aes(y=sample_c),size=5)
 #     #geom_point(size = 5)+
 #     #geom_text( hjust = 0, nudge_x = 0.05,label = a$sample)
 #     #geom_line(aes(y=MI.P, color="red"))
 #   }else{
 #     ggplot(a, aes(x = level, y = sample)) + geom_point(size = 5)+geom_text( hjust = 0, nudge_x = 0.05,label = a$sample)
 #   }
 #   
 #   
 #   #ggplot(a, aes(x = level, y = sample)) + geom_point(size = 5)+geom_text( hjust = 0, nudge_x = 0.05,label = a$sample)
 #   
 # })
 
 output$pi <- renderPlotly({
   a = data.frame(allChoices())
   format(a, digits = 3)
   if("assurance"  %in% names(a)){
      plot_ly(a, x = ~assurance, y = ~sample) %>%
     add_lines()
   }else if("power"  %in% names(a)){
     plot_ly(a, x = ~power, y = ~sample) %>%
       add_lines()
   }else{
     plot_ly(a, x = ~level, y = ~sample) %>%
       add_lines()
   }
 })
  
 output$mytable <- DT::renderDataTable({
   print("55")
   print(menuItemSelected())
   print(rvSelected())
   format(data.frame(allChoices()[1:length(allChoices())-1]), digits = 3)
 })
 #output$table <- renderTable(data.frame(allChoices()))
  #allChoices() %>% 
   # ggvis(~sample, ~level) %>%
    #layer_points() %>%
    #bind_shiny("plot1")
  
 
})