## ui.R ##
library(shinydashboard)
library(ggvis)
library(DT)
library(plotly)
xx <- iconv('非序贯', 'GB2312', 'UTF-8') 
dashboardPage(
    dashboardHeader(title = "贝叶斯样本量计算", 
                    titleWidth = 300,
                    dropdownMenu(type = "messages",
                                 messageItem(
                                   from = "Sales Dept",
                                   message = "Sales are steady this month."
                                 ),
                                 messageItem(
                                   from = "New User",
                                   message = "How do I register?",
                                   icon = icon("question"),
                                   time = "13:45"
                                 ),
                                 messageItem(
                                   from = "Support",
                                   message = "The new server is ready.",
                                   icon = icon("life-ring"),
                                   time = "2014-12-01"
                                 )
                    )
                    ),
    #布局：侧边栏
    dashboardSidebar(
      width = 300,
      sidebarMenu(
        id = "sbMenu",
        #导航页
        #sidebarMenuOutput("menu"),
        #sidebarSearchForm(textId = "searchText", buttonId = "searchButton",
                        #  label = "Search..."),
        #menuItem("样本量计算工具", tabName = "dashboard", icon = icon("dashboard")),
        #menuItem("Widgets", tabName = "widgets", icon = icon("th")),
        menuItem(
          "选择变量",
          tabName = "normal",
          icon = icon("address-card"),
          radioButtons("trial", "试验类型:",c("非序贯" = "one","序贯" = "two"),inline = TRUE),
          conditionalPanel(
            condition = "input.trial == 'one'",
            selectInput(
              "datatype", "数据类型",
              c("连续型"="normal", "二分类"="binary", "生存"="survival")
            ),
            conditionalPanel(
              condition = "input.datatype == 'normal'",
              selectInput("stats", "方法类型:",
                          c("单组方差未知" = "mu.varunknow.seq",
                            "单组方差已知" = "mu.varknown.seq",
                            "两组方差不齐"="mudiff.unequalvar.seq",
                            "两组方差齐" = "mudiff.equalvar.seq",
                            "assurance/两组方差已知"="inv.ANDks",
                            "assurance/两组方差未知"="inv.ANDus",
                            "贝叶斯因子法/两独立样本"= "BFt"
              )
            )
            ),
            conditionalPanel(
              condition = "input.datatype == 'binary'",
              selectInput("stats1", "方法类型:",
                          c(
                            "二分类/区间积分法" = "propdiff.seq",
                            "assurance/两组二分类"="inv.AWBD",
                            "贝叶斯因子法/二分类"= "BFpp")
              )
           ),
           conditionalPanel(
             condition = "input.datatype == 'survival'",
             selectInput("stats2", "方法类型:",
                         c( "assurance/两组生存指数分布"="inv.ASDExp")
             )
           ),
           ),
          conditionalPanel(
            condition = "input.trial == 'two'",
            selectInput("stats3", "方法类型:",
                        c("基于预测概率的SSR"="ppSSR",
                          "基于预测概率的SSD"="PSSD")
            )
            
          )
          
          # selectInput("trialType", "shuju类型:",
          #             c("差异性" = "test",
          #               "非劣效" = "no-inferior",
          #               "等效"="equivance")
          # ),
          # conditionalPanel(
          #   condition = "input.plotType == 'hist'",
          #   selectInput(
          #     "breaks", "Breaks",
          #     c("Sturges", "Scott", "Freedman-Diaconis", "[Custom]" = "custom")
          #   )
          # ),
          # conditionalPanel(
          #   condition = "input.plotType == 'scatter'",
          #   selectInput(
          #     "breaks", "dd",
          #     c("Sturges", "Scott", "Freedman-Diaconis", "[Custom]" = "custom")
          #   )
          # ),
          # selectInput("stats", "方法类型:",
          #             c("单组方差未知" = "mu.varunknow.seq",
          #               "单组方差已知" = "mu.varknown.seq",
          #               "两组方差不齐"="mudiff.unequalvar.seq",
          #               "两组方差齐" = "mudiff.equalvar.seq",
          #               "二分类/区间积分法" = "propdiff.seq",
          #               #"二分类/混合贝叶斯法" = "propmbl.seq",
          #               "assurance/两组方差已知"="inv.ANDks",
          #               "assurance/两组方差未知"="inv.ANDus",
          #               "assurance/两组二分类"="inv.AWBD",
          #               "assurance/两组生存指数分布"="inv.ASDExp",
          #               "贝叶斯因子法/两独立样本"= "BFt",
          #               "贝叶斯因子法/二分类"= "BFpp",
          #               "贝叶斯因子法/生存分析"="BFsur",
          #               "基于预测概率的SSR"="ppSSR",
          #               "基于预测概率的SSD"="PSSD")
          # )
          #div(class ="form-group shiny-input-container", submitButton("提交", icon("refresh")))
        )
        # menuItem(
        #   "生存结局变量",
        #   tabName = "normal",
        #   icon = icon("address-card"),
        #   selectInput("group", "组别:",
        #               c("单组" = "single",
        #                 "两组" = "two")
        #   ),
        #   checkboxGroupInput(
        #     inputId = "genderInput",
        #     label = "",
        #     choices = "",
        #     selected = "",
        #     inline = TRUE
        #   ),
        #   sliderInput(
        #     inputId = "ageInput",
        #     label = "Age",
        #     value = c(10,20),
        #     min = min(10),
        #     max = max(30),
        #     step = 1,
        #     sep = ""
        #   )
        # )
      )
    ),
    #布局：右侧面板
    dashboardBody(
      #导航栏对应的target
      
      fluidRow(
        column(8,
               tabBox(
                 title = "结果",
                 # The id lets us use input$tabset1 on the server to find the current tab
                 id = "tabset1", width = "1000px",
                 tabPanel("样本量",height = "250px",
                          weight="100%",
                          #box(plotOutput("plot1", height = 250)),
                          
                          #box(
                          # title = "Controls",
                          # sliderInput("slider", "Number of observations:", 1, 100, 50)
                          #)
                          #plotOutput("plot1"),
                          plotlyOutput( "pi"),
                          hr(),
                          DT::dataTableOutput("mytable")
                          #dataTableOutput('table'),
                          #tableOutput('table'),
                          #wellPanel(
                          # span("Number of movies selected:",
                          #     textOutput("n_movies")
                          #)
                          #)
                 ),
                 tabPanel("帮助文档", 
                          fluidRow(column(10,htmlOutput("pdfviewer")
                          )
                          )
                 )
               )
        ),
        column(4,
               tabItem(tabName = "adjustparam",
                       
                       shinydashboard::box(width = "100%",
                           title = "参数表", status = "primary", solidHeader = TRUE,
                           "Box content here", br(), "More box content", collapsible = TRUE,
                           uiOutput("param")
                       )
               )
        ),
        column(4,
               tabItem(tabName = "summary",
                       
                       shinydashboard::box(width = "100%",
                                           title = "结果总结", status = "primary", solidHeader = TRUE,
                                           collapsible = TRUE,
                                           textOutput("summary")
                       )
               )
        ),
      )
  )
)