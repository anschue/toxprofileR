#' Run fingerprint browser
#'
#' @return opens a shiny app
#' @export
fp_browser <- function(){

# load libraries --------------------------------------------------------------
library("shiny")
library("plotly")

###
colvec_special<-c(0,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,1)
scale_range<-c(0.1,2.5)

breaksfunction<-function(xlim){
    df<-(xlim[2]/xlim[1])^(1/4)
    breaks<-xlim[2]/df^(c(1,2,3))
    return(sort(round(breaks,digits = 1)))
}
###


# ui --------------------------------------------------------------------------

ui <- fixedPage(
    titlePanel("Toxicogenomic Fingerprint Browser"),
    
    fixedRow(column(12,
                    fixedRow(
                        column(
                            2,
                            
                            # select input file -----------------------------------
                            #Selector for file upload
                            fileInput('file', 'Choose RDS file'),
                            
                            # Select substance ------------------------------------
                            uiOutput("substancechoice"),
                            
                            # Select concentration --------------------------------
                            uiOutput("concslider"),
                            
                            # Select time point -----------------------------------
                            uiOutput("timeslider")
                        ),
                        
                        column(6,
                               # plot map --------------------------------------------
                               plotlyOutput(
                                   outputId = "plot1",
                                   height = 500,
                                   width = 500
                               )),
                        
                        column(4,
                               # plot node -------------------------------------------
                               plotlyOutput(
                                   outputId = "nodeplot",
                                   height = 400,
                                   width = 750
                               ))
                    ))),
    
    # display node information ----------------------------------------------------------
    dataTableOutput("clicknode")
    )

# server ----------------------------------------------------------------------

server <- function(input, output) {
    
    options(shiny.maxRequestSize=1000*1024^2)
    
    observeEvent(input$file,{
        infile <- input$file
        if (is.null(infile)) {
            # User has not uploaded a file yet
            return(NULL)
        }
    
    shiny_data <- readRDS(file = infile$datapath)
    

    
    
    output$substancechoice <- renderUI({
            argsselector <-
                list(
                    inputId = "selectedsubstance1",
                    label = h3("Select Compound"),
                    choices = as.list(unique(shiny_data$targets_all$substance)),
                    selected = 1
                )
            
            selector <- do.call('selectizeInput', argsselector)
            selector
            
        })
  
    # set default concentration and time point ------------------------------------------
    concvalues <- reactiveValues(defaultconc = 0)
    timevalues <- reactiveValues(defaulttime = 0)
    
    eventReactive(input$concselect1, {
        concvalues$defaultconc <- as.numeric(input$concselect1)
    })
    
    eventReactive(input$timeselect1, {
        timevalues$defaulttime <- as.numeric(input$timeselect1)
    })
    
    # make concentration slider ---------------------------------------------------------
    
    output$concslider <- renderUI({
        concentrations <-
            sort(log10(unique(shiny_data$targets_all$concentration_umol_l[
                shiny_data$targets_all$substance == input$selectedsubstance1 #&
                #shiny_data$targets_all$time_hpe == sort(unique(shiny_data$targets_all$time_hpe[shiny_data$targets_all$substance == input$selectedsubstance1]))[timevalues$defaulttime + 1]
                ])))
        
        concentrations[is.infinite(concentrations)] <- 0
        
        args       <-
            list(
                inputId = "concselect1",
                label = "concentration :",
                ticks = concentrations,
                value = concvalues$defaultconc
            )
        
        args$min   <- 1
        args$max   <- length(args$ticks)
        
        ticks <- paste0(round(args$ticks, digits = 2), collapse = ',')
        args$ticks <- T
        html  <- do.call('sliderInput', args)
        html$children[[2]]$attribs[['data-values']] <- ticks
        
        html
    })

    # make time point slider ------------------------------------------------------------
    output$timeslider <- renderUI({
        times <- sort(unique(
                shiny_data$targets_all$time_hpe[shiny_data$targets_all$substance == input$selectedsubstance1 #&
                                                #shiny_data$targets_all$concentration_umol_l == sort(unique(shiny_data$targets_all$concentration_umol_l[shiny_data$targets_all$substance == input$selectedsubstance1]))[concvalues$defaultconc + 1]
                                                ]))
        
        argstime <- list(
                inputId = "timeselect1",
                label = "time [hpf]:",
                ticks = times,
                value = timevalues$defaulttime
            )
        
        argstime$min <- 1
        argstime$max <- length(argstime$ticks)
        
        tickstime <- paste0(round(argstime$ticks, digits = 2), collapse = ',')
        argstime$ticks <- T
        htmltime  <- do.call('sliderInput', argstime)
        htmltime$children[[2]]$attribs[['data-values']] <- tickstime
        
        htmltime
    })
    
    # plot map --------------------------------------------------------------------------
    output$plot1 <- renderPlotly({
        
        if(is.null(input$concselect1)|is.null(input$timeselect1)){
            ggplot2::ggplot(data = data.frame(x = 1, y = 1, label = "Please select treatment"))+
                geom_text(aes(x = x, y = y, label = label)) +
                theme_bw()+
                theme(plot.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.text = element_blank(),
                      axis.title = element_blank(),
                      axis.ticks = element_blank(),
                      legend.position = "none")
        }else{
        
        conc_id <- as.numeric(input$concselect1) + 1
        time_id <- as.numeric(input$timeselect1) + 1
      
        
        maptoplot1 <-
            shiny_data$targets_all$name[shiny_data$targets_all$substance == input$selectedsubstance1 &
                                        shiny_data$targets_all$concentration_umol_l == sort(unique(shiny_data$targets_all$concentration_umol_l[shiny_data$targets_all$substance == input$selectedsubstance1]))[conc_id] &
                                        shiny_data$targets_all$time_hpe == sort(unique(shiny_data$targets_all$time_hpe[shiny_data$targets_all$substance == input$selectedsubstance1]))[time_id]]
        
        plot_1 <- shiny_data$plotlist[[maptoplot1]]
        
        plot_1 %>%
            ggplotly(tooltip = "all") %>%
            layout(dragmode = "select")
        
        }
        
    })
    
    
    # plot node measurements --------------------------------------------------
    
    output$nodeplot <- renderPlotly({
        d <- event_data("plotly_click")
        
        if(is.null(d)){
            ggplot2::ggplot(data = data.frame(x = 1, y = 1, label = "Click on fingerprint to select toxnode"))+
                geom_text(aes(x = x, y = y, label = label)) +
                theme_bw()+
                theme(plot.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                legend.position = "none")
        }else{
            
        nodenr <-
            which(shiny_data$grid$pts[, "x"] == d$x &
                  shiny_data$grid$pts[, "y"] == d$y)
        
        substance <- input$selectedsubstance1
       
        
    
        poiplot_hill <-
            ggplot2::ggplot(data = shiny_data[["D_measured_all"]][shiny_data[["D_measured_all"]][, "substance"] ==
                                                                        substance &
                                                                        shiny_data[["D_measured_all"]][, "nodeID"] == nodenr, ], aes(x = concentration_umol_l, y = logFC)) +
            geom_point(aes(colour = probe_id),
                       size = 1,
                       lwd = 0) + 
            facet_wrap( ~ factor(time_hpe), nrow = 1, scales = "fixed") +
            geom_line(data = shiny_data[["D_fit_all"]][shiny_data[["D_fit_all"]][, "substance"] ==
                                                             substance &
                                                             shiny_data[["D_fit_all"]][, "node"] == nodenr, ],
                      aes(x = concentration_umol_l, y = logFC_hill),
                      lwd = 1.5) +
            geom_ribbon(
                data = shiny_data[["D_fit_all"]][shiny_data[["D_fit_all"]][, "substance"] ==
                                                       substance &
                                                       shiny_data[["D_fit_all"]][, "node"] == nodenr, ],
                aes(
                    x = concentration_umol_l,
                    ymin = lwr_hill,
                    ymax = upr_hill
                ),
                alpha = 0.3,
                inherit.aes = F
            ) +
            #theme_bw()+
            theme(
                legend.position = "right",
                axis.title.y = element_text(size = 16, face = "bold"),
                axis.title.x = element_text(size = 16, face = "bold"),
                axis.text.x = element_text(
                    angle = 70,
                    vjust = 1,
                    hjust = 1,
                    size = 14
                ),
                axis.text.y = element_text(size = 14),
                strip.text = element_text(size = 11, face = "bold")
            ) +
            ylab("logFC") +
            #xlab(expression("exposure concentration [" * mu * "M]")) +
            xlab("\n\nexposure concentration [µM]") +
            geom_hline(aes(yintercept = 0)) +
            scale_x_log10(breaks = breaksfunction) +
            #ggtitle(paste(POI,"hill"))+
            geom_hline(
                data = shiny_data[["Control_CIs_all"]][shiny_data[["Control_CIs_all"]][, "substance"] ==
                                                            substance &
                                                            shiny_data[["Control_CIs_all"]][, "node"] == nodenr, ],
                aes(yintercept = min),
                lty = 2,
                lwd = 0.75
            ) +
            geom_hline(
                data = shiny_data[["Control_CIs_all"]][shiny_data[["Control_CIs_all"]][, "substance"] ==
                                                            substance &
                                                            shiny_data[["Control_CIs_all"]][, "node"] == nodenr, ],
                aes(yintercept = max),
                lty = 2,
                lwd = 0.75
            )
        
        
        
        poiplot_hill %>%
            ggplotly(tooltip = "shape") %>%
            layout(dragmode = "select")
        }
    })
    
    # select node on map ----------------------------------------------------------------
    output$clicknode <- renderDataTable({
        d <- event_data("plotly_click")
        if (is.null(d)) {
            data.frame(Start = "Click on plot to select toxnode")
        } else{
            nodenr <-
                which(shiny_data$grid$pts[, "x"] == d$x &
                      shiny_data$grid$pts[, "y"] == d$y)
            shiny_data$nodeannotation[shiny_data$nodeannotation$toxnode == nodenr,]
        }
    })
    
    
    output$selections <- renderPrint({
        paste(input$selectedsubstance1,
              input$concselect1,
              input$timeselect1)
    })
})
}

shinyApp(ui = ui, server = server)

}
