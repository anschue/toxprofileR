#' Run fingerprint browser
#'
#' @return opens a shiny app
#' @export
tfp_browser <- function() {

  # load libraries --------------------------------------------------------------
  library("shiny")
  library("plotly")

  ###
  colvec_special <- c(0, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 1)
  scale_range <- c(0.1, 2.5)

  breaksfunction <- function(xlim) {
    df <- (xlim[2] / xlim[1])^(1 / 4)
    breaks <- xlim[2] / df^(c(1, 2, 3))
    return(sort(round(breaks, digits = 1)))
  }
  ###


  # ui --------------------------------------------------------------------------

  ui <- fluidPage(
      titlePanel(
          title = div(
              "Toxicogenomic Fingerprint Browser"
          ),
          windowTitle = "Toxicogenomic Fingerprint Browser"
      ),



      fluidRow(
          column(
              2,
              style = "background-color:#92B6D5;",

              tabsetPanel(
                  type = "tabs",

                  tabPanel(
                      "Upload",
                      br(),
                      # Selector for file upload
                      fileInput("file", "Choose RDS file")
                  ),

                  tabPanel(
                      "Treatment",
                      br(),
                      # Select substance ------------------------------------
                      uiOutput("substancechoice"),

                      # Select concentration --------------------------------
                      uiOutput("concslider"),

                      # Select time point -----------------------------------
                      uiOutput("timeslider")
                  ),
                  tabPanel(
                      "Gene search",
                      br(),
                      # Search for Gene -------------------------------------
                      uiOutput("genesearcher")
                  )


              )


          ),

          column(1,
                 # Map legend ------------------------------------------
                 plotOutput("maplegend", height = "auto")),

          column(
              5,
              #plot map --------------------------------------------
              h3("Fingerprint"),
              plotlyOutput(
                  outputId = "plot1",
                  inline = T,
                  height = "auto"
              )
          ),



          column(
              4,
              style = "background-color:#F5F5F5;",
              # plot node -------------------------------------------
              h3("Dynamic response of selected toxnode"),
              plotOutput(outputId = "nodeplot",
                         height = "auto")
          )
      ),

      # display node information ----------------------------------------------------------
      fluidRow(
          column(
              12,
              style = "background-color:#F5F5F5;",
              htmlOutput("tableheader"),
              DT::dataTableOutput("nodeInfo", width = "auto")
          )
      )
  )


  # server ----------------------------------------------------------------------

  server <- function(input, output, session) {
    options(shiny.maxRequestSize = 1000 * 1024^2)

    # get session data for plot scaling
    cdata <- session$clientData


    observeEvent(input$file, {
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
                  label = "Compound",
                  choices = as.list(unique(
                      shiny_data$targets_all$substance
                  )),
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
              sort(log10(
                  unique(shiny_data$targets_all$concentration_umol_l[shiny_data$targets_all$substance == input$selectedsubstance1])
              ))

          concentrations[is.infinite(concentrations)] <- 0

          args       <-
              list(
                  inputId = "concselect1",
                  label = "Concentration [log(µmol/L)] :",
                  ticks = concentrations,
                  value = concvalues$defaultconc
              )

          args$min   <- 1
          args$max   <- length(args$ticks)

          ticks <-
              paste0(round(args$ticks, digits = 2), collapse = ',')
          args$ticks <- T
          html  <- do.call('sliderInput', args)
          html$children[[2]]$attribs[['data-values']] <- ticks

          html
      })

      # make time point slider ------------------------------------------------------------
      output$timeslider <- renderUI({
          times <-
              sort(unique(shiny_data$targets_all$time_hpe[shiny_data$targets_all$substance == input$selectedsubstance1]))

          argstime <- list(
              inputId = "timeselect1",
              label = "Time [hpe]:",
              ticks = times,
              value = timevalues$defaulttime
          )

          argstime$min <- 1
          argstime$max <- length(argstime$ticks)

          tickstime <-
              paste0(round(argstime$ticks, digits = 2), collapse = ',')
          argstime$ticks <- T
          htmltime  <- do.call('sliderInput', argstime)
          htmltime$children[[2]]$attribs[['data-values']] <-
              tickstime

          htmltime
      })

      # make Gene-Searcher -----------------------------------------------------------------
      output$genesearcher <- renderUI({
          selectList <-
              data.table::data.table(
                  genenames = shiny_data$nodeannotation$external_gene_name,
                  toxnode = shiny_data$nodeannotation$toxnode
              )

          argssearcher <-
              list(
                  inputId = "geneselect",
                  label = "Search gene name",
                  choices = c(Choose = '', selectList),
                  selected = NULL,
                  multiple = FALSE
              )

          searcher <- do.call('selectizeInput', argssearcher)
          searcher

      })

      # select node ------------------------------------------------------------------------
      nodecoord <- reactiveValues(x = 0, y = 0)

      click_data <- reactive({
          d <- event_data("plotly_click")
          if (is.null(d)) {
              d <- list(x = 0, y = 0)
          }
          d
      })

      observeEvent(input$geneselect, {
          if (is.null(input$geneselect)) {
              nodecoord$x <- 0
              nodecoord$y <- 0
          } else{
              message(input$geneselect)
              nodecoord$x <-
                  shiny_data$grid$pts[shiny_data$nodeannotation$toxnode[shiny_data$nodeannotation$external_gene_name == input$geneselect][1], "x"]
              nodecoord$y <-
                  shiny_data$grid$pts[shiny_data$nodeannotation$toxnode[shiny_data$nodeannotation$external_gene_name == input$geneselect][1], "y"]
          }
      })

      observeEvent(click_data(), {
          nodecoord$x <- click_data()$x
          nodecoord$y <- click_data()$y
      })

      observeEvent(session, {
          # get map legend ------------------------------------------------------------
          output$maplegend <- renderPlot({
              shiny_data$maplegend
          },
          width = min(100, cdata$output_maplegend_width),
          height = min(400, cdata$output_maplegend_width * 4))

          # plot map --------------------------------------------------------------------------
          output$plot1 <- renderPlotly({
              if (is.null(input$concselect1) | is.null(input$timeselect1)) {
                  ggplot2::ggplot(data = data.frame(
                      x = 1,
                      y = 1,
                      label = "Please select treatment"
                  )) +
                      geom_text(aes(
                          x = x,
                          y = y,
                          label = label
                      )) +
                      theme_bw() +
                      theme(
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          axis.text = element_blank(),
                          axis.title = element_blank(),
                          axis.ticks = element_blank(),
                          legend.position = "none"
                      )
              } else{
                  conc_id <- as.numeric(input$concselect1) + 1
                  time_id <- as.numeric(input$timeselect1) + 1


                  maptoplot1 <-
                      shiny_data$targets_all$name[shiny_data$targets_all$substance == input$selectedsubstance1 &
                                                      shiny_data$targets_all$concentration_umol_l == sort(unique(
                                                          shiny_data$targets_all$concentration_umol_l[shiny_data$targets_all$substance == input$selectedsubstance1]
                                                      ))[conc_id] &
                                                      shiny_data$targets_all$time_hpe == sort(unique(shiny_data$targets_all$time_hpe[shiny_data$targets_all$substance == input$selectedsubstance1]))[time_id]]


                  if (is.na(nodecoord$x) | nodecoord$x == 0) {
                      plot_1 <- shiny_data$plotlist[[maptoplot1]]
                  } else{
                      plot_1 <-
                          shiny_data$plotlist[[maptoplot1]] + geom_point(
                              data = data.frame(x = nodecoord$x, y = nodecoord$y),
                              aes(x = x, y = y),
                              size = 4,
                              stroke = 2,
                              pch = 1,
                              col = "black"
                          )
                  }

                  plot_1 %>%
                      ggplotly(
                          tooltip = "all",
                          width =  min(cdata$output_plot1_width, 750),
                          height = min(cdata$output_plot1_width, 750)
                      ) %>%
                      layout(dragmode = "select") %>%
                      config(displayModeBar = FALSE, displaylogo = FALSE)

              }

          })

          output$mapheight <- reactive({
              cdata$output_plot1_width
          })


          # plot node measurements --------------------------------------------------

          output$nodeplot <- renderPlot({
              if (is.na(nodecoord$x) | nodecoord$x == 0) {
                  ggplot2::ggplot(data = data.frame(
                      x = 1,
                      y = 1,
                      label = "Click on fingerprint to select toxnode"
                  )) +
                      geom_text(aes(
                          x = x,
                          y = y,
                          label = label
                      ), size = 8) +
                      theme_bw() +
                      theme(
                          panel.background = element_rect(fill = "#F5F5F5"),
                          plot.background = element_rect(fill = "#F5F5F5"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          axis.text = element_blank(),
                          axis.title = element_blank(),
                          axis.ticks = element_blank(),
                          legend.position = "none"
                      )
              } else{
                  nodenr <-
                      which(
                          shiny_data$grid$pts[, "x"] == nodecoord$x &
                              shiny_data$grid$pts[, "y"] == nodecoord$y
                      )

                  substance <- input$selectedsubstance1

                  D_measured <-
                      shiny_data[["D_measured_all"]][shiny_data[["D_measured_all"]][, "substance"] ==
                                                         substance &
                                                         shiny_data[["D_measured_all"]][, "nodeID"] == nodenr, ]

                  if (dim(D_measured)[1] > 0) {
                      poiplot_hill <-
                          ggplot2::ggplot(data = D_measured,
                                          aes(x = concentration_umol_l, y = logFC)) +
                          geom_point(aes(colour = probe_id),
                                     size = 1,
                                     lwd = 0) +
                          facet_wrap( ~ factor(time_hpe),
                                      nrow = 1,
                                      scales = "fixed") +
                          geom_line(
                              data = shiny_data[["D_fit_all"]][shiny_data[["D_fit_all"]][, "substance"] ==
                                                                   substance &
                                                                   shiny_data[["D_fit_all"]][, "node"] == nodenr, ],
                              aes(x = concentration_umol_l, y = logFC_hill),
                              lwd = 1.5
                          ) +
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
                          theme_bw() +
                          theme(
                              axis.title.y = element_text(size = 16, face = "bold"),
                              axis.title.x = element_text(size = 16, face = "bold"),
                              axis.text.x = element_text(
                                  angle = 70,
                                  vjust = 1,
                                  hjust = 1,
                                  size = 14
                              ),
                              axis.text.y = element_text(size = 14),
                              strip.text = element_text(size = 11, face = "bold"),
                              legend.position = "top",
                              panel.background = element_rect(fill = "#F5F5F5"),
                              plot.background = element_rect(fill = "#F5F5F5")
                          ) +
                          scale_colour_discrete(name = "ProbeID") +
                          ylab("logFC") +
                          xlab("\n\nexposure concentration [µM]") +
                          geom_hline(aes(yintercept = 0)) +
                          scale_x_log10(breaks = breaksfunction) +
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

                      poiplot_hill
                  } else {
                      ggplot2::ggplot(
                          data = data.frame(
                              x = 1,
                              y = 1,
                              label = "No measured data for selected toxnode.\n Please select different toxnode."
                          )
                      ) +
                          geom_text(aes(
                              x = x,
                              y = y,
                              label = label
                          ), size = 8) +
                          theme_bw() +
                          theme(
                              panel.background = element_rect(fill = "#F5F5F5"),
                              plot.background = element_rect(fill = "#F5F5F5"),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              axis.text = element_blank(),
                              axis.title = element_blank(),
                              axis.ticks = element_blank(),
                              legend.position = "none"
                          )

                  }
              }
          },
          height = cdata$output_nodeplot_width)

          # select node on map ----------------------------------------------------------------
          output$nodeInfo <- DT::renderDataTable({
              if (is.na(nodecoord$x) | nodecoord$x == 0) {
                  data.frame(Start = "Click on plot to select toxnode")
              } else{
                  nodenr <-
                      which(
                          shiny_data$grid$pts[, "x"] == nodecoord$x &
                              shiny_data$grid$pts[, "y"] == nodecoord$y
                      )
                  DT::datatable(
                      shiny_data$nodeannotation[shiny_data$nodeannotation$toxnode == nodenr, c(
                          "ProbeID",
                          "ensembl_gene_id",
                          "external_gene_name",
                          "description",
                          "name_1006",
                          "interpro_description",
                          "gene_biotype"
                      )],
                      colnames = c(
                          "ProbeID",
                          "Ensembl gene-id",
                          "Gene name",
                          "Description",
                          "GO Annotation",
                          "Interpro description",
                          "Gene biotype"
                      )
                  )
              }
          },
          options = list(dom = 't'))

          # render nodeID as table header ----------------------------------------------------
          output$tableheader <- renderPrint({
              if (is.na(nodecoord$x) | nodecoord$x == 0) {
                  h3(paste0("Gene table"))
              } else{
                  nodenr <-
                      which(
                          shiny_data$grid$pts[, "x"] == nodecoord$x &
                              shiny_data$grid$pts[, "y"] == nodecoord$y
                      )

                  h3(paste0("Gene table for toxnode #", nodenr))
              }
          })


      })

})
}
  shinyApp(ui = ui, server = server)
}
