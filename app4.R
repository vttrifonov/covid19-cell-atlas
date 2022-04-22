{
    library(data.table)
    library(magrittr)
    library(shiny)
    library(ggplot2)
    library(DT)
    library(reticulate)
    library(glue)
    library(kableExtra)
    source('helpers.R')
    source('common/module.R')
    load.module('common')
    py_run_string('from covid19_cell_atlas._app4 import app4 as app')
    py_app <- py_eval('app')
}


{
    app4 <- obj()

    app4 %>%
        lazy_prop(cytokines, py_app$cytokines %>% data.table)

    app4 %>%
        lazy_prop(plot1_data, reactive({
            row <- req(this$input$cytokines_rows_selected)
            row <- this$cytokines[row, ]
            py_app$plot1_data(row$cytokine)
        }))

    app4 %>%
        lazy_prop(donors, reactive({
            x <- this$plot1_data()
            donors <- x[[2]] %>% data.table %>% 
                {.[, list(donor, group, resp)]} %>%
                unique(by=NULL) %>%
                {.[, cytokine:=x[[4]]]}
            donors[]
        }))

    app4 %>%
        lazy_prop(plot2_data, reactive({
            donors <- this$donors()
            row <- req(this$input$donors_rows_selected)
            row <- donors[row, ]
            data <- py_app$plot2_data(row$cytokine, row$donor)
            data[[3]] <- row$cytokine
            data[[4]] <- row$donor
            data
        }))

    app4 %>%
        lazy_prop(deg, {
            py_app$fit_level4_deg
        })

    app4 %>%
        lazy_prop(plot4_data, {
            x2 <- py_app$fit_level4_pca %>% data.table
            x4 <- x2$pc1 %>% range
            x5 <- x2$pc2 %>% range
            x3 <- x2 %>% split(x2$t) %>%
                lapply(function(x) {
                    ggplot(x)+
                        aes(
                            pc1, pc2,
                            color=dsm_severity_score_group
                        )+
                        geom_point()+
                        lims(x=x4, y=x5)
                })
        })

    app4$ui <- {
        fluidPage(
            fluidRow(
                column(6,
                    DTOutput('cytokines'),
                )
            ),
            fluidRow(
                column(3,
                    htmlOutput('plot1_stat1')
                ),
                column(3,
                    htmlOutput('plot1_stat2')
                ),
                column(3,
                    htmlOutput('plot1_stat3')
                ),
                column(3,
                    htmlOutput('plot1_stat4')
                )
            ),
            fluidRow(
                column(8,
                    plotOutput('plot1', height=400)
                )
            ),
            fluidRow(
                DTOutput('donors')
            ),
            fluidRow(
                column(8,
                    plotOutput('plot2', height=200)
                )
            ),
           fluidRow(
                DTOutput('deg')
            ),
            sliderInput('t', 'T:', min=1, max=50, value=0),
            fluidRow(
                column(8,
                    plotOutput('plot4', height=200)
                )
            )
        )
    }

    app4 %>%
        lazy_prop(server, function(input, output, session) {
            this$input <- input

            output$cytokines <- renderDT({
                this$cytokines %>%
                    datatable_scroller(
                        filter = 'top', selection = 'single',
                        options = list(dom = 't')
                    )
            })

            output$plot1_stat1 <- renderText({
                x <- this$plot1_data()[[1]]
                x[[1]] %>% py_to_r %>% kbl %>% kable_classic
            })

            output$plot1_stat2 <- renderText({
                x <- this$plot1_data()[[1]]
                x[[2]] %>% py_to_r %>% kbl %>% kable_classic
            })

            output$plot1_stat3 <- renderText({
                x <- this$plot1_data()[[1]]
                x[[3]] %>% py_to_r %>% kbl %>% kable_classic

            })

            output$plot1_stat4 <- renderText({
                x <- this$plot1_data()[[1]]
                sprintf('pval = %.2e', x[[4]])
            })

            output$plot1 <- renderPlot({
                x <- this$plot1_data()
                ggplot()+
                    geom_line(
                        data=x[[2]],
                        mapping=aes(t, level, group=donor, color=group),
                        alpha=0.5
                    )+
                    geom_line(
                        data=x[[3]],
                        mapping=aes(t, level, group=donor, color=group),
                        size=2
                    )+
                    facet_grid(resp~.)+
                    scale_color_manual(values=c('black', 'black', 'red', 'blue'))
            })

            output$donors <- renderDT({
                this$donors() %>%
                    datatable_scroller(
                        filter = 'top', selection = 'single',
                        options = list(dom = 't')
                    )
            })

            output$plot2 <- renderPlot({
                x <- this$plot2_data()
                x1 <- x[[1]]
                x2 <- x[[2]] %>% data.table
                d <- x[[4]]

                ggplot(x1)+
                    aes(days_since_onset, level)+
                    geom_point(aes(color=donor==d))+
                    geom_line(aes(group=donor, color=donor==d))+
                    geom_line(data=x2[variable=='pred1'], linetype='dotted')+
                    geom_line(data=x2[variable=='pred2'], linetype='dashed')+
                    geom_line(data=x2[variable=='pred3'], linetype='solid')+
                    theme(legend.position='none')+
                    scale_color_manual(values=c('lightgray', 'red'))
            })

            output$deg <- renderDT({
                this$deg %>%
                    datatable_scroller(
                        filter = 'top', selection = 'single',
                        options = list(dom = 't')
                    )
            })

            output$plot4 <- renderPlot({
                t <- req(input$t)
                this$plot4_data[[t]]
            })
        })


    app4 %>%
        lazy_prop(shiny, shinyApp(
            ui=this$ui,
            server=this$server
        ))
}

if (FALSE) {
    this <- app4

    this$input <- reactiveValues(
        cytokines_rows_selected=1,
        donors_rows_selected=1
    )
    this$input$donors_rows_selected <- 2
    x <- isolate(this$plot2_data())
    x[[4]]
}


runApp(app4$shiny)