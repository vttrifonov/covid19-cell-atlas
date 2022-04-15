{
    library(data.table)
    library(magrittr)
    library(shiny)
    library(ggplot2)
    library(DT)
    library(reticulate)
    library(glue)
    source('helpers.R')
    source('common/module.R')
    load.module('common')
    py_run_string('from covid19_cell_atlas._app3 import app3 as app')
    py_app <- py_eval('app')
}

{
    app3 <- obj()

    app3 %>%
        lazy_prop(fit_level2, py_app$fit_level2 %>% data.table)

    app3 %>%
        lazy_prop(fit_level2_row, reactive({
            row <- req(this$input$fit_level2_rows_selected)
            row <- this$fit_level2[row, ]
            row
        }))

    app3 %>%
        lazy_prop(plot1_data, reactive({
            row <- this$fit_level2_row()
            lapply(py_app$plot1_data(row$cytokine, row$donor), data.table)
        }))

    app3 %>%
        lazy_prop(plot2_data, {
            x2 <- py_app$data2_pca %>% data.table
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

    app3$ui <- {
        fluidPage(
            DTOutput('fit_level2'),
            fluidRow(
                column(8,
                    plotOutput('plot1', height=200)
                )
            ),
            selectInput(
                'cytokine', 'Cytokine:',
                choices=py_app$plot3_cytokines
            ),
            fluidRow(
                column(8,
                    plotOutput('plot3', height=200)
                )
            ),
            sliderInput('t', 'T:', min=0, max=28, value=0),
            fluidRow(
                column(8,
                    plotOutput('plot2', height=200)
                )
            )
        )
    }

    app3 %>%
        lazy_prop(server, function(input, output, session) {
            this$input <- input

            output$fit_level2 <- renderDT({
                this$fit_level2 %>%
                    datatable_scroller(
                        filter = 'top', selection = 'single',
                        options = list(dom = 't')
                    )
            })

            output$plot1 <- renderPlot({
                row <- this$fit_level2_row()
                x <- this$plot1_data()
                x1 <- x[[1]] %>% data.table
                x2 <- x[[2]] %>% data.table
                x1 <- x1[order(donor==row$donor)]
                ggplot(x1)+
                    aes(days_since_onset, level)+
                    geom_point()+
                    geom_line(data=x2[variable=='pred1'], linetype='dotted', color='red')+
                    geom_line(data=x2[variable=='pred2'], linetype='dashed', color='red')+
                    geom_line(data=x2[variable=='pred3'], linetype='solid')+
                    geom_line(aes(group=donor, color=donor==row$donor))+
                    theme(legend.position='none')+
                    scale_color_manual(values=c('gray', 'red'))
            })

            output$plot2 <- renderPlot({
                t <- req(input$t) %>% as.character
                this$plot2_data[[t]]
            })

            output$plot3 <- renderPlot({
                cytokine <- req(input$cytokine)
                ggplot(py_app$plot3_data(cytokine))+
                    aes(t, level, color=dsm_severity_score_group)+
                    geom_line(aes(group=donor))+
                    facet_wrap(~dsm_severity_score_group)
            })
        })


    app3 %>%
        lazy_prop(shiny, shinyApp(
            ui=this$ui,
            server=this$server
        ))
}

test1 <- function() {
    this <- app3

    this$input <- reactiveValues(
        cytokine='IL-6'
    )
    c <- isolate(req(this$input$cytokine))
    x <- isolate(py_app$plot3_data(c))
}


runApp(app3$shiny)