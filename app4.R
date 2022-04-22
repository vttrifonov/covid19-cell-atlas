{
    library(data.table)
    library(magrittr)
    library(shiny)
    library(ggplot2)
    library(DT)
    library(reticulate)
    library(glue)
    library(kableExtra)
    library(gridExtra)
    library(cowplot)
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
        lazy_prop(plot3_data, reactive({
            row <- req(this$input$deg_rows_selected)
            row <-  this$deg[row, ]
            data <- py_app$plot3_data(row$cytokine, row$t)
            data
        }))

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
            fluidRow(
                column(8,
                    plotOutput('plot3', height=300)
                )
            ),

            sliderInput(
                'pca_n', 'PCA number of top cytokines',
                min=2, max=nrow(app4$cytokines), value=nrow(app4$cytokines)
            ),
            sliderInput('pca_t', 'PCA timepoint', min=2, max=37, value=2),
            fluidRow(
                column(8,
                    uiOutput('plot4_container')
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

            output$plot3 <- renderPlot({
                x <- this$plot3_data()
                ggplot(x)+
                    aes(dsm_severity_score_group, level)+
                    geom_boxplot()+
                    geom_point()
            })

            output$plot4_container <- renderUI({
                 n <- req(this$input$pca_n) %>% as.integer
                 plotOutput('plot4', height=200+100+20*(n+1.5))
            })

            output$plot4 <- renderPlot({
                t <- req(this$input$pca_t)
                n <- req(this$input$pca_n) %>% as.integer

                x1 <- py_app$plot4_data(n)

                .get <- function(l, t) {
                    t1 <- names(l) %>% as.numeric
                    i <- (t1-t) %>% abs %>% which.min
                    l[[i]]
                }

                x2 <- x1[[1]]
                x3 <- .get(x1[[2]], t) %>% py_to_r
                x4 <- .get(x1[[3]], t) %>% py_to_r

                x5 <- .get(x2[[1]], t) %>% py_to_r
                x6 <- x2[[2]]
                x7 <- x2[[3]]

                x8 <- ggplot(x5)+
                    aes(
                        pc1, pc2,
                        color=dsm_severity_score_group
                    )+
                    geom_point()+
                    lims(x=x6, y=x7)

                clip <- function(x, l, h)
                    ifelse(x<l, l, ifelse(x>h, h, x))
                x9 <- ggplot(x3)+
                    aes(
                        donor, cytokine,
                        fill=clip(level, -5, 5)
                    )+
                    geom_tile()+
                    scale_fill_gradient2(
                        low='blue', mid='white', high='red', midpoint=0
                    )+
                    theme(
                        axis.text.x=element_blank(),
                        legend.position='right'
                    )+
                    labs(x='', y='')

                x10 <- ggplot(x4)+
                    aes(donor, y=1, fill=dsm_severity_score_group)+
                    geom_tile()+
                    theme(
                        axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        legend.position='none'
                    )+
                    labs(y='')

                plot_grid(
                    x8,
                    plot_grid(
                        x9, x10, axis='lr', align='v',
                        ncol=1, rel_heights=c(n, 1.5)
                    ),
                    ncol=1, rel_heights=c(200, 100+20*(n+1.5))
                )
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
        donors_rows_selected=1,
        pca_t = 10.0,
        pca_n = 10
    )

    t <- isolate(this$input$pca_t)
    n <- isolate(this$input$pca_n) %>% as.integer
}


runApp(app4$shiny)