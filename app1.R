{
    library(data.table)
    library(magrittr)
    library(shiny)
    library(ggplot2)
    library(DT)
    library(reticulate)
    source('helpers.R')
    source('common/module.R')
    load.module('common')
    py_run_string('from covid19_cell_atlas._app1 import app1 as app')
    py_app <- py_eval('app')
}

{
    app1 <- obj()

    app1 %>%
        lazy_prop(enrich_table, reactive({
            py_app$enrich2_table %>% data.table
        }))

    app1 %>%
        lazy_prop(enrich_table_selected, reactive({
            table <- this$enrich_table()
            rows <- req(this$input$enrich_rows_selected)
            table[rows, ]
        }))

    app1 %>%
        lazy_prop(enrich, reactive({
            this$enrich_table() %>%
                datatable_scroller(
                    filter = 'top', selection = 'single',
                    options = list(dom = 't')
                )
        }))

    app1 %>%
        lazy_prop(genes_table, reactive({
            .enrich_table <- this$enrich_table_selected()
            sigs <- py_app$enrich2_leading_edge %>% data.table
            sigs <- sigs[sig==.enrich_table$sig]
            sigs <- sigs[subset==.enrich_table$subset]
            genes <- py_app$genes_table %>% data.table
            genes <- genes[subset==.enrich_table$subset]
            genes <- genes[gene %in% sigs$gene]
            genes
        }))

    app1 %>%
        lazy_prop(genes_table_selected, reactive({
            table <- this$genes_table()
            rows <- req(this$input$genes_rows_selected)
            table[rows, ]
        }))

    app1 %>%
        lazy_prop(genes, reactive({
            this$genes_table() %>%
                datatable_scroller(
                    filter = 'top', selection = 'single',
                    options = list(dom = 't')
                )
        }))

    app1 %>%
        lazy_prop(plot2, reactive({
            .genes_table <- this$genes_table_selected()
            gene <- .genes_table$gene
            data <- py_app$plot2_data(gene) %>% data.table

            ggplot(data)+
                aes(days_since_onset, y=log2(X+1))+
                geom_point(aes(fill=dsm_severity_score_group), alpha=0.5, size=2, pch=21)+
                geom_line(aes(group=donor), alpha=0.1)+
                geom_smooth(method='lm', aes(color=dsm_severity_score_group), alpha=0.5)+
                facet_grid(subset+gene~., scales='free')+
                labs(y='log2RPM')
        }))

    app1 %>%
        lazy_prop(plot3, reactive({
            .genes_table <- this$genes_table_selected()
            gene <- .genes_table$gene
            data <- py_app$plot3_data(gene) %>% data.table

            ggplot(data)+
                aes(days_since_onset, y=log(X+1))+
                geom_point(aes(fill=dsm_severity_score_group), alpha=0.5)+
                geom_line(aes(group=donor), alpha=0.1)+
                geom_smooth(method='lm', aes(color=dsm_severity_score_group), alpha=0.5)+
                facet_grid(subset+cell_type+gene~., scales='free_x')+
                labs(y='log2RPM')
        }))

    app1$ui <- {
    fluidPage(
        fluidRow(
            column(6,
                DTOutput('enrich')
            ),
            column(6,
                DTOutput('genes')
            )
        ),
        fluidRow(
            column(8,
                plotOutput('plot2', height=200*2)
            )
        ),
        fluidRow(
            column(8,
                plotOutput('plot3', height=200*24)
            )
        )
    )
    }

    app1 %>%
        lazy_prop(server, function(input, output, session) {
            this$input <- input
            output$enrich <- renderDT(this$enrich())
            output$genes <- renderDT(this$genes())
            output$plot2 <- renderPlot(this$plot2())
            output$plot3 <- renderPlot(this$plot3())
        })

    app1 %>%
        lazy_prop(shiny, shinyApp(
            ui=this$ui,
            server=this$server
        ))
}

test1 <- function() {
    this <- app1

    this$input <- reactiveValues(
        enrich_rows_selected=1,
        genes_rows_selected=1
    )

    isolate(this$plot2())
}

runApp(app1$shiny)
