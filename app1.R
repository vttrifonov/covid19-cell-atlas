
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


app1 <- obj()

app1 %>%
    lazy_prop(enrich, {
        py_app$enrich1_table %>%
            datatable_scroller(
                filter = 'top', selection = 'single',
                options = list(dom = 't')
            )
    })

app1 %>%
    lazy_prop(genes, {
        py_app$genes_table %>%
            datatable_scroller(
                filter = 'top', selection = 'single',
                options = list(dom = 't')
            )
    })

 app1 %>%
    lazy_prop(plot2, reactive({
        .genes_table <- this$genes_table()
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
        .genes_table <- this$genes_table()
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
        this$genes_table <- reactive({
            rows <- req(input$genes_rows_selected)
            py_app$genes_table[rows,]
        })

        output$enrich <- renderDT(this$enrich)

        output$genes <- renderDT(this$genes)

        output$plot2 <- renderPlot(this$plot2())

        output$plot3 <- renderPlot(this$plot3())
    })

app1 %>%
    lazy_prop(shiny, shinyApp(
        ui=this$ui,
        server=this$server
    ))

test1 <- function() {
    x <- py_app$genes_table %>% data.table
    x <- x[gene=='IL6' & subset=='innate']
    app1$genes_table <- function() x

    isolate(app1$plot2())
}

runApp(app1$shiny)
