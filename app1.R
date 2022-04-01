library(data.table)
library(magrittr)
library(shiny)
library(ggplot2)
library(DT)
library(reticulate)
source('helpers.R')

ui <- {
  fluidPage(
    fluidRow(
        DTOutput('cytokines')
    ),
    fluidRow(
        column(7,
            uiOutput('plot1.ui')
        )
    ),
    fluidRow(
        DTOutput('genes')
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

server <- function(input, output, session) {
    if (FALSE) {
        input <- list(cytokines=c('IL-6'), genes=c('IL6'))
    }

    py_run_string('from covid19_cell_atlas._app1 import app1 as app')
    py_app <- py_eval('app')

    output$genes <- renderDT(
        py_app$genes_table %>%
            datatable_scroller(
                filter = 'top', selection = 'single',
                options = list(dom = 't')
            )
    )

    output$plot2 <- renderPlot({
        rows <- req(input$genes_rows_selected)
        gene <- py_app$genes_table$gene[rows]
        data <- py_app$plot2_data(gene) %>% data.table

        ggplot(data)+
            aes(days_since_onset, y=log2(X+1))+
            geom_point(aes(fill=dsm_severity_score_group), alpha=0.5, size=2, pch=21)+
            geom_line(aes(group=donor), alpha=0.1)+
            geom_smooth(method='lm', aes(color=dsm_severity_score_group), alpha=0.5)+
            facet_grid(subset+gene~., scales='free')+
            labs(y='log2RPM')
    })

    output$plot3 <- renderPlot({
        row <- req(input$genes_rows_selected)
        gene <- py_app$genes_table$gene[row]
        data <- py_app$plot3_data(gene) %>% data.table

        ggplot(data)+
            aes(days_since_onset, y=log(X+1))+
            geom_point(aes(fill=dsm_severity_score_group), alpha=0.5)+
            geom_line(aes(group=donor), alpha=0.1)+
            geom_smooth(method='lm', aes(color=dsm_severity_score_group), alpha=0.5)+
            facet_grid(subset+cell_type+gene~., scales='free_x')+
            labs(y='log2RPM')
    })
}

app <- shinyApp(
    ui=ui,
    server=server
)

runApp(app)
