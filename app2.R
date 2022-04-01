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
    )
  )
}

server <- function(input, output, session) {
    if (FALSE) {
        input <- list(cytokines=c('IL-6'))
    }

    py_run_string('from covid19_cell_atlas._app2 import app2 as app')
    py_app <- py_eval('app')

    output$cytokines <- renderDT(
        py_app$cytokines_table %>%
            datatable_scroller(
                filter = 'top', selection = 'multiple',
                options = list(dom = 't')
            )
    )

    output$plot1.ui <- renderUI({
        cytokines <- req(input$cytokines_rows_selected)
        plotOutput('plot1', height=200*length(cytokines))
    })

    output$plot1 <- renderPlot({
        row <- req(input$cytokines_rows_selected)
        cytokines <- py_app$cytokines_table$cytokine[row]
        data <- py_app$plot1_data(cytokines) %>% data.table
        data <- data[days_since_onset<=40]
        ggplot(data)+
            aes(x=days_since_onset, y=log2(level))+
            geom_line(aes(group=donor), alpha=0.1)+
            geom_point(aes(fill=dsm_severity_score_group), alpha=0.5, size=2, pch=21)+
            geom_smooth(aes(color=dsm_severity_score_group), alpha=0.5, formula='y~x', method='lm')+
            facet_grid(cytokine~., scales='free')+
            labs(y='log2(pg/mL)')
    })
}

app <- shinyApp(
    ui=ui,
    server=server
)

runApp(app)
