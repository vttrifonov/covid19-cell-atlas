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
        selectizeInput(
            'cytokines',
            'Select cytokines',
            choices=NULL,
            multiple=TRUE
        )
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
    {
        if (FALSE) {
            input <- list(cytokines=c('IL-6'), genes=c('IL6'))
        }

        py_run_string('from covid19_cell_atlas._analysis3 import analysis3 as analysis')
        analysis <- py_eval('analysis')

        round_mantisa <- function(x, n=1) {
            x1 <- floor(log10(x))
            x2 <- ceiling(x*10^(n-x1))/10^n
            x2*10^x1
        }

        genes <- py_eval('analysis.fit1.to_dataframe().reset_index()') %>%
                py_to_r() %>%
                data.table
        genes <- genes[!is.na(`Pr(>F)`)]
        genes <- genes[, list(
            subset, gene,
            `Pr(>F)`, q, F,
            Intercept,
            days_since_onset,
            `dsm_severity_score_group[T.DSM_low]`,
            `days_since_onset:dsm_severity_score_group[T.DSM_low]` 
        )]
        genes[, `Pr(>F)` := round_mantisa(`Pr(>F)`, 1)]
        genes[, q := round_mantisa(q, 1)]
        for(c in setdiff(names(genes), c('q', 'Pr(>F)', 'subset', 'gene'))) {
            genes[[c]] <- round(genes[[c]], 2)
        }

        plot1 <- reactive({
            cytokines <- req(input$cytokines)
            data <- analysis$app1_plot1_data(cytokines) %>% data.table
            data <- data[days_since_onset<=40]
            ggplot(data)+
                aes(x=days_since_onset, y=log2(level))+
                geom_line(aes(group=donor), alpha=0.1)+
                geom_point(aes(fill=dsm_severity_score_group), alpha=0.5, size=2, pch=21)+
                geom_smooth(aes(color=dsm_severity_score_group), alpha=0.5, formula='y~x', method='lm')+
                facet_grid(cytokine~., scales='free')+
                labs(y='log2(pg/mL)')
        })


        plot2 <- reactive({
            rows <- req(input$genes_rows_selected)
            gene <- genes$gene[rows]
            data <- analysis$app1_plot2_data(gene) %>% data.table

            ggplot(data)+
                aes(days_since_onset, y=log2(X+1))+
                geom_point(aes(fill=dsm_severity_score_group), alpha=0.5, size=2, pch=21)+
                geom_line(aes(group=donor), alpha=0.1)+
                geom_smooth(method='lm', aes(color=dsm_severity_score_group), alpha=0.5)+
                facet_grid(subset+gene~., scales='free')+
                labs(y='log2RPM')
        })

        plot3 <- reactive({
            row <- req(input$genes_rows_selected)
            gene <- genes$gene[row]
            data <- analysis$app1_plot3_data(gene) %>% data.table

            ggplot(data)+
                aes(days_since_onset, y=log(X+1))+
                geom_point(aes(fill=dsm_severity_score_group), alpha=0.5)+
                geom_line(aes(group=donor), alpha=0.1)+
                geom_smooth(method='lm', aes(color=dsm_severity_score_group), alpha=0.5)+
                facet_grid(subset+cell_type+gene~., scales='free_x')+
                labs(y='log2RPM')
        })
    }

    {
        updateSelectizeInput(
            session, 'cytokines',
            choices = analysis$app1_cytokines,
            server = TRUE
        )

        output$plot1.ui <- renderUI({
            cytokines <- req(input$cytokines)
            plotOutput('plot1', height=200*length(cytokines))
        })

        output$plot1 <- renderPlot(plot1())

        output$genes <- renderDT(
            genes %>%
                datatable_scroller(
                    filter = 'top', selection = 'single',
                    options = list(dom = 't')
                )
        )

        output$plot2 <- renderPlot(plot2())
        output$plot3 <- renderPlot(plot3())
    }
}

app <- shinyApp(
    ui=ui,
    server=server
)

runApp(app)
