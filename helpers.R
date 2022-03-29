  datatable_scroller<-function(data, options=list(), extensions=c(), ...)  {
    default<-list(
      scrollY=200
    )
    required<-list(
      scrollX=TRUE,
      deferRender=TRUE,
      scroller=TRUE
    )
    
    for(name in names(required)) options[[name]]<-required[[name]]
    for(name in names(default)) 
      if (!name %in% names(options)) options[[name]]<-default[[name]]
    
    DT::datatable(data, 
                  rownames=FALSE,
                  extensions = c(extensions, 'Scroller'),
                  options = options,
                  ...
    )
  }