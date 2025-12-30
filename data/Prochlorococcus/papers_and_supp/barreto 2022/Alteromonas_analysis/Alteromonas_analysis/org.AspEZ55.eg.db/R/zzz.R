datacache <- new.env(hash=TRUE, parent=emptyenv())

org.AspEZ55.eg <- function() showQCData("org.AspEZ55.eg", datacache)
org.AspEZ55.eg_dbconn <- function() dbconn(datacache)
org.AspEZ55.eg_dbfile <- function() dbfile(datacache)
org.AspEZ55.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.AspEZ55.eg_dbInfo <- function() dbInfo(datacache)

org.AspEZ55.egORGANISM <- "Alteromonas spEZ55"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.AspEZ55.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.AspEZ55.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.AspEZ55.eg_dbconn())
}

