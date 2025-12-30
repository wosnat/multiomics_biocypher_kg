datacache <- new.env(hash=TRUE, parent=emptyenv())

org.SspWH8102.eg <- function() showQCData("org.SspWH8102.eg", datacache)
org.SspWH8102.eg_dbconn <- function() dbconn(datacache)
org.SspWH8102.eg_dbfile <- function() dbfile(datacache)
org.SspWH8102.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.SspWH8102.eg_dbInfo <- function() dbInfo(datacache)

org.SspWH8102.egORGANISM <- "Synechococcus spWH8102"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.SspWH8102.eg.sqlite", package=pkgname, lib.loc=libname)
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
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.SspWH8102.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.SspWH8102.eg_dbconn())
}

