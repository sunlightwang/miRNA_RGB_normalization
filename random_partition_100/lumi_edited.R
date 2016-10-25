lumiR <- function (fileName, sep = NULL, detectionTh = 0.01, na.rm = TRUE, 
    convertNuID = TRUE, lib.mapping = NULL, dec = ".", parseColumnName = FALSE, 
    checkDupId = TRUE, QC = TRUE, columnNameGrepPattern = list(exprs = "AVG_SIGNAL", 
        se.exprs = "BEAD_STD", detection = "DETECTION", beadNum = "Avg_NBEADS"), 
    inputAnnotation = TRUE, annotationColumn = c("ACCESSION", 
        "SYMBOL", "PROBE_SEQUENCE", "PROBE_START", "CHROMOSOME", 
        "PROBE_CHR_ORIENTATION", "PROBE_COORDINATES", "DEFINITION"), 
    verbose = TRUE, ...) 
{
    if (is.null(columnNameGrepPattern$exprs)) 
        columnNameGrepPattern$exprs <- "AVG_SIGNAL"
    if (is.null(columnNameGrepPattern$se.exprs)) 
        columnNameGrepPattern$se.exprs <- "BEAD_STD"
    if (is.null(columnNameGrepPattern$detection)) 
        columnNameGrepPattern$detection <- "DETECTION"
    if (is.null(columnNameGrepPattern$beadNum)) 
        columnNameGrepPattern$beadNum <- "Avg_NBEADS"
    if (is.na(columnNameGrepPattern$exprs)) {
        columnNameGrepPattern$exprs <- "AVG_SIGNAL"
        warning("exprs slot is required and default pattern will be used!\n")
    }
    if (is.na(columnNameGrepPattern$se.exprs) && checkDupId) {
        warning("se.exprs slot is required for the VST transformation!\n We strongly suggest to include BEAD_STD columns!\n")
    }
    if (is.na(columnNameGrepPattern$beadNum)) {
    }
    history.submitted <- as.character(Sys.time())
    oldSetting <- options("stringsAsFactors")[[1]]
    options(stringsAsFactors = FALSE)
    version <- 2
    if (!file.exists(fileName)) 
        stop("The file does not exist! Please check your file path!\n")
    info <- readLines(fileName, n = 20)
    nMetaDataLines <- grep(columnNameGrepPattern$exprs, info, 
        ignore.case = TRUE) - 1
    if (is.null(sep)) {
        titleLine <- info[nMetaDataLines + 1]
        dataLine1 <- info[nMetaDataLines + 2]
        dataLine2 <- info[nMetaDataLines + 3]
        sepNum1 <- gregexpr("\t", dataLine1)[[1]]
        sepNum2 <- gregexpr("\t", dataLine2)[[1]]
        if (sepNum1[1] > 0 && length(sepNum1) == length(sepNum2)) {
            sep <- "\t"
        }
        else if (dec != ",") {
            sepNum1 <- gregexpr(",", dataLine1)[[1]]
            sepNum2 <- gregexpr(",", dataLine2)[[1]]
            if (sepNum1[1] > 0 && length(sepNum1) == length(sepNum2)) {
                sep <- ","
            }
            else {
                stop("The seperator is not Tab or comma!\n Please specify the seperator used in the file!\n")
            }
        }
        else {
            stop("Please specify the seperator used in the file!\n")
        }
    }
    dataLine1 <- strsplit(info[nMetaDataLines + 2], sep)[[1]]
    quoteCount1 <- gregexpr("\"", dataLine1[1])[[1]]
    quoteCount2 <- gregexpr("'", dataLine1[1])[[1]]
    if (length(quoteCount1) == 2) {
        quote <- "\""
    }
    else if (length(quoteCount2) == 2) {
        quote <- "'"
    }
    else {
        quote <- ""
    }
    header <- strsplit(info[nMetaDataLines + 1], sep)[[1]]
    probeId.pos <- grep("Probe.?ID", header, ignore.case = TRUE)
    if (length(probeId.pos) > 0) {
        colClasses <- rep(NA, length(header))
        colClasses[probeId.pos] <- "character"
    }
    else {
        colClasses <- NA
    }
    if (nMetaDataLines > 0) {
        info <- readLines(fileName, n = nMetaDataLines)
        markerInd <- grep("^\\[.*\\]", info, ignore.case = TRUE)
        if (length(markerInd) > 0) {
            if (length(grep("^\\[Header\\]", info[markerInd[1]], 
                ignore.case = TRUE)) == 0) 
                warning("The data file may not be in the Illumina BeadStudio or GenomeStudio output format!\n")
            if (length(markerInd) > 1) {
                if (length(grep("^\\[.*\\Profile]", info[markerInd[2]], 
                  ignore.case = TRUE)) == 0) 
                  warning("The data file may not be in the Illumina BeadStudio or GenomeStudio output format!\n")
            }
            version <- 3
            info <- info[-markerInd]
        }
        if (length(info) > 0) {
            info <- sub("[[:blank:]]+$", "", info)
            info <- sub(paste(sep, "+$", sep = ""), "", info)
            ind <- grep("Normalization", info, ignore.case = TRUE)
            if (length(ind) > 0) {
                if (version == 2) {
                  normalization <- strsplit(info, split = "=")[[ind]][2]
                  normalization <- gsub(pattern = " |,", replacement = "", 
                    normalization)
                }
                else {
                  normalization <- strsplit(info, split = sep)[[ind]][2]
                }
                if (length(grep("none", normalization, ignore.case = TRUE)) == 
                  0) {
                  warning("We recommend the raw data not to be normalized in BeadStudio or GenomeStudio.\n")
                }
            }
        }
        else {
            info <- NULL
        }
    }
    else {
        info <- NULL
    }
    allData <- read.table(file = fileName, header = TRUE, sep = sep, 
        dec = dec, skip = nMetaDataLines, row.names = NULL, colClasses = colClasses, 
        quote = quote, as.is = TRUE, check.names = FALSE, strip.white = TRUE, 
        comment.char = "", fill = TRUE, ...)
    sectionInd <- grep("^\\[.*\\]", allData[, 1], ignore.case = TRUE)
    if (length(sectionInd) > 0) {
        if (is.na(version)) 
            version <- 3
        otherData <- allData[sectionInd[1]:nrow(allData), ]
        allData <- allData[1:(sectionInd[1] - 1), , drop = FALSE]
        naCol <- apply(allData, 2, function(x) all(is.na(x) | 
            x == ""))
        allData <- allData[, !naCol]
        sectionInd <- sectionInd - sectionInd[1] + 1
        sectionName <- otherData[sectionInd, 1]
        controlInd <- grep("^\\[Control.*\\]", sectionName, ignore.case = TRUE)
        if (length(controlInd) > 1) {
            ind <- grep("^\\[Control probe.*\\]", sectionName[controlInd], 
                ignore.case = TRUE)
            if (length(ind) > 0) {
                controlInd <- controlInd[ind[1]]
            }
            else {
                controlInd <- controlInd[1]
            }
        }
        if (length(controlInd) > 0) {
            startRow <- sectionInd[controlInd] + 1
            if (length(sectionInd) > controlInd) 
                endRow <- sectionInd[controlInd + 1] - 1
            else endRow <- nrow(otherData)
            controlData <- otherData[startRow:endRow, ]
            naCol <- apply(controlData, 2, function(x) all(is.na(x) | 
                x == ""))
            controlData <- controlData[, !naCol, drop = FALSE]
            tmpFile <- tempfile(pattern = "file", tmpdir = tempdir())
            write.table(controlData, tmpFile, sep = "\t", col.names = FALSE, 
                row.names = FALSE)
            controlData <- getControlData(tmpFile, type = "data.frame")
            unlink(tmpFile)
        }
        else {
            controlData <- data.frame()
        }
        summaryInd <- grep("^\\[Sample.*Table\\]", sectionName, 
            ignore.case = TRUE)
        if (length(summaryInd) > 0) {
            startRow <- sectionInd[summaryInd] + 1
            if (length(sectionInd) > summaryInd) 
                endRow <- sectionInd[summaryInd + 1] - 1
            else endRow <- nrow(otherData)
            sampleSummary <- otherData[startRow:endRow, ]
            naCol <- apply(sampleSummary, 2, function(x) all(is.na(x) | 
                x == ""))
            sampleSummary <- sampleSummary[, !naCol]
            colnames(sampleSummary) <- sampleSummary[1, ]
            sampleSummary <- sampleSummary[-1, ]
        }
        else {
            sampleSummary <- data.frame()
        }
    }
    else {
        controlData <- sampleSummary <- data.frame()
    }
    header <- names(allData)
    targetID <- as.character(as.vector(allData[, 1]))
    id <- targetID
    idName <- header[1]
    if (length(grep("Probe.?ID", header[2], ignore.case = TRUE)) > 
        0) {
        id <- as.character(as.vector(allData[, 2]))
        idName <- header[2]
    }
    else if (!is.null(lib.mapping)) {
        probeId.pos <- grep("Probe.?ID", header, ignore.case = TRUE)
        if (length(probeId.pos) > 0) {
            id <- as.character(as.vector(allData[, probeId.pos[1]]))
            idName <- header[probeId.pos[1]]
        }
    }
    ind <- grep(columnNameGrepPattern$exprs, header, ignore.case = TRUE)
    if (length(ind) == 0) {
        stop("Input data format unrecognizable!\nThere is no column name contains \"AVG_SIGNAL\"!\n")
    }
    else {
        ind2 <- grep(paste("\\(", columnNameGrepPattern$exprs, 
            "\\)", sep = ""), header, ignore.case = TRUE)
        if (length(ind2) == length(ind)/2) 
            ind <- ind[!(ind %in% ind2)]
    }
    exprs <- as.matrix(allData[, ind])
    if (!is.double(exprs[1])) {
        exprs <- matrix(as.double(exprs), nrow = nrow(allData))
    }
    colnames(exprs) <- header[ind]
    if (is.na(columnNameGrepPattern$se.exprs)) {
        ind <- NULL
    }
    else {
        ind <- grep(columnNameGrepPattern$se.exprs, header, ignore.case = TRUE)
    }
    if (length(ind) == 0) {
        se.exprs <- NULL
    }
    else {
        se.exprs <- as.matrix(allData[, ind])
        if (!is.double(se.exprs[1])) {
            se.exprs <- matrix(as.double(se.exprs), nrow = nrow(allData))
        }
        colnames(se.exprs) <- header[ind]
    }
    if (is.na(columnNameGrepPattern$detection)) {
        ind <- NULL
    }
    else {
        ind <- grep(columnNameGrepPattern$detection, header, 
            ignore.case = TRUE)
    }
    if (length(ind) == 0) {
        detection <- NULL
    }
    else {
        detection <- as.matrix(allData[, ind])
        if (!is.double(detection[1])) {
            detection <- matrix(as.double(detection), nrow = nrow(allData))
        }
        colnames(detection) <- header[ind]
        if (length(grep("Detection Pval", header, ignore.case = TRUE)) == 
            0) {
            detection <- 1 - detection
        }
    }
    if (is.na(columnNameGrepPattern$beadNum)) {
        ind <- NULL
    }
    else {
        ind <- grep(columnNameGrepPattern$beadNum, header, ignore.case = TRUE)
    }
    if (length(ind) == 0) {
        beadNum <- NULL
    }
    else {
        beadNum <- as.matrix(allData[, ind])
        if (!is.double(beadNum[1])) {
            beadNum <- matrix(as.double(beadNum), nrow = nrow(allData))
        }
        colnames(beadNum) <- header[ind]
    }
    if (!is.null(beadNum) && length(grep("BEAD_STDERR", colnames(se.exprs), 
        ignore.case = T)) == ncol(se.exprs)) {
        se.exprs <- se.exprs * sqrt(beadNum)
    }
    annotationInfo <- NULL
    if (inputAnnotation) {
        annotationColumn <- header[toupper(header) %in% toupper(annotationColumn)]
        if (length(annotationColumn) == 0) {
            cat("Annotation columns are not available in the data.\n")
        }
        else {
            annotationInfo <- allData[, annotationColumn, drop = FALSE]
        }
    }
    if (checkDupId) {
        dupId <- unique(id[duplicated(id)])
        if (length(dupId) > 0) {
            cat("Duplicated IDs found and were merged!\n")
            rmInd <- NULL
            for (dupId.i in dupId) {
                selInd.i <- which(id == dupId.i)
                exprs[selInd.i[1], ] <- colMeans(exprs[selInd.i, 
                  , drop = FALSE])
                if (is.null(beadNum)) {
                  if (!is.null(se.exprs)) 
                    se.exprs[selInd.i[1], ] <- colMeans(se.exprs[selInd.i, 
                      , drop = FALSE])
                }
                else {
                  totalBead.i <- colSums(beadNum[selInd.i, , 
                    drop = FALSE])
                  beadNum[selInd.i[1], ] <- totalBead.i
                  if (!is.null(se.exprs)) {
                    temp <- colSums(se.exprs[selInd.i, , drop = FALSE]^2 * 
                      (beadNum[selInd.i, , drop = FALSE] - 1))
                    temp <- temp/(totalBead.i - length(selInd.i))
                    se.exprs[selInd.i[1], ] <- sqrt(temp * (colSums(1/beadNum[selInd.i, 
                      , drop = FALSE])))
                  }
                }
                if (!is.null(detection)) {
                  detection[selInd.i[1], ] <- apply(detection[selInd.i, 
                    , drop = FALSE], 2, max)
                }
                rmInd <- c(rmInd, selInd.i[-1])
            }
            exprs <- exprs[-rmInd, , drop = FALSE]
            if (!is.null(se.exprs)) 
                se.exprs <- se.exprs[-rmInd, , drop = FALSE]
            id <- id[-rmInd]
            targetID <- targetID[-rmInd]
            if (!is.null(detection)) 
                detection <- detection[-rmInd, , drop = FALSE]
            if (!is.null(beadNum)) 
                beadNum <- beadNum[-rmInd, , drop = FALSE]
            if (!is.null(annotationInfo)) 
                annotationInfo <- annotationInfo[-rmInd, , drop = FALSE]
        }
    }
    if (na.rm) {
        keepInd <- apply(is.na(exprs), 1, sum) == 0
        exprs <- exprs[keepInd, , drop = FALSE]
        se.exprs <- se.exprs[keepInd, , drop = FALSE]
        if (!is.null(beadNum)) 
            beadNum <- beadNum[keepInd, , drop = FALSE]
        if (!is.null(detection)) 
            detection <- detection[keepInd, , drop = FALSE]
        if (!is.null(annotationInfo)) 
            annotationInfo <- annotationInfo[keepInd, , drop = FALSE]
        id <- id[keepInd]
        targetID <- targetID[keepInd]
    }
    if (!any(duplicated(id))) {
        rownames(exprs) <- id
        if (!is.null(se.exprs)) 
            rownames(se.exprs) <- id
        if (!is.null(beadNum)) 
            rownames(beadNum) <- id
        if (!is.null(detection)) 
            rownames(detection) <- id
    }
    pattern <- paste("[\\.\\:][^[:alnum:]]*", columnNameGrepPattern$exprs, 
        "[^[:alnum:]]*", sep = "")
    if (length(grep(pattern, colnames(exprs), ignore.case = TRUE)) == 
        0) 
        pattern <- paste("[^[:alnum:]]*", columnNameGrepPattern$exprs, 
            "[^[:alnum:]]*", sep = "")
    sampleID <- sub(pattern, "", colnames(exprs), ignore.case = TRUE)
    if (any(duplicated(sampleID))) {
        warning("Duplicated column names found in the raw data! \n A suffix number is added to the duplicated column names.\n")
        dupId <- which(duplicated(sampleID))
        dupName <- unique(sampleID[dupId])
        for (dupName.i in dupName) {
            dupInd.i <- which(sampleID == dupName.i)
            sampleID[dupInd.i] <- paste(sampleID[dupInd.i], 1:length(dupInd.i), 
                sep = ".")
        }
    }
    if (parseColumnName) {
        sampleIDInfo <- strsplit(sampleID, split = "_")
        label <- NULL
        newID <- NULL
        temp <- lapply(sampleIDInfo, function(x) {
            label <<- c(label, x[length(x)])
            newID <<- c(newID, paste(x[1:2], collapse = "_"))
        })
        if (!any(duplicated(newID))) 
            sampleID <- newID
        if (length(unique(label)) != length(label) || length(label) == 
            0 || any(is.na(label))) 
            label <- sampleID
    }
    else {
        label <- sampleID
    }
    reporterInfo <- data.frame(id)
    names(reporterInfo) <- idName
    varMetadata <- data.frame(labelDescription = "The Illumina microarray identifier")
    varName <- idName
    if (idName != header[1]) {
        reporterInfo <- data.frame(reporterInfo, TargetID = targetID)
        varMetadata <- data.frame(labelDescription = c(varMetadata$labelDescription, 
            "The Illumina TargetID"))
        varName <- c(varName, "TargetID")
    }
    if (!is.null(annotationInfo)) {
        reporterInfo <- data.frame(reporterInfo, annotationInfo)
        varMetadata <- data.frame(labelDescription = c(varMetadata$labelDescription, 
            names(annotationInfo)))
        varName <- c(varName, names(annotationInfo))
    }
    rownames(varMetadata) <- make.names(varName, unique = T)
    if (!any(duplicated(id))) 
        rownames(reporterInfo) <- id
    featureData <- new("AnnotatedDataFrame", data = reporterInfo, 
        varMetadata = varMetadata)
    colnames(exprs) <- label
    if (!is.null(se.exprs)) {
        if (ncol(exprs) != ncol(se.exprs)) 
            stop("Different column numbers of exprs and se.exprs! Please check the input data format.\n")
        colnames(se.exprs) <- label
    }
    if (!is.null(beadNum)) {
        if (ncol(beadNum) == length(label)) {
            colnames(beadNum) <- label
        }
        else {
            warning("The number of beadNum columns does not match! Please check the input data format.\n")
            if (ncol(beadNum) > length(label)) 
                beadNum <- beadNum[, 1:length(label)]
            if (ncol(beadNum) < length(label)) {
                for (i in 1:(length(label) - ncol(beadNum))) beadNum <- cbind(beadNum, 
                  rep(NA, nrow(beadNum)))
            }
        }
    }
    if (!is.null(detection)) {
        if (ncol(detection) == length(label)) {
            colnames(detection) <- label
        }
        else {
            warning("The number of detection columns does not match! Please check the input data format.\n")
            if (ncol(detection) > length(label)) 
                beadNum <- beadNum[, 1:length(label)]
            if (ncol(detection) < length(label)) {
                for (i in 1:(length(label) - ncol(detection))) detection <- cbind(detection, 
                  rep(NA, nrow(detection)))
            }
        }
    }
    if (is.null(se.exprs)) {
        cmd <- "x.lumi <- new(\"ExpressionSet\", exprs=exprs"
    }
    else {
        cmd <- "x.lumi <- new(\"LumiBatch\", exprs=exprs, se.exprs=se.exprs"
    }
    if (!is.null(detection)) 
        cmd <- paste(cmd, ", detection=detection")
    if (!is.null(beadNum)) 
        cmd <- paste(cmd, ", beadNum=beadNum")
    cmd <- paste(cmd, ", featureData=featureData")
    if (!all(sampleID == label)) {
        pData <- data.frame(sampleID = sampleID, label = label)
        rownames(pData) <- label
        varMetadata <- data.frame(labelDescription = c("The unique Illumina microarray Id", 
            "The label of the sample"))
        rownames(varMetadata) <- c("sampleID", "label")
    }
    else {
        pData <- data.frame(sampleID = sampleID)
        rownames(pData) <- sampleID
        varMetadata <- data.frame(labelDescription = c("The unique Illumina microarray Id"))
        rownames(varMetadata) <- c("sampleID")
    }
    pdata <- new("AnnotatedDataFrame", data = pData, varMetadata = varMetadata)
    cmd <- paste(cmd, ", phenoData=pdata")
    cmd <- paste(cmd, ")")
    eval(parse(text = cmd))
    if (is.null(se.exprs)) {
        if (checkDupId && convertNuID) {
            x.lumi <- addNuID2lumi(x.lumi, lib.mapping = lib.mapping)
        }
        options(stringsAsFactors = oldSetting)
        return(x.lumi)
    }
    controlData(x.lumi) <- controlData
    x.lumi@QC <- list(BeadStudioSummary = sampleSummary)
    sampleNames(x.lumi) <- label
    if (!is.null(info)) {
        info <- gsub("\t+", "\t", info)
    }
    notes(x.lumi) <- list(`Data File Information` = info)
    history.finished <- as.character(Sys.time())
    history.command <- capture.output(print(match.call(lumiR)))
    if (length(grep(",", history.command)) > 0) {
        history.command <- sub("\\(.+,", paste("(\"", fileName, 
            "\",", sep = ""), history.command)
    }
    else {
        history.command <- sub("\\(.+\\)", paste("(\"", fileName, 
            "\")", sep = ""), history.command)
    }
    lumiVersion <- packageDescription("lumi")$Version
    x.lumi@history <- data.frame(submitted = history.submitted, 
        finished = history.finished, command = history.command, 
        lumiVersion = lumiVersion)
    if (any(toupper(header) == "SPECIES")) {
        species <- as.character(allData[1, header[toupper(header) == 
            "SPECIES"]])
        annotation(x.lumi) <- switch(tolower(species), `homo sapiens` = "lumiHumanAll.db", 
            `mus musculus` = "lumiMouseAll.db", `rattus norvegicus` = "lumiRatAll.db", 
            human = "lumiHumanAll.db", mouse = "lumiMouseAll.db", 
            rat = "lumiRatAll.db", "NA")
    }
    if (QC) 
        x.lumi <- lumiQ(x.lumi, detectionTh = detectionTh, verbose = verbose)
    if (!convertNuID) 
        lib.mapping <- NULL
    if (convertNuID && !is.null(lib.mapping)) {
        x.lumi <- addNuID2lumi(x.lumi, lib.mapping = lib.mapping, 
            verbose = verbose)
    }
    options(stringsAsFactors = oldSetting)
    return(x.lumi)
}


lumiR.batch <- function (fileList, convertNuID = TRUE, lib.mapping = NULL, detectionTh = 0.01, 
    QC = TRUE, transform = c("none", "vst", "log2", "cubicRoot"), 
    sampleInfoFile = NULL, verbose = TRUE, ...) 
{
    oldDir <- getwd()
    dirMode <- FALSE
    transform <- match.arg(transform)
    if (file.exists(fileList[1])) {
        if (file.info(fileList[1])[1, "isdir"]) {
            dirMode <- TRUE
            setwd(fileList)
        }
    }
    if (dirMode && length(fileList) == 1) {
        fileList <- dir(fileList, pattern = ".csv")
        if (length(fileList) == 0) 
            stop("No data files were found!\n")
    }
    history.submitted <- as.character(Sys.time())
    if (verbose) {
        cat("Inputting the data ...\n")
        if (transform != "none") 
            cat(paste("Transformation", transform, "will be performed for each data file ...\n"))
    }
    for (i in 1:length(fileList)) {
        file.i <- fileList[i]
        if (transform != "none") {
            x.lumi.i <- lumiR(file.i, parseColumnName = FALSE, 
                convertNuID = FALSE, verbose = FALSE, ...)
            x.lumi.i <- lumiT(x.lumi.i, method = transform, simpleOutput = TRUE, 
                verbose = FALSE)
        }
        else {
            x.lumi.i <- lumiR(file.i, parseColumnName = FALSE, 
                convertNuID = FALSE, QC = FALSE, verbose = FALSE, 
                ...)
        }
        if (i == 1) {
            x.lumi <- x.lumi.i
        }
        else {
            x.lumi <- combine(x.lumi, x.lumi.i)
        }
    }
    if (!convertNuID) 
        lib.mapping <- NULL
    if (!is.null(lib.mapping) || convertNuID) {
        if (verbose) 
            cat("\nAdding nuID to the data ...\n")
        x.lumi <- addNuID2lumi(x.lumi, lib.mapping = lib.mapping)
    }
    if (!is.null(sampleInfoFile)) {
        if (is.character(sampleInfoFile) || class(sampleInfoFile)[1] == 
            "file") {
            if (file.exists(sampleInfoFile)) {
                sampleInfo <- read.table(sampleInfoFile, header = TRUE, 
                  sep = "\t", colClasses = "character", comment.char = "", 
                  check.names = FALSE)
            }
            else {
                warning("The provided sampleInfoFile does not exist\n!")
                setwd(oldDir)
                return(x.lumi)
            }
        }
        else if (is.data.frame(sampleInfoFile)) {
            sampleInfo <- sampleInfoFile
        }
        colName <- toupper(names(sampleInfo))
        ind <- grep("ID$", colName, ignore.case = TRUE)
        if (length(ind) == 0) {
            ID <- sampleInfo[, 1]
            if (any(duplicated(ID))) {
                warning("In sampleInfoFile, the ID column is required or the first column should be unique!\n")
                setwd(oldDir)
                return(x.lumi)
            }
            ind <- 1
        }
        else {
            ID <- sampleInfo[, ind[1]]
        }
        rownames(sampleInfo) <- ID
        colnames(sampleInfo)[ind[1]] <- "sampleID"
        sampleName <- sampleNames(x.lumi)
        ID <- ID[ID %in% sampleName]
        if (nrow(sampleInfo) != length(ID)) {
            warning("Some IDs provided in the sampleInfoFile do not exist the data file!\n")
            if (length(ID) == 0) {
                stop("The IDs provided in the sampleInfoFile do not match the data file!\n")
            }
        }
        x.lumi <- x.lumi[, ID]
        if (is.null(pData(phenoData(x.lumi)))) {
            pData <- sampleInfo[ID, ]
        }
        else {
            pData <- data.frame(pData(phenoData(x.lumi))[!(toupper(names(pData(phenoData(x.lumi)))) %in% 
                c(toupper(names(sampleInfo)), "ID", "SAMPLEID"))], 
                sampleInfo[ID, ])
        }
        label <- sampleInfo[ID, colName == "LABEL"]
        if (length(label) == length(ID)) {
            rownames(pData) <- label
            sampleNames(x.lumi) <- label
        }
        pdata <- new("AnnotatedDataFrame", data = pData)
        phenoData(x.lumi) <- pdata
    }
    if (is(x.lumi, "LumiBatch")) {
        history.finished <- as.character(Sys.time())
        history.command <- capture.output(print(match.call(lumiR.batch)))
        if (length(fileList) > 1) {
            fileList <- paste("c(\"", paste(fileList, collapse = "\",\"", 
                sep = ""), "\")", sep = "")
        }
        else {
            fileList <- paste("\"", fileList, "\"", sep = "")
        }
        if (length(grep(",", history.command)) > 0) {
            history.command <- sub("\\(.+,", paste("(", fileList, 
                ",", sep = ""), history.command)
        }
        else {
            history.command <- sub("\\(.+\\)", paste("(", fileList, 
                ")", sep = ""), history.command)
        }
        lumiVersion <- packageDescription("lumi")$Version
        x.lumi@history <- data.frame(submitted = history.submitted, 
            finished = history.finished, command = history.command, 
            lumiVersion = lumiVersion)
    }
    if (QC) 
        x.lumi <- lumiQ(x.lumi, detectionTh = detectionTh, verbose = verbose)
    setwd(oldDir)
    return(x.lumi)
}
