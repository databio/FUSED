# PACKAGE DOCUMENTATION
#' FUSED R package
#'
#'
#' FUSED is a collection of functions to determine the presence of hallmarks
#' of common DNA repair mechanisms including NHEJ, MMEJ, SSA, BIR, and MMBIR.
#'
#' @docType package
#' @name fused
#' @author Jason Smith, Ben Morris
#'
#' @references \url{http://github.com/databio/fused}
NULL

################################################################################
# FUNCTIONS

#' Checks if provided genome is an available BSgenome object,
#' and if so, returns it. If the library is not installed, it issues a warning
#' and returns NULL.
#'
#' @param bsg A BSgenome compatible genome string or BSgenome object.
#' @return A BSgenome object if installed.
.validBSgenome = function(bsg) {
    if(inherits(bsg, "BSgenome")){
        return(bsg)
    }
    if (requireNamespace(bsg))
        return(utils::getAnywhere(bsg)$objs[[1]])
    else
        warning(bsg, " is not installed")
    return(NULL)
}

#' Splits a data.table object into a list of data.tables based on an
#' individual column. Returns a list of data.tables if column to split
#' is present, otherwise returns NULL.
#'
#' @param DT A data.table object.
#' @param split_factor A column name in DT to use as a factor to split and
#'                     name list elements. Can be a character vector
#'                     or an integer.
#' @return A list of data.tables.
splitDataTable <- function(dt, split_factor) {
    stopifnot(any(class(dt) == c("data.table")))
    if (split_factor %in% colnames(DT)) {
        factor_order = unique(DT[, get(split_factor)])
        if (is.numeric(split_factor)) {
            split_factor = colnames(DT)[split_factor]
            message("Integer split_factor, changed to: ", split_factor)
        }
        l = lapply(split(seq_len(nrow(DT)), DT[, get(split_factor)]),
                   function(x) DT[x])
        return(l[factor_order])
    } else {
        warning(split_factor, " was not found as a column in the data.table.")
        return(NULL)
    }
}

#' Checks if provided object is an ORIEN derived fusion prediction data.table
#'
#' @param x A data.table object containing ORIEN predicted fusions.
#' @return A warning if the object does not appear to be from ORIEN
.validORIEN <- function(x) {
    stopifnot(!is.null(x))
    stopifnot(any(class(x) == "data.frame" | class(x) == "data.table"))
    ORIENnames <- c("Patient_ID",
                    "Sequencing_ID",
                    "#FusionName",
                    "strand2(gene/fusion)",
                    "LeftBreakpoint",
                    "RightBreakpoint",
                    "JunctionReadCount",
                    "SpanningFragCount",
                    "FFPM",
                    "direction1",
                    "direction2",
                    "PROT")
    minORIENnames <- c("#FusionName",
                       "Sequencing_ID",
                       "LeftBreakpoint",
                       "RightBreakpoint",
                       "direction1",
                       "direction2")
    if(!all(ORIENnames %in% colnames(x))) {
        warning("Cannot validate input as an ORIEN fusions.tsv file.")
    }
    if (!all(minORIENnames %in% colnames(x))) {
        stop("Cannot validate minimal required columns for ORIEN fusion predictions.")
    }
}

#' Checks if provided object is a STAR derived fusion prediction data.table
#'
#' @param x A data.table object containing STAR predicted fusions.
#' @return A warning if the object does not appear to be from STAR.
.validSTAR <- function(x) {
    stopifnot(!is.null(x))
    stopifnot(any(class(x) == "data.frame" | class(x) == "data.table"))
    STARnames <- c("#FusionName",
                   "JunctionReadCount",
                   "SpanningFragCount",
                   "SpliceType",
                   "LeftGene",
                   "LeftBreakpoint",
                   "RightGene",
                   "RightBreakpoint",
                   "LargeAnchorSupport",
                   "FFPM",
                   "LeftBreakDinuc",
                   "LeftBreakEntropy",
                   "RightBreakDinuc",
                   "RightBreakEntropy",
                   "annots")
    minSTARnames <- c("#FusionName",
                      "LeftBreakpoint",
                      "RightBreakpoint")
    if(!all(STARnames %in% colnames(x))) {
        warning("Cannot validate input as a STAR star-fusion.fusion_predictions.tsv file.")
    }
    if (!all(minSTARnames %in% colnames(x))) {
        stop("Cannot validate minimal required columns for STAR fusion predictions.")
    } else {
        # Validate breakpoint values
        if (length(tstrsplit(x$LeftBreakpoint, ":", fixed=TRUE)) != 3) {
            stop("Cannot validate LeftBreakpoint column as a STAR fusion prediction.")
        }
        if (length(tstrsplit(x$RightBreakpoint, ":", fixed=TRUE)) != 3) {
            stop("Cannot validate RightBreakpoint column as a STAR fusion prediction.")
        }
    }
}

#' Checks if provided object is an ARRIBA derived fusion prediction data.table
#'
#' @param x A data.table object containing ARRIBA predicted fusions.
#' @return A warning if the object does not appear to be from ARRIBA
.validARRIBA <- function(x) {
    stopifnot(!is.null(x))
    stopifnot(any(class(x) == "data.frame" | class(x) == "data.table"))
    ARRIBAnames <- c("#gene1",
                     "gene2",
                     "strand1(gene/fusion)",
                     "strand2(gene/fusion)",
                     "breakpoint1",
                     "breakpoint2",
                     "site1",
                     "site2",
                     "type",
                     "split_reads1",
                     "split_reads2",
                     "discordant_mates",
                     "coverage1",
                     "coverage2",
                     "confidence",
                     "reading_frame",
                     "tags",
                     "retained_protein_domains",
                     "closest_genomic_breakpoint1",
                     "closest_genomic_breakpoint2",
                     "gene_id1",
                     "gene_id2",
                     "transcript_id1",
                     "transcript_id2",
                     "direction1",
                     "direction2",
                     "filters",
                     "fusion_transcript",
                     "peptide_sequence",
                     "read_identifiers")
    minARRIBAnames <- c("#gene1",
                        "gene2",
                        "breakpoint1",
                        "breakpoint2",
                        "direction1",
                        "direction2")
    if(!all(ARRIBAnames %in% colnames(x))) {
        warning("Cannot validate input as an ARRIBA fusions.tsv file.")
    }
    if (!all(minARRIBAnames %in% colnames(x))) {
        stop("Cannot validate minimal required columns for ARRIBA fusion predictions.")
    }
}

#' Modify a valid ORIEN fusion data.table for downstream analysis. ORIEN uses
#' a merged output derived from both STAR and ARRIBA fusion predictions.
#'
#' @param x A valid data.table object containing ORIEN predicted fusions.
#' @param window A numeric value by which to expand the region around the
#'               breakpoint
#' @return A fusion data.table object
constructORIEN <- function(x, window = 99) {
    # Create custom id by combining Sequencing_ID and #FusionName
    x[,id:=paste0(x$Sequencing_ID, "_", x$`#FusionName`)]
    # Split out chromosomes and coordinates
    x[, c("lChr", "lStart") := tstrsplit(LeftBreakpoint, ":", fixed=TRUE)]
    x[, c("rChr", "rStart") := tstrsplit(RightBreakpoint, ":", fixed=TRUE)]
    # Split out gene names
    x[, c("lGene", "rGene") := tstrsplit(`#FusionName`, "--", fixed=TRUE)]
    # Determine end coordinates
    x$lStart <- as.numeric(x$lStart)
    x$rStart <- as.numeric(x$rStart)
    x[,lEnd:=ifelse(x$direction1 == "upstream",
                    (lStart - window),
                    (lStart + window))]
    x[,rEnd:=ifelse(x$direction2 == "upstream",
                    (rStart - window),
                    (rStart + window))]
    # Must swap start and ends for upstream scenario
    x[,lEnd:=ifelse(x$direction1 == "upstream",
                    (lStart - window),
                    (lStart + window))]
    x[,rEnd:=ifelse(x$direction2 == "upstream",
                    (rStart - window),
                    (rStart + window))]
    # Make fusion name unique
    x$uniqueID <- make.unique(x$id)
    return(x)
}

#' Modify a valid STAR fusion data.table for downstream analysis
#'
#' @param x A valid data.table object containing STAR predicted fusions.
#' @param window A numeric value by which to expand the region around the
#'               breakpoint
#' @return A fusion data.table object
constructSTAR <- function(x, window = 99) {
    # Create fusion name
    x[,fusion_name:=x$`#FusionName`]
    # Split out chromosomes and coordinates
    x[, c("lChr", "lStart", "direction1") := tstrsplit(LeftBreakpoint, ":", fixed=TRUE)]
    x[, c("rChr", "rStart", "direction2") := tstrsplit(RightBreakpoint, ":", fixed=TRUE)]
    # Split out gene names
    x[, c("lGene", "rGene") := tstrsplit(`#FusionName`, "--", fixed=TRUE)]
    # Determine end coordinates
    x$lStart <- as.numeric(x$lStart)
    x$rStart <- as.numeric(x$rStart)
    x[,lEnd:=ifelse(x$direction1 == "upstream" | x$direction1 == "+",
                    (lStart - window),
                    (lStart + window))]
    x[,rEnd:=ifelse(x$direction2 == "upstream" | x$direction2 == "+",
                    (rStart - window),
                    (rStart + window))]
    # Must swap start and ends for upstream scenario
    x[(direction1 == "upstream" | direction1 == "+"), c("lStart", "lEnd") := .(lEnd, lStart)]
    x[(direction2 == "upstream" | direction2 == "+"), c("rStart", "rEnd") := .(rEnd, rStart)]
    # Make fusion name unique
    x$uniqueID <- make.unique(x$fusion_name)
    return(x)
}

#' Modify a valid ARRIBA fusion data.table for downstream analysis
#'
#' @param x A valid data.table object containing ARRIBA predicted fusions.
#' @param window A numeric value by which to expand the region around the
#'               breakpoint
#' @return A fusion data.table object
constructARRIBA <- function(x, window = 99) {
    # Create fusion name
    x[,fusion_name:=paste0(x$`#gene1`, "--", x$`gene2`)]
    # Assign gene names
    x[,lGene:=x$`#gene1`]
    x[,rGene:=x$`gene2`]
    # Split out chromosomes and coordinates
    x[, c("lChr", "lStart") := tstrsplit(breakpoint1, ":", fixed=TRUE)]
    x[, c("rChr", "rStart") := tstrsplit(breakpoint2, ":", fixed=TRUE)]
    # Determine end coordinates
    x$lStart <- as.numeric(x$lStart)
    x$rStart <- as.numeric(x$rStart)
    x[,lEnd:=ifelse(x$direction1 == "upstream",
                    (lStart - window),
                    (lStart + window))]
    x[,rEnd:=ifelse(x$direction2 == "upstream",
                    (rStart - window),
                    (rStart + window))]
    # Must swap start and ends for upstream scenario
    # Must swap start and ends for upstream scenario
    x[(direction1 == "upstream" | direction1 == "+"), c("lStart", "lEnd") := .(lEnd, lStart)]
    x[(direction2 == "upstream" | direction2 == "+"), c("rStart", "rEnd") := .(rEnd, rStart)]
    # Make fusion name unique
    x$uniqueID <- make.unique(x$fusion_name)
    return(x)
}

#' Helper function to calculate the optimal pairwise alignment of a fusion
#' of type, 'mechanism'.
#'
#' Returns a data.table for the fusion with the best scoring alignment
#' (or no alignment) for each repair pathway.
#'
#' @param DT A data.table object containing fusion predictions.
#' @param bsg A BSgenome string or BSgenome object.
#' @param mechanism A character vector representing which repair pathway to
#'                  investigate in the given fusion data.table.
#' @param window A numeric value representing how upstream and downstream of
#'               a fusion breakpoint within which to perform sequence alignment.
#' @param gapOpening A numeric value representing the cost for opening a gap
#'                   in the alignment.
#' @param gapExtension A numeric value representing the incremental cost for
#'                     extending the gap in alignment.
#' @return A data.table object with identified homologies.
calcPairwiseAlignment <- function(DT,
                                  bsg,
                                  mechanism = c("NHEJ", "SSA", "MMEJ", "BIR", "MMBIR"),
                                  window=4,
                                  gapOpening=Inf,
                                  gapExtension=Inf,
                                  ...) {
    # NHEJ: +/- 4 nt (can only be in 4 nt window around or either side of break)
    # SSA: +/- 70 nt
    # MMEJ: +/- 20 nt (+/- 50 albreakpointsable)
    # BIR: +/- 100 nt
    # requires UCSC compatible chr name
    geneA <- data.table(gene_name = DT$lGene,
                        chr = ifelse(stringr::str_detect(DT$lChr, "chr"),
                                     DT$lChr,
                                     paste0("chr", DT$lChr)),
                        start = DT$lEnd - window + 1,
                        end = DT$lEnd,
                        fusion_name = DT$uniqueID)
    # requires UCSC compatible chr name
    geneB <- data.table(gene_name = DT$rGene,
                        chr = ifelse(stringr::str_detect(DT$rChr, "chr"),
                                     DT$rChr,
                                     paste0("chr", DT$rChr)),
                        start = DT$rStart,
                        end = DT$rStart + window - 1,
                        fusion_name = DT$uniqueID)
    geneGR <- GenomicRanges::makeGRangesFromDataFrame(rbind(geneA, geneB),
                                                      keep.extra.columns=T)
    bsg       <- .validBSgenome(bsg)
    geneSeq   <- Biostrings::getSeq(bsg, geneGR)
    seqAlign  <- Biostrings::pairwiseAlignment(geneSeq[[1]], geneSeq[[2]],
                                               type = "local",
                                               gapOpening=gapOpening,
                                               gapExtension=gapExtension)
    hDT <- data.table(type = mechanism,
                      uniqueID = DT$uniqueID,
                      fusionName = DT$fusion_name,
                      lGene = DT$lGene,
                      lChr = DT$lChr,
                      lStart = DT$lStart,
                      lEnd = DT$lEnd,
                      lSeq = as.character(geneSeq[[1]]),
                      rGene = DT$rGene,
                      rChr = DT$rChr,
                      rStart = DT$rStart,
                      rEnd = DT$rEnd,
                      rSeq = as.character(geneSeq[[2]]),
                      pStart = Biostrings::start(Biostrings::pattern(seqAlign)),
                      pEnd = Biostrings::end(Biostrings::pattern(seqAlign)),
                      sStart = Biostrings::start(Biostrings::subject(seqAlign)),
                      sEnd = Biostrings::end(Biostrings::subject(seqAlign)),
                      hLength = Biostrings::nchar(seqAlign), # length of homology
                      hScore = Biostrings::score(seqAlign),
                      gapCount = stringr::str_count(Biostrings::pattern(seqAlign), "-")
    )
    return(hDT)
}

#' Checks a fusion for hallmarks of microhomology error prone repairs.
#' Returns a data.table for the fusion with the best scoring alignment
#' (or no alignment) for each repair pathway.
#'
#' @param DT A data.table object representing a single fusion.
#' @param bsg A BSgenome string or BSgenome object.
#' @param mechanism A character vector representing which repair pathway to
#'                  investigate in the given fusion data.table.
#' @param window A numeric value representing how upstream and downstream of
#'               a fusion breakpoint within which to perform sequence alignment.
#' @param gapOpening A numeric value representing the cost for opening a gap
#'                   in the alignment.
#' @param gapExtension A numeric value representing the incremental cost for
#'                     extending the gap in alignment.
#' @param source A character vector representing the fusion prediction source.
#' @return A data.table object with fusion alignments for specified mechanism.
calcHomology <- function(DT,
                         bsg,
                         mechanism = c("NHEJ", "SSA", "MMEJ", "BIR", "MMBIR"),
                         window = NA,
                         gapOpening = NA,
                         gapExtension = NA,
                         source = c("STAR", "ARRIBA", "ORIEN"),
                         ...) {
    options(warn=0)
    DT <- as.data.table(DT)
    if (nrow(DT) > 1) {
        rbindlist(lapply(split(DT,seq(nrow(DT))), calcHomology,
                         bsg=bsg, mechanism=mechanism, window=window,
                         gapOpening=gapOpening, gapExtension=gapExtension,
                         source=source))
    } else {
        if (source == "STAR") {
            .validSTAR(DT)
            DT <- constructSTAR(DT)
        } else if (source == "ARRIBA") {
            .validARRIBA(DT)
            DT <- constructARRIBA(DT)
        } else if (source == "ORIEN") {
            .validORIEN(DT)
            DT <- constructORIEN(DT)
        } else {
            stop("Could not identify source of fusion predictions.")
        }
        if (toupper(mechanism) == "NHEJ") {
            # Default window is (+/-) 4 nt upstream and downstream of breakpoint
            if (is.na(window)) {
                window <- 4
            }
            if (is.na(gapOpening)) {
                gapOpening <- Inf
            }
            if (is.na(gapExtension)) {
                gapExtension <- 4
            }
            homo <- calcPairwiseAlignment(DT,
                                          bsg=bsg,
                                          mechanism = "NHEJ",
                                          window = window,
                                          gapOpening = gapOpening,
                                          gapExtension = gapExtension,
                                          ...)
            # Score must be positive, and the template must include
            # the breakpoint
            if (homo$hScore > 0 & (homo$hLength <= 4 &
                                   (homo$pEnd == 4 | homo$sStart == 1))) {
                homo[, criteria := "pass"]
            } else {
                homo[, criteria := "fail"]
            }
            homo[] # Ensure returned data.table prints on first return
            summary(warnings())
            options(warn=1)
            return(homo)
        } else if (toupper(mechanism) == "SSA") {
            # Default window is (+/-) 70 nt upstream and downstream
            # of breakpoint
            if (is.na(window)) {
                window <- 70
            }
            if (is.na(gapOpening)) {
                gapOpening <- Inf
            }
            if (is.na(gapExtension)) {
                gapExtension <- 4
            }
            homo <- calcPairwiseAlignment(DT,
                                          bsg=bsg,
                                          mechanism = "SSA",
                                          window = window,
                                          gapOpening = gapOpening,
                                          gapExtension = gapExtension,
                                          ...)
            # Score must be positive, and the template must be between
            # 30 and 70 nt in length.
            if (homo$hScore > 0 & (homo$hLength >= 30 & homo$hLength <= 70)) {
                homo[, criteria := "pass"]
            } else {
                homo[, criteria := "fail"]
            }
            homo[] # Ensure returned data.table prints on first return
            summary(warnings())
            options(warn=1)
            return(homo)
        } else if (toupper(mechanism) == "MMEJ") {
            # Default window is (+/-) 20 nt upstream and downstream
            # of breakpoint
            if (is.na(window)) {
                window <- 20
            }
            if (is.na(gapOpening)) {
                gapOpening <- Inf
            }
            if (is.na(gapExtension)) {
                gapExtension <- 4
            }
            homo <- calcPairwiseAlignment(DT,
                                          bsg=bsg,
                                          mechanism = "MMEJ",
                                          window = window,
                                          gapOpening = gapOpening,
                                          gapExtension = gapExtension,
                                          ...)
            # Score must be positive, and the template must be less than 20 nt
            if (homo$hScore > 0 & homo$hLength <= 20) {
                homo[, criteria := "pass"]
            } else {
                homo[, criteria := "fail"]
            }
            homo[] # Ensure returned data.table prints on first return
            summary(warnings())
            options(warn=1)
            return(homo)
        } else if (toupper(mechanism) == "BIR") {
            # Default window is (+/-) 100 nt upstream and downstream
            # of breakpoint
            if (is.na(window)) {
                window <- 100
            }
            if (is.na(gapOpening)) {
                gapOpening <- Inf
            }
            if (is.na(gapExtension)) {
                gapExtension <- 4
            }
            homo <- calcPairwiseAlignment(DT,
                                          bsg=bsg,
                                          mechanism = "BIR",
                                          window = window,
                                          gapOpening = gapOpening,
                                          gapExtension = gapExtension,
                                          ...)
            # Score must be positive, and the template must be between
            # 70 and 100 nt in length.
            if (homo$hScore > 0 & (homo$hLength >= 70 & homo$hLength <= 100)) {
                homo[, criteria := "pass"]
            } else {
                homo[, criteria := "fail"]
            }
            homo[] # Ensure returned data.table prints on first return
            summary(warnings())
            options(warn=1)
            return(homo)
        } else if (toupper(mechanism) == "MMBIR") {
            # Default window is (+/-) 20 nt upstream and downstream
            # of breakpoint
            if (is.na(window)) {
                window <- 20
            }
            if (is.na(gapOpening)) {
                gapOpening <- 0
            }
            if (is.na(gapExtension)) {
                gapExtension <- 0
            }
            # MMBIR allows for gaps. Default is no cost for opening or
            # extending gap.
            homo <- calcPairwiseAlignment(DT,
                                          bsg=bsg,
                                          mechanism = "MMBIR",
                                          window = window,
                                          gapOpening = gapOpening,
                                          gapExtension = gapExtension,
                                          ...)
            # Score must be positive, and the template must include
            # the breakpoint and be less than 20 nt in length.
            if (homo$hScore > 0 & homo$hLength <= 20 & homo$gapCount > 0 &
                (homo$gapCount <= (0.2*homo$hLength)) &
                (homo$pEnd == 10 | homo$sStart == 1)) {
                homo[, criteria := "pass"]
            } else {
                homo[, criteria := "fail"]
            }
            homo[] # Ensure returned data.table prints on first return
            summary(warnings())
            options(warn=1)
            return(homo)
        } else {
            warning(mechanism, "is not recognized as a valid mechanism.")
            summary(warnings())
            options(warn=1)
            return(NULL)
        }
    }
}

