#' Read an ONT modkit bedMethyl file
#'
#' Create an \code{RBedMethyl} object backed by HDF5Array from a nanoporetech
#' modkit bedMethyl file (headerless).
#'
#' @param bedmethyl Path to a nanoporetech modkit bedMethyl file (optionally gzipped).
#' @param mod Modification code to retain (\code{"m"} or \code{"h"}).
#' @param chunk_size Reserved for future use.
#' @param h5file Path to the HDF5 file to create. Defaults to a deterministic
#'   path in \code{tempdir()} derived from the input \code{bedmethyl} filename,
#'   so subsequent calls reuse the same file.
#' @param check_sorted Logical, check that records are sorted by chrom and chromStart.
#' @param fields Character vector of numeric fields to load. Defaults to \code{c("coverage", "mod_reads")}.
#'
#' @return An \code{RBedMethyl} object.
#' @export
#' @importFrom stats setNames
#' @importFrom methods new
#' @examples
#' lines <- c(
#'   paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
#'   paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
#' )
#' tmp <- tempfile(fileext = ".bed")
#' writeLines(lines, tmp)
#' bm <- readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))
#' bm
readBedMethyl <- function(bedmethyl, mod = "m", chunk_size = 5e6,
                          h5file = NULL, check_sorted = TRUE,
                          fields = c("coverage", "mod_reads")) {
  if (!mod %in% c("m", "h", "a")) {
    stop("mod must be one of: 'm', 'h', 'a'.")
  }
  if (is.null(h5file)) {
    base <- basename(bedmethyl)
    base <- sub("\\.gz$", "", base, ignore.case = TRUE)
    base <- sub("\\.[^.]+$", "", base)

    hash <- tryCatch(as.character(tools::md5sum(bedmethyl)),
      error = function(e) NA_character_
    )
    if (is.na(hash)) {
      key <- tryCatch(normalizePath(bedmethyl, winslash = "/", mustWork = FALSE),
        error = function(e) bedmethyl
      )
      hash <- paste0("path_", gsub("[^A-Za-z0-9]", "_", key))
    }

    h5file <- file.path(tempdir(), paste0(base, "_", hash, ".h5"))
  }
  required_cols <- c(
    "chrom", "chromStart", "chromEnd", "mod", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "coverage", "pct",
    "mod_reads", "unmod_reads", "other_reads",
    "del_reads", "fail_reads", "diff_reads", "nocall_reads"
  )

  stream_cmd <- function(path) {
    is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
    base_cmd <- if (is_gz) "zcat" else "cat"
    sprintf("%s %s", base_cmd, shQuote(path))
  }

  check_sorted_chunk <- function(dt, chrom_levels, last_chr, last_start) {
    if (nrow(dt) == 0L) {
      return(list(last_chr = last_chr, last_start = last_start))
    }
    codes <- match(dt$chrom, chrom_levels)
    if (any(is.na(codes))) stop("Unknown chromosome encountered during sorting check.")

    if (!is.null(last_chr)) {
      if (codes[1] < last_chr || (codes[1] == last_chr && dt$chromStart[1] < last_start)) {
        stop("bedMethyl records are not sorted by chrom and chromStart.")
      }
    }

    if (nrow(dt) > 1L) {
      prev_codes <- codes[-length(codes)]
      next_codes <- codes[-1L]
      prev_start <- dt$chromStart[-nrow(dt)]
      next_start <- dt$chromStart[-1L]
      bad <- next_codes < prev_codes | (next_codes == prev_codes & next_start < prev_start)
      if (any(bad)) stop("bedMethyl records are not sorted by chrom and chromStart.")
    }

    list(last_chr = codes[nrow(dt)], last_start = dt$chromStart[nrow(dt)])
  }

  numeric_fields <- c(
    "coverage", "pct", "mod_reads", "unmod_reads", "other_reads",
    "del_reads", "fail_reads", "diff_reads", "nocall_reads"
  )
  invalid_compat <- c("score", "thickStart", "thickEnd", "itemRgb")
  if (any(fields %in% invalid_compat)) {
    stop("Unsupported fields requested (legacy compatibility columns): ",
         paste(intersect(fields, invalid_compat), collapse = ", "))
  }
  invalid_fields <- setdiff(fields, numeric_fields)
  if (length(invalid_fields) > 0L) {
    stop("Unsupported fields requested: ", paste(invalid_fields, collapse = ", "))
  }

  fields <- unique(c(fields, "coverage"))
  strand_levels <- c("+", "-")

  ds_name <- function(name) paste0(mod, "_", name)

  # Check if HDF5 file already exists with valid datasets and metadata - reuse if so
  if (file.exists(h5file)) {
    existing_ok <- tryCatch({
      h5_contents <- rhdf5::h5ls(h5file)$name
      required_ds <- c(ds_name("chrom"), ds_name("chromStart"), ds_name("chromEnd"), ds_name("strand"))
      field_ds <- vapply(fields, ds_name, character(1))
      datasets_ok <- all(c(required_ds, field_ds) %in% h5_contents)

      # Also verify metadata attributes exist
      if (datasets_ok) {
        attrs <- rhdf5::h5readAttributes(h5file, ds_name("chrom"))
        attrs_ok <- !is.null(attrs$chrom_levels) && !is.null(attrs$chr_index)
        datasets_ok && attrs_ok
      } else {
        FALSE
      }
    }, error = function(e) FALSE,
    finally = rhdf5::h5closeAll())

    if (existing_ok) {
      # Load metadata from HDF5 attributes
      attrs <- rhdf5::h5readAttributes(h5file, ds_name("chrom"))
      chrom_levels <- as.character(attrs$chrom_levels)
      chr_index <- attrs$chr_index
      if (!is.matrix(chr_index)) chr_index <- as.matrix(chr_index)
      if (!is.null(attrs$chr_index_rownames)) {
        rownames(chr_index) <- as.character(attrs$chr_index_rownames)
      }
      if (!is.null(attrs$chr_index_colnames)) {
        colnames(chr_index) <- as.character(attrs$chr_index_colnames)
      }
      rhdf5::h5closeAll()

      assays <- S4Vectors::SimpleList(
        chrom = HDF5Array::HDF5Array(h5file, ds_name("chrom")),
        chromStart = HDF5Array::HDF5Array(h5file, ds_name("chromStart")),
        chromEnd = HDF5Array::HDF5Array(h5file, ds_name("chromEnd")),
        strand = HDF5Array::HDF5Array(h5file, ds_name("strand"))
      )
      for (field in fields) {
        assays[[field]] <- HDF5Array::HDF5Array(h5file, ds_name(field))
      }

      return(new("RBedMethyl",
        assays = assays,
        chrom_levels = chrom_levels,
        strand_levels = strand_levels,
        chr_index = chr_index,
        index = seq_along(assays$coverage),
        mod = mod
      ))
    } else {
      # File exists but is invalid/incomplete - close handles and remove it
      rhdf5::h5closeAll()
      unlink(h5file)
    }
  }

  preview <- data.table::fread(
    cmd = stream_cmd(bedmethyl),
    header = FALSE,
    sep = "\t",
    quote = "",
    nrows = 1
  )
  if (ncol(preview) < length(required_cols)) {
    stop("Unexpected number of columns in bedMethyl file.")
  }

  core_fields <- c("chrom", "chromStart", "chromEnd", "strand", "mod")
  select_fields <- unique(c(core_fields, fields))
  col_map <- setNames(seq_along(required_cols), required_cols)
  select_idx <- unname(col_map[select_fields])

  dt <- data.table::fread(
    cmd = stream_cmd(bedmethyl),
    header = FALSE,
    sep = "\t",
    quote = "",
    select = select_idx
  )
  data.table::setnames(dt, select_fields)

  int_cols <- intersect(c(
    "chromStart", "chromEnd", "score", "thickStart", "thickEnd",
    "coverage", "mod_reads", "unmod_reads", "other_reads",
    "del_reads", "fail_reads", "diff_reads", "nocall_reads"
  ), names(dt))
  num_cols <- intersect(c("pct"), names(dt))

  for (nm in int_cols) dt[[nm]] <- as.integer(dt[[nm]])
  for (nm in num_cols) dt[[nm]] <- as.numeric(dt[[nm]])

  if (check_sorted) {
    chrom_levels_sort <- unique(dt$chrom)
    check_sorted_chunk(dt, chrom_levels_sort, NULL, NULL)
  }

  dt <- dt[dt[["mod"]] == mod, , drop = FALSE]
  if (nrow(dt) == 0) stop("No records found for requested modification.")

  chrom_levels <- unique(dt$chrom)

  rle_chr <- rle(dt$chrom)
  ends <- cumsum(rle_chr$lengths)
  starts <- ends - rle_chr$lengths + 1L
  chr_index <- cbind(starts, ends)
  rownames(chr_index) <- rle_chr$values

  HDF5Array::writeHDF5Array(as.matrix(match(dt$chrom, chrom_levels)),
    filepath = h5file,
    name = ds_name("chrom")
  )
  HDF5Array::writeHDF5Array(as.matrix(dt$chromStart),
    filepath = h5file,
    name = ds_name("chromStart")
  )
  HDF5Array::writeHDF5Array(as.matrix(dt$chromEnd),
    filepath = h5file,
    name = ds_name("chromEnd")
  )
  HDF5Array::writeHDF5Array(as.matrix(match(dt$strand, strand_levels)),
    filepath = h5file,
    name = ds_name("strand")
  )

  # Store metadata as HDF5 attributes for reuse
  fid <- rhdf5::H5Fopen(h5file)
  did <- rhdf5::H5Dopen(fid, ds_name("chrom"))
  rhdf5::h5writeAttribute(chrom_levels, did, "chrom_levels")
  rhdf5::h5writeAttribute(chr_index, did, "chr_index")
  rhdf5::h5writeAttribute(rownames(chr_index), did, "chr_index_rownames")
  if (!is.null(colnames(chr_index))) {
    rhdf5::h5writeAttribute(colnames(chr_index), did, "chr_index_colnames")
  }
  rhdf5::H5Dclose(did)
  rhdf5::H5Fclose(fid)

  assays <- S4Vectors::SimpleList(
    chrom = HDF5Array::HDF5Array(h5file, ds_name("chrom")),
    chromStart = HDF5Array::HDF5Array(h5file, ds_name("chromStart")),
    chromEnd = HDF5Array::HDF5Array(h5file, ds_name("chromEnd")),
    strand = HDF5Array::HDF5Array(h5file, ds_name("strand"))
  )

  for (field in fields) {
    HDF5Array::writeHDF5Array(as.matrix(dt[[field]]),
      filepath = h5file,
      name = ds_name(field)
    )
    assays[[field]] <- HDF5Array::HDF5Array(h5file, ds_name(field))
  }

  new("RBedMethyl",
    assays = assays,
    chrom_levels = chrom_levels,
    strand_levels = strand_levels,
    chr_index = chr_index,
    index = seq_along(assays$coverage),
    mod = mod
  )
}

#' List retrievable bedMethyl fields
#'
#' Returns a data.frame describing retrievable bedMethyl fields and their types.
#'
#' @return A \code{data.frame} with columns \code{field}, \code{type}, and \code{description}.
#' @export
#' @examples
#' bedMethylFields()
bedMethylFields <- function() {
  data.frame(
    field = c(
      "coverage", "pct", "mod_reads", "unmod_reads", "other_reads",
      "del_reads", "fail_reads", "diff_reads", "nocall_reads"
    ),
    type = c(
      "int", "float", "int", "int", "int",
      "int", "int", "int", "int"
    ),
    description = c(
      "Valid coverage (Nvalid_cov).",
      "Fraction modified (Nmod / Nvalid_cov).",
      "Number of modified calls (Nmod).",
      "Number of canonical calls (Ncanonical).",
      "Number of other-mod calls (Nother_mod).",
      "Number of deletions at reference position (Ndelete).",
      "Number of failed calls below threshold (Nfail).",
      "Number of non-canonical base calls (Ndiff).",
      "Number of no-call reads with canonical base (Nnocall)."
    ),
    stringsAsFactors = FALSE
  )
}
