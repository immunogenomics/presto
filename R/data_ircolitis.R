#' Download and load the ircolitis tissue CD8 demo dataset
#'
#' Downloads `GSE206299_ircolitis-tissue-cd8.h5ad.gz` from GEO on first call
#' and caches it under [tools::R_user_dir()]. Subsequent calls reuse the cache.
#' Used to demonstrate presto on a realistic single-cell dataset.
#'
#' Source: Thomas et al. (Nat. Med. 2024), GSE206299. 25,341 colon-tissue
#' CD8 T cells from immune checkpoint colitis patients and controls.
#'
#' Requires the suggested package `rhdf5` (Bioconductor).
#'
#' @param cache_dir Directory for the cached `.h5ad` file. Defaults to
#'   `tools::R_user_dir("presto", which = "cache")`.
#' @param overwrite If `TRUE`, re-download even if the file is cached.
#' @param verbose Print download / parsing progress.
#'
#' @return A list with elements:
#' \itemize{
#'   \item `counts` - gene-by-cell `dgCMatrix` of raw UMI counts.
#'   \item `obs` - cell metadata `data.frame` (categorical fields decoded
#'     to character; `UMAP1`/`UMAP2` and `PC1..PC24` kept as numeric).
#'   \item `var` - gene metadata `data.frame`.
#' }
#'
#' @examples
#' \dontrun{
#'   d <- load_ircolitis_cd8()
#'   dim(d$counts)
#'   table(d$obs$cluster)
#'   res <- wilcoxauc(d$counts, d$obs$cluster)
#' }
#'
#' @export
load_ircolitis_cd8 <- function(
    cache_dir = tools::R_user_dir("presto", which = "cache"),
    overwrite = FALSE,
    verbose = TRUE
) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop(
            "Package 'rhdf5' is required. Install with:\n",
            "  BiocManager::install('rhdf5')"
        )
    }

    url <- paste0(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206299/",
        "suppl/GSE206299_ircolitis-tissue-cd8.h5ad.gz"
    )
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    gz_path <- file.path(cache_dir, "ircolitis-tissue-cd8.h5ad.gz")
    h5_path <- file.path(cache_dir, "ircolitis-tissue-cd8.h5ad")

    if (overwrite || !file.exists(h5_path)) {
        if (overwrite || !file.exists(gz_path)) {
            if (verbose) message("Downloading: ", url)
            utils::download.file(url, gz_path, mode = "wb")
        }
        if (verbose) message("Decompressing: ", gz_path)
        R.utils::gunzip(gz_path, destname = h5_path, remove = FALSE,
                        overwrite = TRUE)
    } else if (verbose) {
        message("Using cached: ", h5_path)
    }

    parse_ircolitis_h5ad(h5_path, verbose = verbose)
}

#' @noRd
parse_ircolitis_h5ad <- function(path, verbose = TRUE) {
    if (verbose) message("Reading X (sparse counts)...")
    shape <- rhdf5::h5readAttributes(path, "/X")$shape
    n_cells <- shape[1]
    n_genes <- shape[2]
    x_data    <- rhdf5::h5read(path, "/X/data")
    x_indices <- rhdf5::h5read(path, "/X/indices")
    x_indptr  <- rhdf5::h5read(path, "/X/indptr")

    ## AnnData CSC: rows = cells, cols = genes. Build cell-by-gene dgCMatrix,
    ## then transpose to the gene-by-cell layout that presto expects.
    cell_by_gene <- methods::new(
        "dgCMatrix",
        i = as.integer(x_indices),
        p = as.integer(x_indptr),
        x = as.numeric(x_data),
        Dim = c(as.integer(n_cells), as.integer(n_genes))
    )
    counts <- Matrix::t(cell_by_gene)

    if (verbose) message("Reading var (gene metadata)...")
    var_df <- read_h5_table(path, "/var")
    ## /var/_index is "ENSG...|SYMBOL", already globally unique. Keep it
    ## verbatim as rownames; expose ensembl_id and symbol as separate columns
    ## so plotting code can strip the prefix at render time.
    gene_ids <- as.character(var_df[["_index"]])
    parts <- strsplit(gene_ids, "|", fixed = TRUE)
    var_df$ensembl_id <- vapply(
        parts,
        function(p) p[[1]],
        character(1)
    )
    var_df$symbol <- vapply(
        parts,
        function(p) if (length(p) >= 2) p[[2]] else NA_character_,
        character(1)
    )
    rownames(counts) <- gene_ids
    rownames(var_df) <- gene_ids

    if (verbose) message("Reading obs (cell metadata)...")
    obs_df <- read_h5_table(path, "/obs")
    cell_ids <- as.character(obs_df[["_index"]])
    if (length(cell_ids) == ncol(counts)) {
        colnames(counts) <- cell_ids
        rownames(obs_df) <- cell_ids
    }

    list(counts = counts, obs = obs_df, var = var_df)
}

#' Read an AnnData /obs or /var group into a data.frame.
#' Handles legacy categorical encoding (/obs/__categories/<field> integer codes)
#' and modern encoding (/obs/<field>/categories + /codes).
#' @noRd
read_h5_table <- function(path, group) {
    all_nodes <- rhdf5::h5ls(path)
    children <- all_nodes[all_nodes$group == group, ]
    legacy_cats_group <- file.path(group, "__categories")
    legacy_cats <- all_nodes[all_nodes$group == legacy_cats_group, "name"]

    cols <- list()
    for (i in seq_len(nrow(children))) {
        name <- children$name[i]
        otype <- children$otype[i]
        if (name == "__categories") next
        full <- file.path(group, name)

        if (otype == "H5I_DATASET") {
            vals <- rhdf5::h5read(path, full)
            ## Legacy categorical: integer codes paired with /obs/__categories/<name>
            if (name %in% legacy_cats) {
                levels_str <- rhdf5::h5read(
                    path, file.path(legacy_cats_group, name)
                )
                vals <- as.character(levels_str)[as.integer(vals) + 1L]
            } else if (is.array(vals) && length(dim(vals)) > 1) {
                ## Skip multi-dim datasets (e.g. embeddings stored under obsm-style)
                next
            }
            cols[[name]] <- as.vector(vals)
        } else if (otype == "H5I_GROUP") {
            ## Modern categorical encoding: <name>/categories + <name>/codes
            sub <- all_nodes[all_nodes$group == full, "name"]
            if (all(c("categories", "codes") %in% sub)) {
                lv <- rhdf5::h5read(path, file.path(full, "categories"))
                cd <- rhdf5::h5read(path, file.path(full, "codes"))
                cd[cd < 0] <- NA_integer_
                cols[[name]] <- as.character(lv)[as.integer(cd) + 1L]
            }
        }
    }

    n <- max(vapply(cols, length, integer(1)))
    cols <- cols[vapply(cols, length, integer(1)) == n]
    as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
}
