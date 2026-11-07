library(rtracklayer)
library(Rsamtools)
library(tidyverse)

# =============================================================================
# Functions (generic, no species or reference specifics)
# =============================================================================

# GFF/GTF features of the given types as a tibble (multi-valued attributes
# become list columns)
read_gff_features <- function(gff_path, feature_types = "gene") {
  gr <- import(gff_path, feature.type = feature_types)
  attrs <- mcols(gr)[setdiff(names(mcols(gr)), c("score", "phase"))]
  attrs <- lapply(attrs, \(col) {
    if (is(col, "List")) as.list(col) else as.character(col)
  })
  tibble(
    seqid = as.character(seqnames(gr)),
    start = start(gr),
    end = end(gr),
    strand = as.character(strand(gr))
  ) |>
    bind_cols(as_tibble(attrs))
}

# e.g. extract_dbxref(genes$Dbxref, "GeneID") -> "81399"
extract_dbxref <- function(dbxref, db) {
  prefix <- paste0(db, ":")
  map_chr(dbxref, \(x) {
    hit <- x[startsWith(x, prefix)]
    if (length(hit) == 0) NA_character_ else str_remove(hit[1], fixed(prefix))
  })
}

# (from, to) lookup, one row per value of a multi-valued `from`
make_lookup <- function(tbl, from, to, sep = ",\\s*") {
  out <- tbl |> transmute(from = {{ from }}, to = {{ to }})
  if (is.character(out$from)) out <- out |> mutate(from = str_split(from, sep))
  out |>
    unnest(from) |>
    mutate(from = str_trim(from)) |>
    filter(!is.na(from), from != "", !is.na(to)) |>
    distinct()
}

# (x -> y) then (y -> z) gives (x -> z)
chain_lookups <- function(a, b) {
  inner_join(
    a, b,
    by = join_by(to == from),
    suffix = c("", "_next"),
    relationship = "many-to-many"
  ) |>
    transmute(from, to = to_next) |>
    distinct()
}

# Map IDs through an ordered, named list of lookups. The first tier with a hit
# decides; later tiers can only break a tie between its candidates.
map_tiered <- function(query, tiers, normalize = identity) {
  queries <- tibble(query = unique(query)) |> mutate(key = normalize(query))
  hits <- imap(tiers, \(lookup, tier) {
    lookup |>
      mutate(from = normalize(from)) |>
      inner_join(queries, by = join_by(from == key), relationship = "many-to-many") |>
      distinct(query, to) |>
      mutate(rank = match(tier, names(tiers)))
  }) |>
    list_rbind()
  first <- hits |> filter(rank == min(rank), .by = query)
  tiebreaks <- first |>
    filter(n() > 1, .by = query) |>
    rename(first_rank = rank) |>
    inner_join(hits, by = join_by(query, to)) |>
    filter(rank > first_rank) |>
    filter(n() == 1, .by = c(query, rank)) |>
    slice_min(rank, by = query) |>
    transmute(query, tiebreak_target = to, tiebreak = names(tiers)[rank])
  queries |>
    select(query) |>
    left_join(
      first |>
        summarise(
          method = names(tiers)[rank[1]],
          n_candidates = n(),
          candidates = str_flatten(sort(to), ";"),
          .by = query
        ),
      by = "query"
    ) |>
    left_join(tiebreaks, by = "query") |>
    mutate(
      n_candidates = replace_na(n_candidates, 0L),
      target = if_else(n_candidates == 1, candidates, tiebreak_target),
      status = case_when(
        !is.na(target) ~ "mapped",
        n_candidates > 1 ~ "ambiguous",
        .default = "unmapped"
      ),
      method = factor(method, levels = names(tiers))
    ) |>
    mutate(n_sharing_target = if_else(is.na(target), NA_integer_, n()), .by = target) |>
    select(query, target, status, method, tiebreak, n_candidates, candidates, n_sharing_target)
}

summarise_mapping <- function(mapping) {
  mapping |>
    count(method, status) |>
    mutate(pct = round(100 * n / sum(n), 2))
}

# =============================================================================
# Human: GRCh37 gene symbols -> T2T-CHM13v2.0
# Decisions are written up in human_dosage_mapping.md
# =============================================================================

human_dir <- "human"
t2t_gff_path <- file.path(human_dir, "T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff")
hgnc_path <- file.path(human_dir, "hgnc_master.tsv")
dosage_path <- file.path(human_dir, "Collins_rCNV_2022.dosage_sensitivity_scores.tsv")
dosage_t2t_path <- file.path(human_dir, "Collins_rCNV_2022.dosage_sensitivity_scores.T2T.tsv")
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gencode_grch37_path <- file.path(human_dir, "GRCh37/gencode.v19.annotation.gtf.gz")

# ---- Load -------------------------------------------------------------------

t2t_genes <- read_gff_features(t2t_gff_path, c("gene", "pseudogene")) |>
  mutate(
    ncbi_gene_id = extract_dbxref(Dbxref, "GeneID"),
    hgnc_id = extract_dbxref(Dbxref, "HGNC")
  )

hgnc <- read_tsv(hgnc_path, quote = "", na = "", col_types = cols(.default = "c")) |>
  rename_with(\(x) str_replace_all(str_to_lower(x), " ", "_"))

dosage <- read_tsv(dosage_path, col_types = "cdd") |>
  rename(gene = `#gene`)

# The dosage table uses GENCODE v19 gene names
gencode_grch37_genes <- read_gff_features(gencode_grch37_path, "gene") |>
  mutate(ensembl_gene_id = str_remove(gene_id, "\\.\\d+$"))

# ---- Lookups ----------------------------------------------------------------

hgnc_id_to_t2t <- make_lookup(t2t_genes, hgnc_id, Name)
hgnc_approved <- hgnc |> filter(status == "Approved")

# In priority order
tiers <- list(
  ensembl_id = make_lookup(gencode_grch37_genes, gene_name, ensembl_gene_id) |>
    chain_lookups(make_lookup(hgnc_approved, ensembl_gene_id, hgnc_id)) |>
    chain_lookups(hgnc_id_to_t2t),
  direct = make_lookup(t2t_genes, Name, Name),
  hgnc_previous = make_lookup(hgnc_approved, previous_symbols, hgnc_id) |>
    chain_lookups(hgnc_id_to_t2t),
  hgnc_alias = make_lookup(hgnc_approved, alias_symbols, hgnc_id) |>
    chain_lookups(hgnc_id_to_t2t),
  gff_synonym = make_lookup(t2t_genes, gene_synonym, Name)
)

# ---- Map symbols, then attach the T2T gene entries ---------------------------

mapping <- map_tiered(dosage$gene, tiers, normalize = toupper)
summarise_mapping(mapping) |> print()

# Flag genes where the symbol tiers alone would have picked a different target
symbol_mapping <- map_tiered(dosage$gene, tiers[names(tiers) != "ensembl_id"], normalize = toupper)
mapping <- mapping |>
  left_join(select(symbol_mapping, query, symbol_target = target), by = "query") |>
  mutate(symbol_based_alt = if_else(symbol_target != target, symbol_target, NA_character_)) |>
  select(-symbol_target)

# One row per T2T gene entry, so a symbol with several entries is repeated
dosage_t2t <- dosage |>
  left_join(mapping, by = join_by(gene == query)) |>
  left_join(
    t2t_genes |>
      select(
        t2t_symbol = Name, t2t_gff_id = ID, ncbi_gene_id, hgnc_id,
        seqid, start, end, strand, gene_biotype
      ) |>
      add_count(t2t_symbol, name = "n_t2t_entries"),
    by = join_by(target == t2t_symbol),
    relationship = "many-to-many"
  ) |>
  select(
    grch37_symbol = gene, t2t_symbol = target, t2t_gff_id, ncbi_gene_id, hgnc_id,
    pHaplo, pTriplo,
    seqid, start, end, strand, gene_biotype,
    status, method, tiebreak, n_candidates, candidates, symbol_based_alt,
    n_sharing_target, n_t2t_entries
  )

write_tsv(dosage_t2t, dosage_t2t_path, na = "")
