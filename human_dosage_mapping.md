# Moving the human dosage sensitivity scores onto the T2T reference

This is a record of what `ref_converter.R` does to the Collins et al. (2022) dosage sensitivity table, and of the decisions taken along the way. It is meant to be read alongside the script, and to explain the choices that the code alone does not.

## What the script does

The dosage table gives two scores per gene, `pHaplo` (how likely the gene is to cause harm when one copy is lost) and `pTriplo` (the same for a gained copy). It lists 18,641 genes, identified only by name, and those names come from an annotation of the older GRCh37 reference genome.

The script works out which gene in the T2T-CHM13v2.0 annotation each of those names refers to. It then writes a new table carrying the original scores next to the T2T gene's name, identifiers and position. The original table is not modified.

Of the 18,641 genes, 18,028 (96.7%) were matched to a T2T gene, 4 could not be settled between two candidates, and 609 (3.3%) could not be matched at all.

## The files involved

| File | What it is |
|---|---|
| `human/Collins_rCNV_2022.dosage_sensitivity_scores.tsv` | The input: gene name, `pHaplo`, `pTriplo`. |
| `human/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff` | NCBI's annotation of the T2T genome (release RS_2025_08). This is the target. |
| `human/hgnc_master.tsv` | The HGNC table of official human gene names, with each gene's earlier names, alternative names and Ensembl identifier. |
| `human/GRCh37/gencode.v19.annotation.gtf.gz` | The GENCODE release 19 annotation of GRCh37. Downloaded on 11 October 2026 from the address noted in the script. |
| `human/Collins_rCNV_2022.dosage_sensitivity_scores.T2T.tsv` | The output. |

## Why this is not a simple lookup

Three different groups name human genes, and the three files above come from all three.

The dosage table uses names from GENCODE, which is the Ensembl project's annotation. The T2T file comes from NCBI, which maintains its own, separate annotation. HGNC is the committee that assigns official names, and both of the others follow it, but not at the same pace and not for every gene.

Two things follow from this. First, a gene's name can change over time, so a name that was current when the dosage table was built may be listed under something else today. Second, GENCODE gives a temporary name to any gene that has no official name yet. These temporary names are built from the label of the piece of sequenced DNA the gene was found on, for example `RP11-156P1.2` or `AC120194.1`. NCBI never uses names of this kind, and HGNC does not record them, so they cannot be found by name anywhere else.

We confirmed that the dosage table was built on GENCODE release 19: every one of its 18,641 names appears in that release.

## How a gene is matched

Each gene name is tried against five sources, in a fixed order. The first source that recognises the name decides the answer.

1. **Ensembl identifier.** GENCODE release 19 gives every gene name a permanent Ensembl identifier. HGNC records which official gene each Ensembl identifier belongs to, and the T2T file records the HGNC identifier of each of its genes. Following that chain leads from the old name to the T2T gene without relying on the name itself.
2. **Same name.** The name is identical to a gene name in the T2T file.
3. **Earlier official name.** HGNC lists the name as a previous official name of a gene.
4. **Alternative name.** HGNC lists the name as an informal alternative for a gene.
5. **Synonym in the T2T file.** The T2T file itself lists the name as a synonym of one of its genes.

| Source that decided | Genes matched | Left undecided |
|---|---|---|
| Ensembl identifier | 17,699 | – |
| Same name | 214 | – |
| Earlier official name | 101 | 1 |
| Alternative name | 9 | – |
| Synonym in the T2T file | 5 | 3 |
| None | 609 unmatched | |

Once a T2T gene name is settled, the script attaches that gene's details from the T2T file: its record identifier, its NCBI gene number, its HGNC identifier, its position and its type.

## Decisions

### The Ensembl identifier is tried first

The first version of the script tried names first and had no Ensembl step at all. The Ensembl step was added later to recover the temporary names, and was then moved to the front.

The reason for moving it is that names get reused. In GRCh37 the name `NAA38` belonged to a gene on chromosome 7 that is now called `LSM8`. Today `NAA38` is the name of a different gene on chromosome 17, which used to be called `LSMD1`. Matching by name put the old `NAA38` scores on the wrong gene. The identifier does not have this problem, because it stays with the gene when the name changes.

This was not a clear-cut win, and it should not be read as one. The two approaches agree for 17,509 genes and disagree for only 24. For those 24 we checked which candidate actually sits where the original gene sat on GRCh37. The identifier's answer was the better fit for 8, the name's answer for 8, and neither was clearly better for the remaining 8. The identifier tends to be right where a name was reused, and wrong for records that span two neighbouring genes (for example `TMEM256-PLSCR3`) and for near-identical gene copies (for example `POTEM` and `POTEG`).

Since neither order fixes those 24, they are marked in the output instead: the column `symbol_based_alt` holds the answer the name-based steps would have given. Any gene with a value in that column should be checked by hand before its scores are used.

### A later source can break a tie but never overrule

Sometimes a source offers two candidates for one name. GENCODE release 19, for instance, has two genes both called `LSP1`: the real one on chromosome 11 and a short fragment on chromosome 13.

When that happens, the script looks at the later sources. If one of them points to exactly one of the two candidates, that candidate is chosen. Eight genes were settled this way, and the output column `tiebreak` says which source did it. A later source is never allowed to replace an answer with something the earlier source did not offer.

### Genes that stay tied are left unmatched

Four genes could not be settled: `ATP6C`, `TP53TG3C`, `POTEB` and `TBC1D3C`. Rather than pick one candidate arbitrarily, the script leaves them without a T2T gene, marks them `ambiguous`, and lists the candidates.

### HGNC names are followed by identifier, not by name

When HGNC says an old name now belongs to a certain gene, the script finds that gene in the T2T file by its HGNC identifier rather than by its current official name. The T2T file and the HGNC table were made at different times, and they disagree on the current name of 115 genes.

### Records marked as pseudogenes are included

The T2T file keeps working genes and pseudogenes (broken or inactive copies of genes) as two separate kinds of record. Both are loaded. If only working genes were loaded, a name that belongs to a pseudogene in T2T would find no match under its own name and could then be matched, wrongly, to a different gene that happens to list it as an alternative name. In the output, 115 rows point to pseudogene records. The `gene_biotype` column makes them easy to filter out.

### Capitalisation is ignored

A few names in the dosage table differ from the T2T names only in capitalisation, such as `C2ORF15` against `C2orf15`. Names are compared without regard to case.

### Two old genes landing on one T2T gene are both kept

Some genes that were separate in GRCh37 have since been merged, so two rows of the dosage table can point to the same T2T gene. This affects 101 dosage genes, which land on 50 T2T genes. The script keeps every row with its own scores and does not combine them. The column `n_sharing_target` shows how many rows share each T2T gene. How to combine the scores, if at all, is left to whoever uses the table.

### One name with two T2T records gives two rows

The name `CLN3` appears on two records in the T2T file: the real gene and a 253-base fragment. Both are kept, so the output has 18,642 rows for 18,641 genes. The column `n_t2t_entries` marks this.

### Which identifiers are attached

The request was to attach RefSeq identifiers from the gene records of the T2T file. Those records carry a record identifier (such as `gene-CACNA1C`) and an NCBI gene number, and these are what the output contains. The transcript-level identifiers that begin with `NM_` or `NR_` are stored on different records and are not included.

### Things tried and left out

Two further steps were tested and dropped because they matched nothing new: following HGNC's notes on withdrawn names, and looking names up among HGNC's current official names.

Matching genes by their position on the genome, rather than by name or identifier, was considered as a way to reach the unmatched genes. It was not done. It would need an extra file describing how GRCh37 positions correspond to T2T positions, and a match by position is less certain than a match by identifier.

## What could not be matched

Of the 609 unmatched genes, 544 have temporary GENCODE names of the kind described above. The remaining 65 are a mix of old informal names and official names.

Looked at another way, 573 of the 609 have an Ensembl identifier that HGNC does not list. Either the gene never received an official name or it has since been withdrawn. The other 36 do have an HGNC record, but the T2T file has no gene carrying that HGNC identifier. `GSTM1`, `CFHR1` and `CFHR3` are examples.

## Where to be cautious

- **The 24 genes with a value in `symbol_based_alt`.** The two approaches disagree and neither is reliably right.
- **Matches made by alternative name or by T2T synonym.** These 14 genes rest on the weakest evidence. One of them, `OR2A7`, was matched to `OR2A4`, which sits on a different chromosome in T2T.
- **`MIR3654`.** Matched by identical name, but it is on a different chromosome in T2T than in GRCh37.
- **`MEGT1` and `LST3`.** Both were matched, by Ensembl identifier, to records that span two neighbouring genes.
- **Gene type.** The dosage table covers protein-coding genes, but 303 of the matched rows are now classed as long non-coding RNA in T2T and 115 as pseudogenes.

Apart from `MIR3654` and `OR2A7`, every matched gene sits on the same chromosome in T2T as it did in GRCh37. The scores themselves were checked to be unchanged from the input.

## Columns in the output

| Column | Meaning |
|---|---|
| `grch37_symbol` | Gene name as given in the dosage table. |
| `t2t_symbol` | Name of the matched gene in the T2T file. Empty if not matched. |
| `t2t_gff_id` | Record identifier of that gene in the T2T file. |
| `ncbi_gene_id` | NCBI gene number. |
| `hgnc_id` | HGNC identifier. |
| `pHaplo`, `pTriplo` | The original scores, unchanged. |
| `seqid`, `start`, `end`, `strand` | Position of the gene on T2T. |
| `gene_biotype` | Gene type according to the T2T file. |
| `status` | `mapped`, `ambiguous` or `unmapped`. |
| `method` | Which of the five sources decided the match. |
| `tiebreak` | The later source that settled a tie, if there was one. |
| `n_candidates`, `candidates` | How many T2T genes the deciding source offered, and which. |
| `symbol_based_alt` | A different answer given by the name-based sources. Check these by hand. |
| `n_sharing_target` | How many dosage genes were matched to this same T2T gene. |
| `n_t2t_entries` | How many records in the T2T file carry this gene name. |

## Running it again

Run `Rscript ref_converter.R` from the `reference_data` directory. It takes about ten seconds and overwrites the output file. It needs the R packages `rtracklayer`, `Rsamtools` and `tidyverse`.

The first half of the script holds general-purpose functions that know nothing about humans or about these particular files. To do the same job for another species, the part to write is a new list of sources in the order they should be tried; the matching itself is reused as it is.
