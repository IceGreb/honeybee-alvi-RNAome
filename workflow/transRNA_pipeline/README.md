# transRNA Pipeline v1.2

Paired-end RNA-seq pipeline: raw reads → filtering → taxonomy → plots.  
Datasets: **RJ** (Royal Jelly) and **ST** (Systemic larval tissue), *Apis mellifera*.

---

## How to run

```bash
# From the CSD3 login node — do NOT sbatch this:
bash run.sh           # fresh run from scratch
bash run.sh --resume  # resume after failure (Nextflow skips completed tasks)
```

To change partition or project, edit the top two lines of `run.sh`:
```bash
PARTITION="icelake"
PROJECT="<your-slurm-project>"
```

---

## Samples

| Sample | Group | Sample | Group |
|--------|-------|--------|-------|
| RJ1    | RJ    | T1GMN  | ST    |
| RJ2    | RJ    | T4GMN  | ST    |
| RJ3    | RJ    | TA1    | ST    |
|        |       | TB1    | ST    |
|        |       | TC1    | ST    |
|        |       | TD1    | ST    |

Edit `samples.csv` to add or remove samples.

---

## Pipeline steps and outputs

```
Step 1   Raw reads (fq.gz)
           → reports/01_raw/{sample}_raw_stats.tsv

Step 2   Trimmed reads (TrimGalore, pre-computed)
           → reports/02_trimmed/{sample}_trimmed_stats.tsv
           → reports/reads_posttrim_tab.tsv            ← produced fresh here

Step 3   STAR unmapped reads (pre-computed)
           → reports/03_star/{sample}_star_unmapped_stats.tsv
           → reports/03_star/{sample}_Log.final.out

Step 4a  BBsplit host/human filtered (pre-computed)
           → reports/04a_bbsplit_host/{sample}_bbsplit_host_stats.tsv

Step 4b  BBsplit antiviral filtered (pre-computed)
           → reports/04b_bbsplit_virus/{sample}_bbsplit_virus_stats.tsv

Step 5   BBsplit MAG-clean reads (pre-computed)
           → reports/05_no_mags/{sample}_noMAGs_stats.tsv

Step 6   KRAKEN annotation + filtering
           annotate_kraken_lineage.py  →  taxonkit appends lineage col 6
           filter_kraken_invertebrates.py:
             - KEEPS C (classified) rows only
             - U rows silently skipped → handled by BLAST path (not lost)
             - Removes: viruses, invertebrates, host genera, length < 18 nt
           → kraken_annotated/{sample}_kraken_annotated.tsv
           → kraken_filtered/{sample}_kraken_invertebrates_filtered.tsv

Step 7   BLAST annotation + filtering (mates 1 and 2 independently)
           annotate_blast_lineage.py   →  taxonkit appends 8 taxonomy cols
           filter_blast_all_conditions.py:
             - Removes: viruses, invertebrates, host species
             - Removes: mismatch≠0, gapopen≠0, not-full-span, length<18 nt
           → blast_annotated/{sample}_{1,2}_blast_annotated.tsv
           → blast_filtered/{sample}_{1,2}_blast_all_lengths_filtered.tsv

Step 8   COLLAPSE (extract from trimmed reads → merge → seqkit rmdup)
           Passing IDs from BLAST and Kraken filtered TSVs are used to fetch
           sequences from trimmed FASTQ files. Mate 1 and mate 2 BLAST-passing
           reads are fetched separately; Kraken-passing IDs are fetched from
           mate 1 only. All three are cat-merged and collapsed with seqkit rmdup.
           → collapsed/{sample}_merged.fq
           → collapsed/{sample}_merged_collapsed_clean.fq
           → collapsed/{sample}_merged_duplicated.detail.txt
           → collapsed/{sample}_collapse_stats.tsv

Step 9   CANDIDATE SELECTION (ge5_18nt_filtering_pipeline_for_both_mates.py)
           Candidate = ≥5 duplicate occurrences AND length ≥ 18 nt.
           Lengths read from collapsed_clean.fq; weights from detail.txt.
           → candidates/{sample}_ge5_detected_ids.txt
           → candidates/{sample}_ge5_detected.hist          (weighted by dup count)
           → candidates/{sample}_ge5_detected_weighted_ids.tsv

Step 10  EXTRACT FASTA (seqkit grep + fq2fa on collapsed_clean.fq)
           → final_transRNAs_fasta/{sample}_final_transRNAs.fasta

Step 11  LENGTH HISTOGRAMS + PLOT
           → final_transRNAs_length_hists/{sample}_ge5_detected.hist
           → plots/{date}_avg_reinflated_length_distribution_15to100nt.png

Step 12  TAXONOMY PARSER  (called ONCE GLOBALLY — all samples together)
           04_05_2026_transRNA_taxonomy_parser.py:
             - Re-inflates using duplicate weights
             - RPM-normalises using reads_posttrim_tab.tsv
             - Priority: blast (configurable)
             - Produces per-sample AND per-group (RJ/ST) tables
           → final_transRNAs_taxonomies/{sample}_{blast,kraken,combined}_{Rank}_summary.tsv
           → final_transRNAs_taxonomies/{sample}_{blast,kraken,combined}_{Rank}_top10.tsv
           → final_transRNAs_taxonomies/{RJ,ST}_{blast,kraken,combined}_{Rank}_summary.tsv
           → final_transRNAs_taxonomies/{RJ,ST}_{blast,kraken,combined}_{Rank}_top10.tsv
           → final_transRNAs_taxonomies/fetch_reinflate_report.tsv

Step 13  TAXONOMY PLOT (plot_top10_taxa_global_colors.py)
           Reads RJ_*_top10.tsv and ST_*_top10.tsv files.
           Multi-panel broken-axis horizontal bar plot; ranks: Domain/Kingdom/Order/Species
           → plots/*.png  (RJ and ST side-by-side panels)

Step 14  AGGREGATE READ-COUNT REPORT
           → reports/pipeline_read_counts_report.tsv
           → reports/pipeline_read_counts_report.html
           → reports/dataset_summary_report.tsv
           → reports/virus_exclusion_report.tsv
```

---

## Filtering criteria — complete specification

### Kraken (`filter_kraken_invertebrates.py`)

| | Criterion | How checked |
|---|---|---|
| **SKIP** | Unclassified reads (U) | C/U flag col 1 — these go to BLAST path |
| **REMOVE** | Either mate < 18 nt | length field `"35\|35"`, both sides checked |
| **REMOVE** | Virus/acellular domain | Lineage[0] (Domain rank) |
| **REMOVE** | Host genera | Lineage[6] (Genus) ∈ {Homo, Mus, Canis, Felis, Apis} |
| **REMOVE** | Invertebrate phyla | Regex on full lineage vs invertebrate_phyla.txt names |

Lineage appended as col 6: `Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species`  
(taxonkit `{k};{K};{p};{c};{o};{f};{g};{s}` with `--fill-miss-rank`)

### BLAST (`filter_blast_all_conditions.py`)

| Criterion | Field |
|---|---|
| Virus/acellular domain | last-8 col 0 (Domain) |
| Invertebrate phylum | last-8 col 2 (Phylum), exact match |
| Host species | last-8 col 7 (Species), substring: homo/mus/canis lupus/felis/apis mellifera |
| mismatch ≠ 0 | col 5 |
| gapopen ≠ 0 | col 6 |
| Not full-span: abs(qend−qstart)+1 ≠ length | cols 4,7,8 |
| Alignment length < 18 nt | col 4 |

---

## Taxonomy lineage — both tools

```bash
echo "$taxid" \
  | taxonkit lineage --data-dir $DB \
  | taxonkit reformat --data-dir $DB \
      --format "{k};{K};{p};{c};{o};{f};{g};{s}" \
      --fill-miss-rank --miss-taxid-repl unassigned
```

8 ranks always present: **Domain ; Kingdom ; Phylum ; Class ; Order ; Family ; Genus ; Species**  
Kingdom is rank index 1. All 8 ranks appear in all taxonomy tables and plots.

---

## taxonkit setup (once only)

```bash
mkdir -p ~/.taxonkit && cd ~/.taxonkit
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
echo "9606" | taxonkit lineage   # test: should print Homo sapiens lineage
```

---

## Software requirements

```
python/3.11  (+ pip install pandas numpy matplotlib seaborn biopython)
seqkit/2.8.0
taxonkit/0.17.0
nextflow/24.04.4
```

---

## HPC tips (CSD3)

- Set `workDir` in `nextflow.config` to a path outside your home quota (scratch or RDS)
- `process.cache = 'lenient'` prevents spurious cache misses from RDS NFS jitter
- For very large BLAST TSVs, switch `ANNOTATE_BLAST` to `icelake-himem` in `cambridge.config`
- Monitor jobs: `squeue -u <your_username>` and `tail -f logs/nextflow_*.log`
