"""
Aim 1: Computational Survey of Dietary Interventions in Brain Transcriptomics
Queries NCBI GEO for RNA-seq datasets involving dietary interventions in brain tissue.
Searches per-condition to stay within URL limits, then aggregates.
"""

import time
import re
import pandas as pd
from Bio import Entrez

Entrez.email = "luutrungduc2909@gmail.com"

# ── Search vocabulary ─────────────────────────────────────────────────────────

# Each entry is (short_key, full label, GEO search term)
DIET_QUERIES = [
    ("ketogenic",         "Ketogenic diet",           '"ketogenic diet"'),
    ("high_fat",          "High fat diet",             '"high fat diet"'),
    ("caloric_restrict",  "Caloric restriction",       '"caloric restriction"'),
    ("intermit_fasting",  "Intermittent fasting",      '"intermittent fasting"'),
    ("fasting",           "Fasting",                   '"fasting"'),
    ("choline_def",       "Choline deficiency",        '"choline deficien"'),
    ("vit_d_def",         "Vitamin D deficiency",      '"vitamin D deficien"'),
    ("zinc_def",          "Zinc deficiency",           '"zinc deficien"'),
    ("iron_def",          "Iron deficiency",           '"iron deficien"'),
    ("folate_def",        "Folate deficiency",         '"folate deficien"'),
    ("b12_def",           "Vitamin B12 deficiency",    '"vitamin B12"'),
    ("magnesium_def",     "Magnesium deficiency",      '"magnesium deficien"'),
    ("selenium_def",      "Selenium deficiency",       '"selenium deficien"'),
    ("fructose",          "Fructose",                  '"fructose"'),
    ("sucrose",           "Sucrose / high sugar",      '"sucrose"'),
    ("omega3_dha",        "Omega-3 / DHA",             '"omega-3" OR "DHA" OR "docosahexaenoic"'),
    ("omega3_epa",        "EPA",                       '"eicosapentaenoic"'),
    ("western_diet",      "Western diet",              '"Western diet"'),
    ("mediterranean",     "Mediterranean diet",        '"Mediterranean diet"'),
    ("protein_restrict",  "Protein restriction",       '"protein restriction"'),
    ("methionine_restr",  "Methionine restriction",    '"methionine restriction"'),
    ("high_salt",         "High salt diet",            '"high salt diet" OR "high sodium"'),
    ("resveratrol",       "Resveratrol",               '"resveratrol"'),
    ("curcumin",          "Curcumin",                  '"curcumin"'),
]

BRAIN_REGION_TERMS = (
    '"hippocampus" OR "hypothalamus" OR "cortex" OR "cerebral cortex" OR '
    '"prefrontal cortex" OR "cerebellum" OR "striatum" OR "amygdala" OR '
    '"brainstem" OR "whole brain" OR "brain tissue" OR "frontal cortex" OR '
    '"olfactory bulb" OR "substantia nigra"'
)

ASSAY_TERMS = '"RNA-seq" OR "RNAseq" OR "transcriptome" OR "gene expression" OR "microarray"'

# ── GEO helpers ───────────────────────────────────────────────────────────────

def search_geo(query, retmax=300):
    try:
        handle = Entrez.esearch(db="gds", term=query, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.35)
        return record.get("IdList", []), int(record.get("Count", 0))
    except Exception as e:
        print(f"    [search error: {e}]")
        return [], 0


def fetch_summaries(uid_list, batch=200):
    results = []
    for i in range(0, len(uid_list), batch):
        chunk = uid_list[i:i+batch]
        try:
            handle = Entrez.esummary(db="gds", id=",".join(chunk))
            records = Entrez.read(handle)
            handle.close()
            results.extend(records)
            time.sleep(0.35)
        except Exception as e:
            print(f"    [fetch error chunk {i}: {e}]")
    return results


# ── Metadata parser ───────────────────────────────────────────────────────────

SPECIES_MAP = {
    "mouse":  ["mus musculus", "mouse", "mice", "c57bl", "balb"],
    "rat":    ["rattus norvegicus", "rat ", "sprague", "wistar"],
    "human":  ["homo sapiens", "human", "patient", "subject"],
}

SEX_MAP = {
    "male":   [" male", "males"],
    "female": ["female", "females"],
    "both":   ["both sexes", "male and female", "female and male"],
}

REGION_LIST = [
    "hippocampus", "hypothalamus", "prefrontal cortex", "frontal cortex",
    "temporal cortex", "cerebral cortex", "cortex", "cerebellum",
    "striatum", "amygdala", "brainstem", "olfactory bulb",
    "substantia nigra", "nucleus accumbens", "whole brain",
]

def detect(text, vocab_map):
    t = text.lower()
    for label, kws in vocab_map.items():
        if any(k in t for k in kws):
            return label
    return "unknown"

def detect_list(text, items):
    t = text.lower()
    return [item for item in items if item.lower() in t]


def parse_record(rec, diet_key, diet_label):
    title    = str(rec.get("title", ""))
    summary  = str(rec.get("summary", ""))
    organism = str(rec.get("taxon", ""))
    full     = (title + " " + summary + " " + organism).lower()

    regions  = detect_list(full, REGION_LIST)
    species  = detect(full, SPECIES_MAP)
    sex      = detect(full, SEX_MAP)
    year     = str(rec.get("PDAT", ""))[:4]
    gds_type = str(rec.get("gdsType", ""))
    is_rnaseq = any(k in (full + gds_type.lower())
                    for k in ["rna-seq", "rnaseq", "transcriptom", "mrna-seq"])

    return {
        "accession":    str(rec.get("Accession", "")),
        "diet_key":     diet_key,
        "diet_label":   diet_label,
        "title":        title,
        "organism":     organism,
        "species":      species,
        "n_samples":    int(rec.get("n_samples", 0)),
        "regions":      "; ".join(regions),
        "sex":          sex,
        "year":         year,
        "assay_rnaseq": is_rnaseq,
        "gds_type":     gds_type,
        "summary":      summary[:400],
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def run_survey():
    print("=" * 65)
    print("Aim 1 — GEO Brain Dietary Transcriptomics Survey")
    print("=" * 65)

    all_rows = []
    seen_uids = set()

    for diet_key, diet_label, diet_term in DIET_QUERIES:
        query = f"({diet_term}) AND ({BRAIN_REGION_TERMS}) AND ({ASSAY_TERMS})"
        uids, total = search_geo(query, retmax=300)

        # de-duplicate across conditions
        new_uids = [u for u in uids if u not in seen_uids]
        seen_uids.update(new_uids)

        print(f"  {diet_label:<35} GEO hits: {total:>4}  |  new unique: {len(new_uids)}")

        if not new_uids:
            # still record the 0-hit condition
            all_rows.append({
                "accession": "", "diet_key": diet_key, "diet_label": diet_label,
                "title": "", "organism": "", "species": "", "n_samples": 0,
                "regions": "", "sex": "unknown", "year": "", "assay_rnaseq": False,
                "gds_type": "", "summary": "", "_total_hits": total,
            })
            continue

        recs = fetch_summaries(new_uids)
        for rec in recs:
            row = parse_record(rec, diet_key, diet_label)
            row["_total_hits"] = total
            all_rows.append(row)

    df = pd.DataFrame(all_rows)
    df_data = df[df["accession"] != ""].copy()

    # Save
    df_data.to_csv("/home/duc/letters/aim1_brain_diet_filtered.csv", index=False)
    print(f"\nTotal unique records fetched: {len(df_data)}")

    # ── Landscape Summary ─────────────────────────────────────────────────────

    print("\n" + "=" * 65)
    print("LANDSCAPE MAP")
    print("=" * 65)

    # 1. Hits per dietary condition
    print("\n-- Datasets per dietary condition (brain + assay filter) --")
    hit_summary = {}
    for _, diet_key, diet_label, diet_term in [(r["diet_key"],r["diet_key"],r["diet_label"],r.get("summary","")) for r in all_rows]:
        pass
    # cleaner approach
    for _, row in df.iterrows():
        k = row["diet_label"]
        if k not in hit_summary:
            hit_summary[k] = {"n_records": 0, "total_hits": int(row.get("_total_hits", 0))}
        if row["accession"]:
            hit_summary[k]["n_records"] += 1

    print(f"  {'Condition':<38} {'GEO hits':>9}  {'Records fetched':>15}")
    print("  " + "-" * 65)
    for _, diet_label, _ in DIET_QUERIES:
        info = hit_summary.get(diet_label, {"n_records": 0, "total_hits": 0})
        gap_flag = "  <<< CANDIDATE GAP" if info["total_hits"] == 0 else ""
        print(f"  {diet_label:<38} {info['total_hits']:>9}  {info['n_records']:>15}{gap_flag}")

    # 2. Brain regions
    print("\n-- Brain regions covered across all filtered records --")
    region_counts = {}
    for _, row in df_data.iterrows():
        for r in row["regions"].split("; "):
            r = r.strip()
            if r:
                region_counts[r] = region_counts.get(r, 0) + 1
    for region, n in sorted(region_counts.items(), key=lambda x: -x[1]):
        print(f"  {region:<40} {n}")

    # 3. Species
    print("\n-- Species --")
    print(df_data["species"].value_counts().to_string())

    # 4. Sex
    print("\n-- Sex representation --")
    print(df_data["sex"].value_counts().to_string())

    # 5. Assay type
    print("\n-- Assay type (RNA-seq vs other) --")
    print(df_data["assay_rnaseq"].value_counts().rename(
        {True: "RNA-seq/transcriptome", False: "other/unspecified"}).to_string())

    # 6. Year
    print("\n-- Year distribution --")
    print(df_data[df_data["year"] != ""]["year"].value_counts().sort_index().to_string())

    # ── Gap Analysis ──────────────────────────────────────────────────────────

    print("\n" + "=" * 65)
    print("GAP ANALYSIS")
    print("=" * 65)
    priority = ["Choline deficiency", "Vitamin D deficiency", "Zinc deficiency",
                "Iron deficiency", "Folate deficiency", "Vitamin B12 deficiency",
                "Magnesium deficiency", "Selenium deficiency",
                "Ketogenic diet", "Intermittent fasting"]
    for label in priority:
        n = hit_summary.get(label, {}).get("total_hits", 0)
        status = "NO DATA" if n == 0 else f"{n} datasets"
        print(f"  {label:<38} {status}")

    # ── Dataset Listing ───────────────────────────────────────────────────────

    print("\n" + "=" * 65)
    print("INDIVIDUAL DATASET LISTING")
    print("=" * 65)
    for _, row in df_data.sort_values(["diet_label", "year"], ascending=[True, False]).iterrows():
        print(f"\n  [{row['diet_label']}]  {row['accession']}  {row['year']}  "
              f"{row['species']} / {row['organism']}")
        print(f"  Title:   {row['title'][:90]}")
        print(f"  Region:  {row['regions'] or 'not specified'}")
        print(f"  Sex: {row['sex']}  |  n={row['n_samples']}  |  RNA-seq: {row['assay_rnaseq']}")

    print("\nOutputs saved to:")
    print("  aim1_brain_diet_filtered.csv")


if __name__ == "__main__":
    run_survey()
