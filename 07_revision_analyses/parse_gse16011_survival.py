#!/usr/bin/env python3
"""Parse Gravendeel 2009 Supplementary Table 1 (GSE16011) into a clean survival table
and match the per-patient 'Database number' to the GEO series-matrix sample titles
('glioma N' / 'control N'). Survival time is in years; status Dead=event1, Alive/Lost=censored.

Inputs:
  data/raw/GEO/00085472can092307-sup-stabs_1-6.pdf   (supplementary tables)
  data/raw/GEO/GSE16011_series_matrix.txt.gz          (for !Sample_title -> number, !Sample_geo_accession)
Outputs:
  data/processed/gse16011_survival_from_supp.csv
"""
import os, re, gzip, csv, sys

PROJ = os.getcwd()
PDF  = os.path.join(PROJ, "data/raw/GEO/00085472can092307-sup-stabs_1-6.pdf")
SM   = os.path.join(PROJ, "data/raw/GEO/GSE16011_series_matrix.txt.gz")
OUT  = os.path.join(PROJ, "data/processed/gse16011_survival_from_supp.csv")

# ---- 1. Extract ONLY Supplementary Table 1 (patient characteristics) ----
# The PDF holds 6 supplementary tables; later ones (LOH / mutation) are also
# integer-keyed and would clobber Table S1. Gate strictly to Table S1: start when
# its header (Survival / Alive / histological) appears, stop at the next table's
# header (LOH / mutation / etc.). Read status & survival-years by column position.
import pdfplumber

def to_num(s):
    s = (s or "").replace(",", ".").strip()
    m = re.search(r"\d+(\.\d+)?", s)
    return float(m.group()) if m else None

STOPKW = ("loh", "mutation", "copy number", "methylation", "idh1", "1p/19q",
          "19q", "probe", "cgh", "primer", "sequence")
records = {}
t1_started = False
stop = False
with pdfplumber.open(PDF) as pdf:
    for page in pdf.pages:
        if stop:
            break
        tbls = page.extract_tables({"vertical_strategy": "text",
                                    "horizontal_strategy": "text"}) or []
        for t in tbls:
            if not t:
                continue
            hdr = " ".join(" ".join((c or "") for c in row) for row in t[:2]).lower()
            if "survival" in hdr or ("histological" in hdr and "alive" in hdr):
                t1_started = True
            elif t1_started and any(k in hdr for k in STOPKW) and "survival" not in hdr:
                stop = True
                break
            for r in t:
                if not r:
                    continue
                dbtok = (r[0] or "").strip()
                if not re.fullmatch(r"\d+", dbtok):
                    continue
                dbn = int(dbtok)
                status_cell = (r[8] if len(r) > 8 else "") or ""
                years_cell  = (r[9] if len(r) > 9 else "") or ""
                hist = (r[2] if len(r) > 2 else "") or ""
                age  = to_num(r[3]) if len(r) > 3 else None
                rowtxt = " ".join((c or "") for c in r).lower()
                sl = status_cell.strip().lower()
                # "Lost to follow-up" may arrive cid-garbled (e.g. 'follow(cid:882)up');
                # detect via the whole-row text as well. Lost/Alive = censored (event 0).
                is_lost = ("lost" in sl) or ("follow" in sl) or ("lost" in rowtxt and "follow" in rowtxt)
                if "dead" in sl:
                    ev, stl = 1, "Dead"
                elif "alive" in sl:
                    ev, stl = 0, "Alive"
                elif is_lost:
                    ev, stl = 0, "Lost"
                elif "control" in rowtxt:
                    ev, stl = None, "control"
                else:
                    ev, stl = None, status_cell.strip()
                # survival time strictly from the 'Survival (years)' column (no fallback
                # scan — that would grab Age/KPS and inflate values).
                yr = to_num(years_cell)
                if yr is not None and yr > 60:
                    yr = None
                records[dbn] = dict(db=dbn, status=stl, event=ev, years=yr,
                                    hist=hist, age=age)

# ---- 2. Parse GEO series matrix: title number -> GSM accession ----
titles, accs = [], []
with gzip.open(SM, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            accs = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if titles and accs:
            break

geo = []  # (number, kind, gsm, title)
for i, t in enumerate(titles):
    m = re.search(r"(control|glioma)\s+(\d+)", t, re.I)
    if m:
        geo.append((int(m.group(2)), m.group(1).lower(), accs[i] if i < len(accs) else "", t))

# ---- 3. Join GEO arrays to supplementary survival by Database number ----
out_rows = []
matched = 0
for num, kind, gsm, title in geo:
    rec = records.get(num)
    if rec is None:
        out_rows.append(dict(geo_title=title, gsm=gsm, db_number=num, kind=kind,
                             status="", OS_event="", OS_years="", OS_days="",
                             histology="", note="no_supp_record"))
        continue
    ev = rec.get("event")
    yr = rec.get("years")
    has_surv = (kind == "glioma" and ev is not None and yr is not None)
    if has_surv:
        matched += 1
    out_rows.append(dict(
        geo_title=title, gsm=gsm, db_number=num, kind=kind,
        status=rec.get("status", ""),
        OS_event=(ev if ev is not None else ""),
        OS_years=(yr if yr is not None else ""),
        OS_days=(round(yr * 365.25) if yr is not None else ""),
        histology=rec.get("hist", ""),
        note=("ok" if has_surv else "control_or_missing"),
    ))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as o:
    w = csv.DictWriter(o, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)

# ---- 4. Compact summary (small output to avoid channel issues) ----
n_geo = len(geo)
n_glioma = sum(1 for g in geo if g[1] == "glioma")
n_event1 = sum(1 for r in out_rows if r["OS_event"] == 1)
n_event0 = sum(1 for r in out_rows if r["OS_event"] == 0)
print(f"GEO_samples={n_geo} glioma={n_glioma} controls={n_geo-n_glioma}")
print(f"supp_records={len(records)}")
print(f"matched_with_survival={matched} (events={n_event1} censored={n_event0})")
print(f"written={OUT}")
