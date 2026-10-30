# -*- coding: utf-8 -*-
"""18 个团：晕椭率/位置角 × FIRST 射电 PA（并行版）"""
import ssl
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

ctx = ssl._create_unverified_context()


def get(url, timeout=12):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=timeout, context=ctx).read().decode("utf-8", "ignore")


CLUSTERS = [
    ("A68", -3.3, 0.552), ("A209", 16.3, 0.428), ("RXJ0142", -53.6, 0.799),
    ("A267", 26.5, 0.769), ("A291", 33.9, 0.771), ("A383", -5.3, 0.814),
    ("A521", -7.4, 0.176), ("A586", -87.7, 0.426), ("ZwCl0740", 19.3, 0.553),
    ("ZwCl0823", -47.3, 0.144), ("A611", 37.1, 0.517), ("A697", -52.1, 0.156),
    ("A1835", 45.5, 0.420), ("ZwCl1454", -16.1, 0.568), ("A2219", -78.2, 0.705),
    ("A2390", -45.0, 0.598), ("A2485", -8.6, 0.808), ("A2631", 78.6, 0.347),
]


def sesame(name):
    u = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp/SNV?" + urllib.parse.quote(name)
    t = get(u)
    import re
    m = re.search(r"<jradeg>([-\d.]+)</jradeg>", t)
    d = re.search(r"<jdedeg>([-\d.]+)</jdedeg>", t)
    return (float(m.group(1)), float(d.group(1))) if m and d else None


def first_cone(ra, dec, radius=3.0):
    q = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=VIII/92/first14"
         f"&-c={ra:.4f}+{dec:.4f}&-c.r={radius}"
         "&-out=RAJ2000,DEJ2000,Fpeak,Fint,Maj,Min,PA&-out.max=100")
    t = get(q, 20)
    rows, started = [], False
    for line in t.split("\n"):
        if line.startswith("RAJ2000"):
            started = True
            continue
        if not started or line.startswith("#") or not line.strip():
            continue
        p = [x.strip() for x in line.split("\t")]
        if len(p) < 7:
            continue
        try:
            rows.append({"ra": p[0], "dec": p[1], "fpeak": float(p[2] or 0),
                         "fint": float(p[3] or 0), "maj": float(p[4] or 0),
                         "min": float(p[5] or 0), "pa": float(p[6]) if p[6] else None})
        except Exception:
            pass
    return rows


def process(item):
    name, te, e = item
    try:
        pos = sesame(name)
        if not pos:
            return (name, te, e, None, [], "no sesame")
        rows = first_cone(*pos)
        cand = [r for r in rows if r["pa"] is not None and r["maj"] > 15 and r["fpeak"] > 5]
        cand.sort(key=lambda r: -r["fint"])
        return (name, te, e, pos, cand, f"{len(rows)} sources")
    except Exception as ex:
        return (name, te, e, None, [], f"ERR {type(ex).__name__}")


results = []
with ThreadPoolExecutor(max_workers=6) as pool:
    futs = {pool.submit(process, c): c for c in CLUSTERS}
    for f in as_completed(futs):
        r = f.result()
        results.append(r)
        name, te, e, pos, cand, note = r
        if pos:
            print(f"{name:10s} RA={pos[0]:8.3f} DEC={pos[1]:8.3f}  晕θe={te:6.1f}  "
                  f"{note}, 延伸候选={len(cand)}", flush=True)
            for c in cand[:3]:
                print(f"      PA={c['pa']:6.1f}  Maj={c['maj']:6.1f}  Min={c['min']:5.1f}  "
                      f"Fpeak={c['fpeak']:7.1f}", flush=True)
        else:
            print(f"{name:10s} {note}", flush=True)

print("\n=== 汇总 ===", flush=True)
pairs = []
for name, te, e, pos, cand, note in sorted(results):
    if cand:
        best = cand[0]
        d = abs(te - best["pa"]) % 180
        d = min(d, 180 - d)
        pairs.append((name, te, best["pa"], d, e))
        print(f"{name:10s} 晕θe={te:7.1f}  射电PA={best['pa']:6.1f}  差角={d:5.1f}  椭率={e:.2f}", flush=True)

print(f"\n匹配到 {len(pairs)} 个团", flush=True)
