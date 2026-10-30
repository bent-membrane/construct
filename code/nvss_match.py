# -*- coding: utf-8 -*-
"""18 团：Okabe 晕椭率/θe × HEASARC NVSS Position_Angle"""
import ssl
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed

ctx = ssl._create_unverified_context()
NS = {"v": "http://www.ivoa.net/xml/VOTable/v1.3"}


def get(u, timeout=25):
    req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=timeout, context=ctx).read().decode("utf-8", "ignore")


CLUSTERS = [
    ("A68", -3.3, 0.552), ("A209", 16.3, 0.428), ("RXJ0142", -53.6, 0.799),
    ("A267", 26.5, 0.769), ("A291", 33.9, 0.771), ("A383", -5.3, 0.814),
    ("A521", -7.4, 0.176), ("A586", -87.7, 0.426), ("ZwCl0740", 19.3, 0.553),
    ("ZwCl0823", -47.3, 0.144), ("A611", 37.1, 0.517), ("A697", -52.1, 0.156),
    ("A1835", 45.5, 0.420), ("ZwCl1454", -16.1, 0.568), ("A2219", -78.2, 0.705),
    ("A2390", -45.0, 0.598), ("A2485", -8.6, 0.808), ("A2631", 78.6, 0.347),
]

FIXED = {
    "A68": (9.2758, 9.1592), "A383": (42.014, -3.529), "A291": (30.429, -2.197),
    "A1835": (210.259, 2.879), "A611": (120.245, 36.047), "A2219": (250.092, 46.706),
}


def sesame(name):
    u = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp/SNV?" + urllib.parse.quote(name)
    import re
    t = get(u, 15)
    m = re.search(r"<jradeg>([-\d.]+)</jradeg>", t)
    d = re.search(r"<jdedeg>([-\d.]+)</jdedeg>", t)
    return (float(m.group(1)), float(d.group(1))) if m and d else None


import base64
import struct


def parse_binary_votable(t):
    """HEASARC TAP 返回 BINARY 格式：char(*)=4字节长度+串，double=8字节大端。"""
    root = ET.fromstring(t)
    stream = root.find(".//v:BINARY/v:STREAM", NS)
    if stream is None or not stream.text:
        return []
    raw = base64.b64decode(stream.text.strip())
    pos, rows = 0, []
    while pos + 4 <= len(raw):
        n = struct.unpack(">I", raw[pos:pos + 4])[0]
        pos += 4
        if n > 200 or pos + n + 48 > len(raw):
            break
        name = raw[pos:pos + n].decode("utf-8", "ignore")
        pos += n
        ra, dec, flux, maj, mi, pa = struct.unpack(">6d", raw[pos:pos + 48])
        pos += 48
        rows.append({"name": name, "ra": ra, "dec": dec, "flux": flux,
                     "maj": maj, "mi": mi, "pa": pa})
    return rows


def nvss(ra, dec, radius=0.06):
    adql = (f"SELECT name,ra,dec,flux_20_cm,major_axis,minor_axis,position_angle "
            f"FROM nvss WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{radius}))")
    u = ("https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync?REQUEST=doQuery&LANG=ADQL"
         "&FORMAT=votable&QUERY=" + urllib.parse.quote(adql))
    return parse_binary_votable(get(u, 30))


def process(item):
    name, te, e = item
    pos = FIXED.get(name)
    try:
        if not pos:
            pos = sesame(name)
    except Exception as ex:
        return (name, te, e, None, [], f"sesame {type(ex).__name__}")
    if not pos:
        return (name, te, e, None, [], "no sesame")
    try:
        rows = nvss(*pos)
    except Exception as ex:
        return (name, te, e, pos, [], f"nvss {type(ex).__name__}")
    cand = [r for r in rows if r["pa"] is not None and r["maj"] > 20 and r["flux"] > 5]
    cand.sort(key=lambda r: -r["flux"])
    return (name, te, e, pos, cand, f"{len(rows)} src")


results = []
with ThreadPoolExecutor(max_workers=6) as pool:
    futs = [pool.submit(process, c) for c in CLUSTERS]
    for f in as_completed(futs):
        r = f.result()
        results.append(r)
        name, te, e, pos, cand, note = r
        if pos:
            print(f"{name:10s} RA={pos[0]:8.3f} DEC={pos[1]:8.3f}  θe={te:6.1f}  {note}  候选={len(cand)}", flush=True)
            for c in cand[:3]:
                print(f"      PA={c['pa']:6.1f}  Maj={c['maj']:6.1f}  Min={c['mi']:5.1f}  Flux={c['flux']:8.1f}", flush=True)
        else:
            print(f"{name:10s} {note}", flush=True)

print("\n=== 汇总（取每个团最亮的延伸源）===", flush=True)
pairs = []
for name, te, e, pos, cand, note in sorted(results):
    if cand:
        b = cand[0]
        d = abs(te - b["pa"]) % 180
        d = min(d, 180 - d)
        pairs.append((name, te, b["pa"], d, e))
        print(f"{name:10s} 晕θe={te:7.1f}  射电PA={b['pa']:6.1f}  差角={d:5.1f}  椭率={e:.2f}", flush=True)
print(f"\n匹配到 {len(pairs)} / 18 个团", flush=True)
if pairs:
    ds = [p[3] for p in pairs]
    import statistics
    print(f"差角：均值={statistics.mean(ds):.1f}  中位={statistics.median(ds):.1f}  "
          f"小于30°的比例={sum(1 for x in ds if x<30)}/{len(ds)}", flush=True)
