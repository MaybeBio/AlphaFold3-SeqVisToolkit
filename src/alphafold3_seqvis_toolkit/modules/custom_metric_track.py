import glob
import json
import os
import csv
from typing import List, Dict, Any, Optional

def _safe_load_json(p: str) -> Dict[str, Any]:
    try:
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        return {"__file__": p, "__error__": str(e)}

def track_metrics(folder: str, pattern: str = "*.json", out_csv: Optional[str] = None) -> List[Dict[str, Any]]:
    paths = sorted(glob.glob(os.path.join(folder, pattern)))
    records = []
    for p in paths:
        data = _safe_load_json(p)
        data["__file__"] = os.path.basename(p)
        records.append(data)
    if out_csv and records:
        # flatten first-level keys only
        keys = sorted({k for r in records for k in r.keys()})
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in records:
                w.writerow(r)
    return records