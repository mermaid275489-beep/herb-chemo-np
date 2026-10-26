"""Download every herb in legacy TCMSP 2.3 (old.tcmsp-e.com): herb list, ingredients
(OB, DL and other ADME columns) and ingredient -> target-name edges.

Resumable: each herb is cached as raw/tcmsp/herbs/<n>.json and skipped on rerun.
Polite: one request at a time with a pause between requests (robots.txt allows all).

  python build/fetch_tcmsp.py [--limit N] [--pause 1.0] [--workers 3]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

BASE = "https://old.tcmsp-e.com"
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw" / "tcmsp"
UA = "Mozilla/5.0 (compatible; academic network pharmacology research/1.0)"


def data_arrays(html: str) -> list[list[dict]]:
    out = []
    for m in re.finditer(r"data\s*:\s*(\[\{.*?\}\])\s*,?\s*\n\s*pageSize", html, re.S):
        try:
            out.append(json.loads(m.group(1)))
        except json.JSONDecodeError:
            continue
    return out


def any_arrays(html: str) -> list[list[dict]]:
    out = []
    for m in re.finditer(r"data\s*:\s*(\[\{.*?\}\])", html, re.S):
        try:
            out.append(json.loads(m.group(1)))
        except json.JSONDecodeError:
            continue
    return out


class Client:
    def __init__(self, pause: float):
        self.s = requests.Session()
        self.s.headers.update({"User-Agent": UA})
        self.pause = pause
        self.token = None

    def refresh(self):
        home = self.s.get(f"{BASE}/tcmsp.php", timeout=60)
        home.raise_for_status()
        m = re.search(r'name=["\']token["\']\s+value=["\']([^"\']+)', home.text)
        if not m:
            raise RuntimeError("TCMSP token not found")
        self.token = m.group(1)

    def get(self, path: str, **params) -> str:
        for attempt in range(4):
            if self.token is None:
                self.refresh()
            try:
                r = self.s.get(f"{BASE}/{path}", params={**params, "token": self.token}, timeout=90)
                r.raise_for_status()
                time.sleep(self.pause)
                return r.text
            except requests.RequestException as e:
                print(f"  retry {attempt + 1}: {e}", flush=True)
                self.token = None
                time.sleep(5 * (attempt + 1))
        raise RuntimeError(f"failed: {path} {params}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--pause", type=float, default=1.0)
    ap.add_argument("--workers", type=int, default=1, help="parallel connections (each has its own session)")
    ap.add_argument("--reverse", action="store_true", help="work from the end of the list")
    a = ap.parse_args(argv)
    (RAW / "herbs").mkdir(parents=True, exist_ok=True)
    c = Client(a.pause)

    herbs_path = RAW / "herbs.json"
    if not herbs_path.exists():
        html = c.get("browse.php", qc="herbs")
        herbs = next(arr for arr in any_arrays(html) if arr and "herb_cn_name" in arr[0])
        herbs_path.write_text(json.dumps(herbs, ensure_ascii=False, indent=0), encoding="utf-8")
    herbs = json.loads(herbs_path.read_text(encoding="utf-8"))
    print(f"[tcmsp] {len(herbs)} herbs in browse list", flush=True)

    todo = [(i, h) for i, h in list(enumerate(herbs))[: a.limit or None]
            if not (RAW / "herbs" / f"{i:03d}.json").exists()]

    if a.reverse:
        todo.reverse()

    def fetch(item, client):
        i, h = item
        if (RAW / "herbs" / f"{i:03d}.json").exists():   # another run may have fetched it meanwhile
            return
        html = client.get("tcmspsearch.php", qr=h["herb_en_name"], qsr="herb_en_name")
        ingredients, edges = [], []
        for arr in data_arrays(html):
            keys = set(arr[0])
            if {"MOL_ID", "molecule_name", "ob", "dl"} <= keys:
                ingredients = arr
            elif {"MOL_ID", "molecule_name", "target_name"} <= keys:
                edges = arr
        (RAW / "herbs" / f"{i:03d}.json").write_text(
            json.dumps({"herb": h, "ingredients": ingredients, "edges": edges}, ensure_ascii=False), encoding="utf-8")
        print(f"  {i:03d} {h['herb_cn_name']} {h['herb_en_name']}: {len(ingredients)} ingredients, "
              f"{len(edges)} target edges", flush=True)

    # The server throttles each connection (~25 KB/s, no compression), so a few parallel
    # connections, each with its own session, raise throughput without hammering it.
    clients = [c] + [Client(a.pause) for _ in range(a.workers - 1)]
    with ThreadPoolExecutor(a.workers) as ex:
        list(ex.map(lambda k: [fetch(item, clients[k]) for item in todo[k::a.workers]], range(a.workers)))
    done = len(list((RAW / "herbs").glob("*.json")))
    print(f"[tcmsp] cached {done}/{len(herbs)} herbs", flush=True)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
