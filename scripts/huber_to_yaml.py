#!/usr/bin/env python3
"""Download Huber reactor antineutrino flux dat files and convert to data/huber.yaml.

Data source: https://www1.phys.vt.edu/~pahuber/reactorfluxes/
Paper: Phys. Rev. C 84, 024617 (2011), arXiv:1106.0687
       P. Huber, "Determination of antineutrino spectra from nuclear reactors"

Dat file column format (14 columns):
  energy  flux  ...  ...  ...  ...  ...  ...  ...  ...  ...  ...  errLo  errHi
"""

import sys
import os
import urllib.request

BASE_URL = "https://www1.phys.vt.edu/~pahuber/reactorfluxes/"
ELEMENTS = ["U235", "Pu239", "Pu241"]
OUT_FILE = os.path.join(os.path.dirname(__file__), "..", "data", "huber.yaml")


def download_dat(element):
    fname = f"{element}-anti-neutrino-flux-250keV.dat"
    url = BASE_URL + fname
    print(f"Downloading {url} ...", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode()


def parse_dat(text):
    energies, fluxes, err_neg, err_pos = [], [], [], []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 14:
            continue
        energies.append(float(cols[0]))
        fluxes.append(float(cols[1]))
        err_neg.append(abs(float(cols[12])))
        err_pos.append(float(cols[13]))
    return energies, fluxes, err_neg, err_pos


def fmt_list(values):
    parts = []
    for v in values:
        if abs(v) < 1e-3 or abs(v) >= 1e4:
            parts.append(f"{v:.4g}")
        else:
            s = f"{v:.6g}"
            parts.append(s)
    return "[" + ", ".join(parts) + "]"


def main():
    os.makedirs(os.path.dirname(os.path.abspath(OUT_FILE)), exist_ok=True)

    entries = []
    for element in ELEMENTS:
        text = download_dat(element)
        energies, fluxes, err_neg, err_pos = parse_dat(text)
        entries.append((element, energies, fluxes, err_neg, err_pos))
        print(f"  {element}: {len(energies)} points, "
              f"E = {energies[0]:.3f} - {energies[-1]:.3f} MeV", file=sys.stderr)

    lines = []
    lines.append('source: "https://journals.aps.org/prc/abstract/10.1103/PhysRevC.84.024617"')
    lines.append('arxiv: "1106.0687"')
    lines.append('reference: "Phys. Rev. C 84, 024617 (2011)"')
    lines.append('authors: "P. Huber"')
    lines.append("note: |")
    lines.append("  Reactor antineutrino flux from conversion of ILL beta spectra.")
    lines.append("  energy: bin center [MeV], bin width = 250 keV")
    lines.append("  flux: N_antineutrino per fission per MeV")
    lines.append("  errNeg: lower 1-sigma absolute uncertainty")
    lines.append("  errPos: upper 1-sigma absolute uncertainty")
    lines.append("  Covers U235, Pu239, Pu241 (U238 is not in this dataset; use mueller.yaml).")
    lines.append("")
    lines.append("data:")
    for element, energies, fluxes, err_neg, err_pos in entries:
        lines.append(f"  - name: {element}")
        lines.append(f"    energy: {fmt_list(energies)}")
        lines.append(f"    flux:   {fmt_list(fluxes)}")
        lines.append(f"    errNeg: {fmt_list(err_neg)}")
        lines.append(f"    errPos: {fmt_list(err_pos)}")
        lines.append("")

    with open(OUT_FILE, "w") as f:
        f.write("\n".join(lines))

    print(f"Written: {OUT_FILE}", file=sys.stderr)


if __name__ == "__main__":
    main()
