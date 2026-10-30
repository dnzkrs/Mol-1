"""Compare Mol-1 predictions with the experimental structures of a UniProt entry.

Usage (Colab or anywhere with internet):
    python scripts/compare_crystal.py PRED_DIR [UNIPROT]            # all PDB entries of UNIPROT
    python scripts/compare_crystal.py PRED_DIR --local a.pdb b.cif ...  # local files
"""
import glob
import os
import re
import sys

import gemmi
import numpy as np
import requests
from Bio.Align import PairwiseAligner

CACHE = os.path.join(os.getcwd(), "pdb_cache")


# ---------- fetching ----------
def pdb_entries(uniprot):
    q = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": uniprot}},
        "return_type": "entry", "request_options": {"return_all_hits": True}}
    r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=q, timeout=60)
    r.raise_for_status()
    ids = sorted(h["identifier"] for h in r.json().get("result_set", []))
    out = []
    for pid in ids:
        info = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}", timeout=60).json()
        res = (info.get("rcsb_entry_info", {}).get("resolution_combined") or [None])[0]
        method = info.get("exptl", [{}])[0].get("method", "?")
        year = (info.get("rcsb_accession_info", {}).get("initial_release_date") or "?")[:4]
        title = info.get("struct", {}).get("title", "")
        out.append(dict(id=pid, resolution=res, method=method, year=year, title=title))
    return out


def fetch_cif(pid):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f"{pid}.cif")
    if not os.path.exists(path):
        r = requests.get(f"https://files.rcsb.org/download/{pid}.cif", timeout=120)
        r.raise_for_status()
        open(path, "w").write(r.text)
    return path


# ---------- structure handling ----------
def protein_chains(path):
    """{chain: [(one_letter, {atom_name: xyz}), ...]} for amino-acid residues, first model."""
    st = gemmi.read_structure(path)
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    chains = {}
    for ch in st[0]:
        residues = []
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid() or res.find_atom("CA", "*") is None:
                continue
            # gemmi maps modified residues to their parent in lower case (MSE -> m)
            aa = (info.one_letter_code.strip() or "X").upper()
            residues.append((aa, {a.name: np.array(a.pos.tolist()) for a in res}))
        if residues:
            chains[ch.name] = residues
    return chains


def align_map(seq_a, seq_b):
    """Index pairs (i, j) of aligned residues between two sequences."""
    al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                         open_gap_score=-5, extend_gap_score=-0.5,
                         end_gap_score=0)
    a = al.align(seq_a, seq_b)[0]
    pairs = []
    for (a0, a1), (b0, b1) in zip(*a.aligned):
        pairs += list(zip(range(a0, a1), range(b0, b1)))
    return pairs


def kabsch(P, Q):
    """Rotation R and translation t that map P onto Q (P @ R + t ~ Q)."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return R, qc - pc @ R


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def tm_score(P, Q, L):
    """TM-score of the matched pairs, optimised over superpositions seeded from
    fragments (as in TM-score), normalised by L."""
    d0 = max(1.24 * np.cbrt(max(L, 19) - 15) - 1.8, 0.5)
    n = len(P)
    best = 0.0
    for frag in (n, n // 2, n // 4, 20, 12):
        frag = max(min(frag, n), 4)
        for start in range(0, n - frag + 1, max(frag // 2, 1)):
            idx = np.arange(start, start + frag)
            for _ in range(20):
                R, t = kabsch(P[idx], Q[idx])
                d = np.linalg.norm(P @ R + t - Q, axis=1)
                best = max(best, float((1 / (1 + (d / d0) ** 2)).sum() / L))
                new = np.where(d < d0 + 1.0)[0]
                if len(new) < 4 or np.array_equal(new, idx):
                    break
                idx = new
    return best


def gdt_ts(P, Q):
    """GDT_TS with superpositions seeded from fragments (close to LGA)."""
    n = len(P)
    best = {c: 0 for c in (1, 2, 4, 8)}
    for frag in (n, 40, 20, 10):
        frag = min(frag, n)
        for start in range(0, n - frag + 1, max(frag // 2, 1)):
            idx = np.arange(start, start + frag)
            for c in (1, 2, 4, 8):
                sub = idx
                for _ in range(10):
                    R, t = kabsch(P[sub], Q[sub])
                    d = np.linalg.norm(P @ R + t - Q, axis=1)
                    best[c] = max(best[c], int((d < c).sum()))
                    new = np.where(d < c)[0]
                    if len(new) < 3 or np.array_equal(new, sub):
                        break
                    sub = new
    return 100 * np.mean([best[c] / n for c in (1, 2, 4, 8)])


def lddt(pred_atoms, ref_atoms, inclusion=15.0):
    """Superposition-free lDDT over the given atoms (any atom set, e.g. CA)."""
    dr = np.linalg.norm(ref_atoms[:, None] - ref_atoms[None], axis=-1)
    dp = np.linalg.norm(pred_atoms[:, None] - pred_atoms[None], axis=-1)
    mask = (dr < inclusion) & ~np.eye(len(dr), dtype=bool)
    diff = np.abs(dr - dp)[mask]
    per_pair = np.mean([(diff < t) for t in (0.5, 1, 2, 4)], axis=0)
    per_res = np.zeros(len(dr))
    rows = np.where(mask)[0]
    np.add.at(per_res, rows, per_pair)
    counts = mask.sum(1)
    return float(per_pair.mean()), np.divide(per_res, counts, out=np.zeros_like(per_res), where=counts > 0)


# ---------- comparison ----------
def compare(model_path, crystal_path, verbose=False):
    m_res = next(iter(protein_chains(model_path).values()))
    m_seq = "".join(a for a, _ in m_res)
    plddt = np.array([float(line.split()[12]) for line in open(model_path)
                      if line.startswith("ATOM") and line.split()[2] == "CA"])
    results = []
    for chain, c_res in protein_chains(crystal_path).items():
        c_seq = "".join(a for a, _ in c_res)
        pairs = align_map(m_seq, c_seq)
        ident = sum(m_seq[i] == c_seq[j] for i, j in pairs)
        if len(pairs) < 50 or ident / len(pairs) < 0.9:
            continue
        mi = np.array([i for i, _ in pairs])
        ci = np.array([j for _, j in pairs])
        P = np.array([m_res[i][1]["CA"] for i in mi])
        Q = np.array([c_res[j][1]["CA"] for j in ci])
        R, t = kabsch(P, Q)
        dev = np.linalg.norm(P @ R + t - Q, axis=1)
        # heavy atoms present in both, same superposition
        pa, qa = [], []
        for i, j in pairs:
            for name, xyz in m_res[i][1].items():
                if name in c_res[j][1] and m_seq[i] == c_seq[j]:
                    pa.append(xyz)
                    qa.append(c_res[j][1][name])
        pa, qa = np.array(pa), np.array(qa)
        l_ca, l_ca_res = lddt(P, Q)
        r = dict(chain=chain, n_crystal=len(c_seq), n_matched=len(pairs),
                 identity=ident / len(pairs),
                 ca_rmsd=rmsd(P @ R + t, Q), allatom_rmsd=rmsd(pa @ R + t, qa),
                 tm=tm_score(P, Q, len(c_seq)), gdt_ts=gdt_ts(P, Q),
                 lddt_ca=l_ca, lddt_all=lddt(pa, qa)[0],
                 frac_ca_lt1=float((dev < 1).mean()), frac_ca_lt2=float((dev < 2).mean()),
                 mi=mi, dev=dev, lddt_res=l_ca_res, plddt=plddt[mi], m_seq=m_seq,
                 m_res=m_res, c_res=c_res, pairs=pairs, R=R, t=t)
        results.append(r)
    return results


def active_site(r, motif=r"HC..G..R"):
    """CA and all-atom RMSD of the PTP-like P-loop (HCxxGxxR...) under the global superposition."""
    m = re.search(motif, r["m_seq"])
    if not m:
        return None
    span = range(m.start(), m.end() + 2)  # include the following T/S
    pmap = dict(r["pairs"])
    pa, qa, ca_p, ca_q = [], [], [], []
    for i in span:
        if i not in pmap:
            continue
        mres, cres = r["m_res"][i][1], r["c_res"][pmap[i]][1]
        for name, xyz in mres.items():
            if name in cres:
                pa.append(xyz)
                qa.append(cres[name])
        ca_p.append(mres["CA"])
        ca_q.append(cres["CA"])
    R, t = r["R"], r["t"]
    return dict(motif=r["m_seq"][span.start:span.stop], start=span.start + 1,
                ca=rmsd(np.array(ca_p) @ R + t, np.array(ca_q)),
                allatom=rmsd(np.array(pa) @ R + t, np.array(qa)))


def segments(idx):
    out = []
    for i in idx:
        if out and i == out[-1][1] + 1:
            out[-1][1] = i
        else:
            out.append([i, i])
    return ", ".join(f"{a}-{b}" if a != b else f"{a}" for a, b in out)


def report(pred_dir, crystals):
    models = sorted(glob.glob(os.path.join(pred_dir, "*_sample_*.cif")))
    top = [m for m in models if m.endswith("_sample_0.cif")][0]
    print("\n=== Top-ranked model vs every chain of every entry ===")
    print(f"{'entry':<8}{'res(A)':>7}{'chain':>6}{'match':>7}{'CA-RMSD':>9}{'allatom':>9}"
          f"{'TM':>7}{'GDT_TS':>8}{'lDDT-CA':>9}{'lDDT':>7}")
    best = None
    for c in crystals:
        for r in compare(top, c["path"]):
            print(f"{c['id']:<8}{(c['resolution'] or 0):>7.2f}{r['chain']:>6}{r['n_matched']:>7}"
                  f"{r['ca_rmsd']:>9.2f}{r['allatom_rmsd']:>9.2f}{r['tm']:>7.3f}{r['gdt_ts']:>8.1f}"
                  f"{r['lddt_ca']:>9.3f}{r['lddt_all']:>7.3f}")
            key = (c["resolution"] or 99, -r["n_matched"])
            if best is None or key < best[0]:
                best = (key, c, r)
    if best is None:
        print("No chain matched the model sequence.")
        return
    _, c, r = best
    print(f"\n=== Detail: best-resolution reference {c['id']} chain {r['chain']} "
          f"({c['resolution']} A, {c['method']}, released {c['year']}) ===")
    print(c.get("title", ""))
    print(f"Residues compared: {r['n_matched']} of {len(r['m_seq'])} modelled, "
          f"{r['n_crystal']} observed in crystal; identity {100 * r['identity']:.1f}%")
    print(f"CA within 1 A: {100 * r['frac_ca_lt1']:.1f}%   within 2 A: {100 * r['frac_ca_lt2']:.1f}%")
    resnum = r["mi"] + 1
    print("Model residues with CA deviation > 2 A:", segments(resnum[r["dev"] > 2]) or "none")
    print("Model residues not observed in crystal:",
          segments(sorted(set(range(1, len(r['m_seq']) + 1)) - set(resnum))) or "none")
    cc = np.corrcoef(r["plddt"], r["lddt_res"])[0, 1] if r["lddt_res"].std() > 0 else float("nan")
    print(f"Correlation of pLDDT with per-residue lDDT-CA: r = {cc:.2f}")
    print(f"Mean pLDDT/100 = {r['plddt'].mean() / 100:.3f}  vs  actual lDDT-CA = {r['lddt_ca']:.3f}")
    a = active_site(r)
    if a:
        print(f"P-loop {a['motif']} (model residues {a['start']}-{a['start'] + len(a['motif']) - 1}): "
              f"CA RMSD {a['ca']:.2f} A, all-atom RMSD {a['allatom']:.2f} A")
    print("\n=== All five samples vs that chain ===")
    for m in models:
        rr = [x for x in compare(m, c["path"]) if x["chain"] == r["chain"]][0]
        print(f"{os.path.basename(m):<28} CA-RMSD {rr['ca_rmsd']:.2f}  TM {rr['tm']:.3f}  "
              f"lDDT-CA {rr['lddt_ca']:.3f}  GDT_TS {rr['gdt_ts']:.1f}")


if __name__ == "__main__":
    pred_dir = sys.argv[1]
    if "--local" in sys.argv:
        crystals = []
        for path in sys.argv[sys.argv.index("--local") + 1:]:
            st = gemmi.read_structure(path)
            crystals.append(dict(id=os.path.splitext(os.path.basename(path))[0][:8],
                                 resolution=st.resolution or None, method="local", year="?",
                                 title=st.info["_struct.title"] if "_struct.title" in st.info else "", path=path))
    else:
        uniprot = sys.argv[2] if len(sys.argv) > 2 else "Q7WUJ1"
        crystals = pdb_entries(uniprot)
        print(f"{len(crystals)} PDB entries for {uniprot}:")
        for c in crystals:
            print(f"  {c['id']}  {c['resolution']} A  {c['method']}  {c['year']}  {c['title'][:70]}")
            c["path"] = fetch_cif(c["id"])
    report(pred_dir, crystals)
