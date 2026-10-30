*Activity 01 · 2026-10-30*

# Inspect the Mol-1 repository

**A read-only survey of the Mol-1 structure-prediction codebase that mapped its layout, its Protenix lineage and its install paths, and found that CI never ran because it watched a branch that does not exist.**

## Aim

The request was a single line: "inspect this repository". The repository, a fork of `molence/Mol-1`, was new to the session, so the aim was an orientation a newcomer could act on. That meant covering what the code does, how it is organised, how it is installed and run, and anything that looked wrong.

Done meant a summary that answered those questions with file paths and numbers, plus a short list of findings worth acting on. Nothing was to be changed. Tests could not be run: PyTorch was not installed in the cloud container, and its network policy blocked most external hosts.

## Tools

| Tool | Version | Used for |
| --- | --- | --- |
| git | container default | History and branch list |
| find, wc, grep, sed | GNU coreutils | File inventory, line counts, reading headers and config |
| Python | 3.13.17 | Checked whether `torch` was importable (it was not) |

## Process

### 1. History, file tree and README

The first pass listed commits, the working tree and the README together. Seven commits make up the history: one "initial release" followed by README, contact-info and QR-code edits. The README describes Mol-1 as an all-atom predictor for proteins, nucleic acids, ligands and ions, driven by a `mol` command-line tool.

*First survey command, as run*

```bash
git log --oneline | head -20 && git status && find . -path ./.git -prune -o -type f -print | head -100 && cat README* | head -100
```

*Output: `d01253d Mol-1 initial release` at the root of history; 35 test files under `tests/`; example inputs, docs and a technical report PDF under `assets/`.*

### 2. Package layout and size

Line counts per subpackage showed where the weight of the code sits. Data handling is the largest part, larger than the model itself.

*Lines of Python per subpackage of `mol/`*

```bash
for d in mol/*/; do echo "$d $(find $d -name '*.py' | xargs cat 2>/dev/null | wc -l)"; done
```

![figure](figures/01-inspect-mol1-repository-fig1.svg)

*Figure 1 — Size of each `mol/` subpackage, in lines of Python (51,531 in the repository overall)*

> `data/` (featurisation, MSA and template handling) outweighs `model/`, which is typical for AlphaFold 3-style systems: turning a JSON job into tensors is more code than the network. `web_service/` is small but matters at run time, because it is the client for the public ColabFold MSA server that `mol pred` calls by default.

### 3. Lineage, configuration and CI

File headers showed the code is a modified fork of ByteDance's Protenix. The model file names the upstream file and lists the changes:

*Header of `mol/model/mol.py`, trimmed*

```python
# Based on Protenix: https://github.com/bytedance/Protenix
# Upstream file: protenix/model/protenix.py
# Changes: Mol-1 trunk: PairMixer blocks ahead of the Pairformer stack; triangle kernels
# selected per call, including FlashPairformer; ESM input removed.
```

Reading `pyproject.toml`, `requirements.txt` and the CI workflow next to `git branch -a` showed the main problem in the repository.

*The CI trigger, and the branches that actually exist*

```yaml
on:
  push:
    branches: [ "master" ]
  pull_request:
    branches: [ "master" ]
---
* ccr-312027ff-1vk8pa
  main
  remotes/origin/main
```

![figure](figures/01-inspect-mol1-repository-fig2.svg)

*Figure 2 — How a `mol pred` job flows through the code, and which steps need the network or extra binaries*

> The two blue boxes are the steps that later caused trouble in activity 03: the MSA step depends on a public server, and the model compiles its own CUDA LayerNorm on first use. At this stage they were only noted as requirements from the README; neither had been exercised.

## Results

The session produced a written orientation of the repository, delivered in chat. The findings:

- Mol-1 is a 576.8 M-parameter Protenix derivative (package `molence`, import name `mol`, Python 3.12–3.13), with five subcommands: `pred`, `json`, `msa`, `prep` and `doctor`.
- CI had never run: `.github/workflows/ci.yml` triggered on `master`, but the only default branch is `main`. This was fixed in activity 02.
- `requirements.txt` pins CUDA 13 wheels and a git-sourced `opt_core`; `pyproject.toml` uses ranges and omits `opt_core`. Only the latter is used by `pip install -e`.
- By default, prediction sends the sequence to the public ColabFold MMseqs2 server for MSA search.

### Limits and what was left

- Nothing was executed: no tests, no `mol doctor`. The findings come from reading code and configuration only.
- The survey did not catch two defects that surfaced later on real hardware: `ninja` is missing from the dependencies, and the LayerNorm build skips the T4's GPU architecture (activity 03).
