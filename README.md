# EDA Protein Consequences: Transcript, RNA Evidence and Haplotype Phasing

**Reproducible Python workflows for EDA-A1/A2 coding consequences, calibrated RNA evidence, longitudinal inference and probabilistic phase reconstruction.**

[View repository](https://github.com/Mugambi24/EDA_Protein_Consequences) · [Task specification](instruction.md) · [Reference implementation](solution/) · [Tests](tests/)

## Overview

This repository contains a scientific-computing task built around the human **EDA (ectodysplasin-A)** gene and its A1/A2 transcript products. It combines transcript-aware sequence reconstruction with several forms of RNA evidence and uncertainty-aware inference. The implementation is designed to test reproducible reasoning across linked biological and computational constraints—not simply classify individual variants.

**Data and clinical scope:** The participant genotypes, molecule counts, mixtures, assay measurements, case labels and planted truth values used in this task are **synthetic**. Real EDA gene information, transcript accessions and reference sequences provide biological context. These simulated cases are **not patients or clinical observations**, and the workflow is **not a validated diagnostic system**.

The task defines **10 inputs and 8 JSON outputs**. The complete command-line interface, file contracts, coordinate conventions, expected behaviour and resource limits are specified in [`instruction.md`](instruction.md), which is authoritative if a short overview here omits a detail.

## Scientific motivation

The EDA gene encodes ectodysplasin-A, a signalling protein involved in the development of ectoderm-derived structures. Pathogenic EDA variants can cause **X-linked hypohidrotic ectodermal dysplasia (XLHED)**. Different EDA transcripts and coding products create an instructive setting for modelling how the same phased genomic edits can affect more than one coding sequence.

Real molecular interpretation also depends on evidence quality, isoform mixtures, assay-specific calibration and the uncertainty of inferred allele or haplotype states. This project brings those computational challenges together in a controlled synthetic benchmark.

## Objectives

1. Reconstruct transcript-specific coding sequences from phased genomic edits and determine their protein consequences.
2. Evaluate quality-filtered RNA junction evidence and estimate isoform contributions using calibration standards.
3. Reconcile variant records against versioned reference and transcript information.
4. Model batch-specific, phased and longitudinal synthetic RNA measurements while representing uncertainty.
5. Reconstruct long-range haplotype phase from noisy and sometimes contradictory molecule-link observations.
6. Check reproducibility, output contracts and scientific consistency through independent implementations and automated tests.

## Computational workflow

| Component | Main computational challenge | Relevant code or data |
| --- | --- | --- |
| Coding consequences | Apply compatible, cis-phased genomic-plus edits to EDA transcript models; reconstruct and translate coding products | `solution/coding.py`, `environment/data/` |
| RNA junction evidence | Filter reads by quality, resolve UMI evidence and distinguish unique from shared support | `solution/evidence.py` |
| Isoform inference | Use spike-in controls and likelihood-based estimates with uncertainty intervals | `solution/inference.py` |
| Variant crosswalk | Check accession/version, coordinates, reference alleles and transcript consequence consistency | `environment/data/public_eda_variants.json`, `authoring/provenance/` |
| Batch and phased RNA | Account for assay calibration, nonlinear detection effects and informative mixture composition | `environment/data/*assay.json`, protocol files |
| Longitudinal RNA | Separate donor, time, treatment, batch capture and cross-mapping effects in matched synthetic panels | `environment/data/longitudinal_assay.json` |
| Phase graph | Infer multi-site phase and posterior uncertainty from calibrated molecule links | `environment/data/phase_graph_assay.json` |
| Output auditing | Check expected structural and scientific properties of generated results | `solution/audit.py`, `tests/` |

These descriptions summarise the supplied task documentation; consult the specification for exact requirements and output filenames.

### 1. Transcript-aware coding consequences

The workflow applies genomic-plus allele edits using **zero-based start positions**; insertions occur immediately before the indicated base. For the relevant transcript representation, exons are reverse-complemented and assembled in transcript order before translation.

Important edge cases include duplicate case identifiers, duplicate edit starts, overlapping variants, incompatible alleles and edits affecting initiation or termination. The task distinguishes a **lost native stop codon** from an intact termination codon merely shifted by an upstream in-frame insertion. Exact coding/protein outputs and associated checks are defined in the task specification.

### 2. Quality-controlled RNA evidence and isoform estimation

RNA junction molecules are assessed using mapping/base-quality criteria and reproducible UMI-level selection rules. Unique versus shared junction support is handled explicitly. Fixed spike-ins inform two-component likelihood inference, including estimates and profile intervals.

Validation addresses not just agreement with expected values but also internal likelihood identities, numerical precision and appropriate probability or interval bounds.

### 3. Calibrated batch, phased and longitudinal assays

The synthetic assay extensions introduce batch-dependent measurement behaviour, imperfect detection, potentially contaminated aliquots and calibration uncertainty. The phased deconvolution task requires **at least four distinct informative mixed compositions**, including supported three-product mixtures.

Longitudinal inference uses **four matched panels per synthetic donor** and seeks to distinguish treatment-related changes from time trends, capture effects and cross-mapping. Protocol-specific rules are documented in the four files below:

- [`batch_assay_protocol.txt`](environment/data/batch_assay_protocol.txt)
- [`phased_rna_protocol.txt`](environment/data/phased_rna_protocol.txt)
- [`longitudinal_protocol.txt`](environment/data/longitudinal_protocol.txt)
- [`phase_graph_protocol.txt`](environment/data/phase_graph_protocol.txt)

### 4. Long-range phase graph inference

The phase-graph extension evaluates joint phase inference across **10–12 synthetic heterozygous sites per case**. Links may be contradictory across cycles, and measurement error varies by batch and link class. In one synthetic batch, the apparent direction of cis/trans evidence is inverted, making naive vote-counting unreliable.

The first site anchors phase orientation. Remaining sites are inferred from the combined evidence, with **marginal posterior probabilities**, a **whole-configuration maximum-a-posteriori (MAP) probability**, and **null/no-call states** where the evidence is insufficient. Global uniqueness of molecule-link replicate identifiers, sensitivity to calibration and input-order invariance are among the relevant validation concerns.

## Implementation

The constrained reference implementation uses **Python standard-library modules** and is organised into five files:

```text
solution/
├── cli.py        # Command-line interface and input/output coordination
├── coding.py     # Transcript reconstruction and coding consequences
├── evidence.py   # Evidence processing
├── inference.py  # Statistical inference
├── audit.py      # Result auditing
└── solve.sh      # Task solution setup/entry point
```

The project also includes an independently developed implementation in [`authoring/provenance/independent_submission/`](authoring/provenance/independent_submission/) for cross-checking selected computational results. Development and validation utilities may have additional requirements; these should not be confused with the standard-library restriction on the constrained solution.

## Repository layout

```text
EDA_Protein_Consequences/
├── README.md
├── instruction.md                    # Authoritative task and CLI specification
├── task.toml                         # Task metadata/configuration
├── solution/                         # Reference Python implementation
│   ├── cli.py
│   ├── coding.py
│   ├── evidence.py
│   ├── inference.py
│   ├── audit.py
│   └── solve.sh
├── environment/
│   ├── Dockerfile
│   ├── data/                         # Reference, practice and assay inputs
│   └── workbench/                    # Supporting workbench code
├── tests/
│   ├── Dockerfile
│   ├── test.sh
│   ├── test_workflow.py
│   ├── run_isolated.py
│   └── ...                            # Test helpers and validation fixtures
└── authoring/
    ├── provenance/                   # Generators, audits and independent checks
    └── evidence/                     # Recorded validation evidence
```

Some verifier-side fixtures in `tests/` contain expected or planted synthetic truth values. Their inclusion in any public distribution should be reviewed against the task provider's rules before release.

## Getting started

### Prerequisites

- Git
- Python 3 with virtual-environment support
- `pytest` for local test execution
- Docker **if** using the containerised task/verifier environment

The solution's core Python modules are designed around the standard library; test tools and development checks can require separate packages. Refer to [`authoring/provenance/requirements.txt`](authoring/provenance/requirements.txt) and the Dockerfiles when running those components. A specific Python version or dependency set should be confirmed from the provided environment before attempting to reproduce a full verifier run.

### Clone and prepare an environment

```bash
git clone https://github.com/Mugambi24/EDA_Protein_Consequences.git
cd EDA_Protein_Consequences
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip pytest
```

### Inspect the exact task interface

Read [`instruction.md`](instruction.md) before invoking the reference solution. It defines the required input arguments, the output JSON contracts, paths, validation rules and resource limits. Do **not** assume that calling `cli.py` without its required arguments will run the complete workflow.

### Run local automated tests

From the repository root, with the required environment active:

```bash
python -m pytest tests/test_workflow.py -v
```

For test discovery without running assertions:

```bash
python -m pytest tests/test_workflow.py --collect-only -q
```

The repository also includes [`tests/test.sh`](tests/test.sh), a task-specific verifier entry point intended for its configured environment. Some tests may depend on container paths, fixture permissions, a prepared solution artifact or runner-specific settings. A local pytest invocation **is not equivalent** to a successful isolated/container verification unless those conditions match.

**Results are not claimed in this README without a reproducible execution record.** Run the relevant commands in the target environment and record the actual pass/fail totals and tool versions when reporting validation.

## Testing and validation approach

The documentation describes a **nine-group, 72-test verifier release**. Verification concerns include:

- Transcript sequence, coding consequences, protein audits and output schemas.
- Quality filtering, UMI handling, calibrated likelihood estimates and numerical precision.
- Batch-specific and phased mixture identifiability.
- Longitudinal estimation with calibrated uncertainty.
- Phase-graph planted truth, posterior calculation, ambiguity/no-call behaviour, permutation invariance and responses to perturbed inputs.
- Rejection of malformed or scientifically inconsistent input records.
- Isolated execution and consistent reporting of later test groups when an earlier assertion fails.

Supporting evidence and reproducible generation/validation scripts are under [`authoring/evidence/`](authoring/evidence/) and [`authoring/provenance/`](authoring/provenance/). Some validation compares a tuple-based reference formulation with an independently written integer-bitmask implementation. Such agreement is useful evidence of computational consistency, but it does not establish clinical validity or guarantee acceptance by an external evaluation platform.

## Data provenance and interpretation

| Data category | Interpretation |
| --- | --- |
| Human EDA reference/accession context | Real biological reference material used to frame sequence and variant interpretation |
| Synthetic coding variants and assay observations | Artificial benchmark inputs, not individual patient measurements |
| Synthetic planted phase and truth labels | Evaluation targets used to check inference behaviour |
| Authoring validation evidence | Records of development checks and simulations, not independent clinical validation |

**Reproducibility note:** The exact meaning of each input field, assay calibration and JSON result is specified in [`instruction.md`](instruction.md) and the accompanying protocol files. Real reference accessions should not be confused with simulated sample identity or disease status.

## Limitations and responsible use

- Results are conditional on the synthetic data-generating assumptions, model choices and documented calibration procedures.
- A benchmark passing its tests does not demonstrate that a variant is pathogenic or that a participant has a genetic condition.
- Real clinical variant interpretation requires appropriate evidence standards, independent review and, where indicated, further laboratory and clinical confirmation.
- External platform evaluation, difficulty assessments and container verification must be reported from actual runs rather than inferred from local documentation.
- This repository may contain verifier-side fixtures or assessment materials. Check permission and remove restricted material **from Git history as well as the latest tree**, where necessary, before public release.

## References and further reading

- [NCBI Gene — *EDA* (Gene ID 1896)](https://www.ncbi.nlm.nih.gov/gene/1896)
- [NCBI Nucleotide — EDA-A1 transcript, NM_001399.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_001399.5)
- [NCBI Nucleotide — EDA-A2 transcript, NM_001005609.2](https://www.ncbi.nlm.nih.gov/nuccore/NM_001005609.2)
- [HGVS sequence variant nomenclature recommendations](https://hgvs-nomenclature.org/)

The above resources supply biological and nomenclature context; they do not validate the synthetic assay outcomes.

## Author

**Mugambi Isaiah**  
*MSc Bioinformatics and Genomic Medicine | Computational Genomics | Molecular Biology*

Areas of interest: next-generation sequencing, variant interpretation, RNA sequencing, genomic medicine and reproducible computational biology.

- [GitHub — Mugambi24](https://github.com/Mugambi24)
- [LinkedIn — Mugambi Isaiah](https://www.linkedin.com/in/mugambi-isaiah-230593150/)

---

**Research-use disclaimer:** This is a scientific-computing and educational project using synthetic case and assay data. It is not intended for clinical diagnosis or independent medical decision-making.