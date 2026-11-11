# EDA A1/A2 transcript consequences and RNA evidence

## Scientific scope

This reproducible scientific-computing task reconciles cis-phased synthetic EDA coding edits, two EDA transcript products, quality-filtered junction molecules, calibrated isoform mixtures, a pinned public EDA variant crosswalk, batch-specific RNA mixtures, and paired-donor longitudinal RNA signals, and calibrated long-range phase reconstruction from noisy multi-batch pair links. All simulated participant genotypes, read counts, RNA mixtures and truth labels are synthetic, not clinical measurements. The real human EDA locus and public EDA accessions are contextual, not a representation of these synthetic participants.

The task must be completed using only five standard-library Python modules named `cli.py`, `coding.py`, `evidence.py`, `inference.py`, and `audit.py`. The solver processes ten inputs and produces eight JSON outputs. See `instruction.md` for the precise CLI, output contracts, sequence-coordinate conventions, stop precedence, and resource limits. The four assay protocol files under `environment/data/` give the batch, phased, longitudinal and phase-graph contracts. The reference implementation is under `solution/`; an independently implemented alternative is under `authoring/provenance/independent_submission/`. Test truth is withheld from the agent environment.

## Scientific methods and key edge cases

* Sequence reconstruction uses the declared zero-based, genomic-plus variant `start`, with insertions immediately before that coordinate; reverse-complemented exons are assembled in transcript order. Duplicate DNA case identifiers, duplicate edit starts (including paired insertions at one cut), overlapping edits and bad alleles have explicit failure/rejection rules. Disrupting a native termination codon is `stop_lost` even if a later stop shifts to the original codon position; an intact stop shifted by an upstream in-frame insertion is not lost.
* Junction evidence filters on mapping and base quality, selects a reproducible best read per UMI, then distinguishes uniquely assigned and shared reads. Fixed spike-ins calibrate two-component likelihood fits and profile intervals. The likelihood verifier checks the internal delta identity, six-decimal rounding and `[0,1]` domain of both estimates and endpoints, not only proximity to expected outputs.
* Curated crosswalk records reconcile the GRCh38 chromosome accession and its exact version with both transcript HGVS accessions/versions, positions, bases and encoded peptide substitution. Gene interval membership is not a claim of real genomic allele validation.
* Replicated RNA standards allow nonlinear detection bias, contaminated aliquots and calibration uncertainty. Phased deconvolution requires at least four *distinct* informative mixed compositions, including valid three-product mixtures. Longitudinal inference uses four matched panels per donor and separates time, treatment, batch capture and cross-mapping effects.


### Long-range calibrated phase graph — substantive difficulty extension

The eighth output `phase_graph.json` requires joint inference of 10–12 synthetic heterozygous sites per case. Each molecule-link replicate identifier must be unique across the entire assay, not just within one case; the generator and the verifier now enforce this together. Each panel contains multiple contradictory cis/trans links across cycles, three imperfect measurement batches, and different link-class misclassification rates. In one batch, cis calls are inversely associated with true cis phase, so naive majority-rule and simple local spanning-tree assignments are unreliable. The physical orientation is anchored at the first site, but the remaining sites must be inferred globally from the link counts. The required marginal posterior and whole-configuration MAP probability demand uncertainty propagation across connected cycles. Weakly identified sites have null phase bits rather than fabricated calls.

Ground truth is planted using `authoring/provenance/generate_phase_graph.py` with explicit independent seeds (2026110911 for practice, 2026110913 for sealed, 2026110917 for hidden). All case phase vectors are restricted to verifier-side `tests/*_phase_graph_truth.json`; agent input has no bit labels. The independently coded integer-bitmask solver, implemented under `authoring/provenance/independent_submission`, provides a separate computational check from the reference's tuple-based exact enumeration. Five dedicated verifier tests compare posterior probabilities with an independent direct likelihood calculation, check sealed and hidden phase truth, test uncertainty/no-call behavior, permutation invariance, response to changed calibration/link counts, and explicit rejection of invalid assay records. A fixed-zero phase baseline fails on a substantial fraction of sites, and fixed 0.5 probabilities fail against the posterior on strongly identified sites. Case-level discrimination is tested too: every graded panel now has at least one nonanchor site whose posterior favors phase bit 1 (probability above 0.5), so no individual panel can be cleared by a near-certain all-zero answer. The generator reproducibly replaces a rare degenerate all-zero planted panel using a separately seeded phase draw and resampled link counts; the remaining panels, calibration standards and agent-visible output schema are unchanged. In particular, hidden case h09 now carries four nonzero planted sites rather than eleven zeros. The verifier checks this per-case property and records it in the independent validation evidence. The independent simulations and planted-truth error statistics are reproducible with `python3 authoring/provenance/validate_phase_graph.py` and retained in `authoring/evidence/phase_graph_validation_2026-11-09.json`. These checks do not adjust the prior RNA gates or tighten existing tolerances.

## Protein audit schema and full-suite fairness

The `audit.json` field `missing_a2_residues` represents the two missing
EDA-A2 residues in order. Both `"Val307-Glu308"` (a JSON string) and
`["Val307", "Glu308"]` (a two-element JSON string array) represent the same
valid result. Reversed, omitted, or incorrect residues are rejected. Both
encodings are documented in the agent instructions; scientific protein lengths,
full hashes and remaining audit fields are still checked exactly.

The two equivalent `missing_a2_residues` representations are checked by the current `tests/test_workflow.py::test_authentic_protein_audit`, while the complete current-run results are recorded in `authoring/evidence/verifier_budget_alignment_2026-11-10.json`. Older, shorter test-run receipts are excluded from this submission to avoid confusing historical observations with the nine-group, 72-test release.

The isolated verifier evaluates all 72 tests even if an early group reports a
failed assertion. A failed test still yields reward 0, but later scientific
results appear in JUnit and CTRF, preventing a single serialization mistake
from hiding accuracy results. Runner timeouts and incomplete groups still
fail verification. No accuracy, uncertainty, coverage or calibration gate was
weakened by this change.

## Verifier coverage

The 72 named scientific tests also include direct adversarial inputs for paired insertions and a substitution at the same genomic start, accession-version-only transcript HGVS changes, and structurally invalid likelihood output ranges/deltas. They check both the sealed and hidden synthetic cases, exact batch molecule totals, partial/no-call nullability, calibration response, 90th-percentile error, independent truth coverage, width bounds, input permutations and perturbation sensitivity. The runner caches baseline sealed/hidden outputs in a private, fresh verifier-owned directory across pytest groups; this avoids ten repeated solver runs and enforces at most 100 invocations and 900 seconds of cumulative solver execution. An unfinished starter should receive no credit. No solver is permitted to access the sealed truths.

A passing local reference suite does not by itself establish acceptance by an external quality, easiness, difficulty or run-audit stage. All validation results below are limited to the synthetic fixtures and generators supplied with this task.

## Recomputed current-fixture evidence — 10 November 2026

All figures below were recomputed from the actual supplied synthetic fixtures; do not reuse numerical claims from earlier versions. Run `python3 authoring/provenance/phased_model_audit.py`, `python3 authoring/provenance/validate_longitudinal_two_methods.py`, and `python3 authoring/provenance/validate_batch_assay.py` to regenerate independently calibrated results.

### Phased mixture results

| Method / fixture | Fraction MAE | Maximum fraction error | Interval coverage | Mean interval width |
|---|---:|---:|---:|---:|
| Reference EM / sealed | 0.015974 | 0.082760 | 91 / 92 | 0.100753 |
| Reference EM / hidden | 0.014987 | 0.057440 | 112 / 112 | 0.096999 |
| Independent least squares / sealed | 0.016164 | 0.079274 | 90 / 92 | 0.101545 |
| Independent least squares / hidden | 0.014707 | 0.053635 | 112 / 112 | 0.098232 |

### Longitudinal calibrated treatment-by-time effects

The reference EM estimator and an independently coded constrained least-squares estimator each passed the seven published accuracy, coverage and width gates on 22/22 seeded synthetic datasets (11 seed settings × sealed/hidden). Truth was generated separately; neither estimator has access to the answer keys.

| Worst/minimum across 22 datasets | Multinomial EM | Independent least squares | Published bound |
|---|---:|---:|---:|
| Largest case MAE | 0.042229 | 0.043211 | < 0.060 |
| Largest case P90 error | 0.079023 | 0.077172 | < 0.105 |
| Largest donor MAE | 0.087355 | 0.085979 | < 0.115 |
| Minimum case interval coverage | 0.869565 | 0.869565 | ≥ 0.800 |
| Minimum donor interval coverage | 0.782609 | 0.785714 | ≥ 0.750 |
| Largest mean case width | 0.187881 | 0.192509 | < 0.310 |
| Largest mean donor width | 0.272418 | 0.272466 | < 0.480 |

Both positive and negative longitudinal effects are represented, with case-level true magnitudes stratified far enough from zero and fixed ±0.25 guesses to expose those simple answers. A global median can still be numerically close to a particular case mean by coincidence: therefore an explicit *within-case treated-followup counterfactual* checks that every sealed case changes when its corresponding junction reads change. The case-response check is in addition to donor-level accuracy and interval calibration, not a replacement for them. These are synthetic benchmarking results, not clinical uncertainty calibration.

Fresh results and planted-truth provenance are in `authoring/evidence/phased_multivariate_audit.json`, `authoring/evidence/longitudinal_two_method_neutrality.json`, and `authoring/evidence/current_submission_validation_2026-11-10.json`.


The phased-RNA allelic contrast is the difference of two fractions in [-1,1]; the `variant_linked_shift` is a difference of case and control contrasts and therefore has the mathematically valid domain [-2,2], not [-1,1]. The longitudinal effect is a treatment-by-time difference of contrasts and has a different [-4,4] theoretical domain.

## Difficulty

Clinical molecular geneticists and genomics bioinformaticians encounter analogous transcript-variant reconciliation, while RNA-seq analysts and statistical geneticists perform batch and donor-effect inference. A successful solution requires integrating five distinct computational steps. The negative-strand genome uses an insertion-before-base coordinate convention and two different coding exon models; naive reference-index edits change the output hashes. UMI-level read-quality selection, cross-mapped splice reads, mismatched capture efficiencies, one contaminated standard and constrained four-product mixtures each introduce real identifiability or calibration decisions. The paired intervention further mixes a donor's two time points with unequal batch assignment and needs a difference-in-differences estimator, not pooled counts. These are substantive modeling requirements, not arbitrarily narrow accuracy cutoffs. A clear implementation must remain reproducible, correctly report no-calls and quantify calibration uncertainty. The scientific simulation deliberately avoids revealing consequence labels in agent-visible variants, and both case and donor effects are graded against concealed independent truth.

## Reference solution

The five Python standard-library modules implement the complete workflow. `coding.py` applies validated cis-phased genomic edits before transcript-specific exon projection and codon translation, tracks whether the native stop survives, and emits SHA-256 digests and classified consequences. `evidence.py` filters qualities and deduplicates UMIs to derive unique, shared and insufficient RNA evidence. `audit.py` checks the anchored protein sequence, synthetic coordinate disclaimer, transcript accessions and public HGVS crosswalk. `inference.py` combines spike-in-calibrated two-component profile likelihood, nonlinear multi-batch standards, four-product allele/isoform EM deconvolution and paired-donor time-by-treatment contrasts, plus exact Bayesian phasing on noisy signed cycles. `cli.py` validates the ten inputs and writes eight named JSON results with fixed schemas. The reference algorithm is one permitted method, not a required implementation; an independently coded constrained least-squares estimator provides a method-neutral check for the RNA stages. The crosswalk is contextual and never claims synthetic alleles were found in patients.

## Verification

The offline verifier checks the outputs using two isolated sets of seeded, synthetic sealed/hidden fixtures and multiple deliberately corrupted or perturbed inputs. There are 72 named tests. Identity and hash comparisons are exact where outputs are deterministic; statistical RNA estimates are checked against independent latent mixture/donor truth, not forced to equal the reference coefficients. Coding answer keys were initially written by the reference oracle `solution/coding.py:compute`, then independently checked by reconstructing edited genomic-plus exons, reverse complementing in transcript order and comparing coding and translated-protein lengths/hashes using Biopython; see `authoring/provenance/validate_coding_keys.py` and `authoring/evidence/coding_key_provenance_2026-11-09.json`. The sealed expected files are pinned by hash in the source manifest. Neither keys nor latent truth appear in the agent input environment.

Batch calibration was evaluated against independently planted fractions across 12 seeded fixture/replicate panels per estimator, with five numeric seed streams reused across the sealed/hidden populations. Both estimation methods passed all 12 under the original worst-error (<0.065), mean-error (<=0.035), ≥85% coverage and interval-width (<0.45 maximum, <0.32 mean) limits. The reference worst maximum/mean errors were 0.050782/0.021722; the independent method's were 0.054635/0.016088. Both achieved 100% minimum coverage. Uncalibrated raw pooling with fixed ±0.11 intervals fails every sealed/hidden case; the smallest worst-panel raw errors were 0.150188 and 0.147921, respectively. The lane-B calibration response maximum was 0.029504. All statistics are measured on planted synthetic truth and reproduced by `authoring/provenance/validate_batch_assay.py`; no accuracy or uncertainty limit was tightened.

Phased assay and longitudinal thresholds are supported by separate EM and constrained least-squares simulations summarized above. These gates jointly require point accuracy, coverage, width and appropriate response to perturbed input; they cannot be met by emitting a reference truth file or a fixed central value. For longitudinal effects the revised sealed/hidden case and donor truths exclude zero and fixed ±0.25 baselines by design; the code asserts these truth margins without exposing them to the agent. Output schema, strict JSON integer phase bits, six-decimal posterior probabilities, no-call semantics and process isolation are tested separately. Calibration claims concern only the supplied synthetic generator, not clinical accuracy. The reference must pass and an unfinished starter must fail. The test runner's source of truth is `tests/test_workflow.py`, and method-neutral evidence is regenerated by the provided authoring scripts.

### Independent calibration of perturbation and robustness gates (10 November 2026)

The verifier also checks whether estimates respond to *specific, documented* changes in the supplied synthetic assay inputs, and whether one corrupted calibration replicate causes excessive drift. These are operational sensitivity/robustness checks, separate from the planted-truth accuracy, coverage and interval-width gates above. Their cutoffs do not describe patient-level accuracy. The seven verifier-enforced perturbation inequalities, plus one separately labeled non-gating hap0 diagnostic, are independently replayed by `python3 authoring/provenance/validate_perturbation_bounds.py`; the measured results, input hashes and reproduction instructions are in `authoring/evidence/perturbation_bounds_validation_2026-11-10.json`. The script directly runs two distinct estimators without importing the verifier or matching one estimator's answers to the other.

| Test and metric, sealed synthetic assay | Required bound | Reference observed | Independent observed | Reason for this check |
| --- | ---: | ---: | ---: | --- |
| Batch lane-B spike-in change, maximum case-fraction response | >0.004 | 0.029504 | 0.018098 | Excludes an implementation ignoring lane-specific standards; both estimators' actual responses exceed the minimum by more than 0.014. |
| Batch single contaminated 0.5 standard, maximum case-fraction drift | <0.030 | 0.029577 | 0.028427 | Limits damage from one deliberately corrupted replicate, while both robust methods pass. The reference margin is only 0.000423, so this is a narrow fixture-specific bound, **not** proof of stability for all possible contaminants. |
| Phased mixed-standard fraction change, maximum hap1 case-fraction response | >0.004 | 0.179941 | 0.179796 | Excludes ignoring mixed-standard compositions; the current measured response is approximately 0.180 in both independently implemented methods. |
| Phased reassignment of h1_A2 junction reads to h1_A1, signed hap1 response | >0.035 | 0.070078 | 0.070616 | A real hap1-specific junction perturbation must increase the inferred hap1 A1 fraction; unchanged estimates fail. |
| Phased one contaminated `a1_rich` replicate, maximum hap0 drift | <0.050 | 0.000016 | 0.000000 | A robust multi-replicate calibration should not be dominated by one damaged standard; both independent methods have much smaller observed drift. |
| Phase-graph lane-A calibration relation reversal, summed marginal response | >0.10 | 5.726118 | 5.726118 | Excludes fixed phase-graph probabilities that ignore the batch-specific reference links. |
| Phase-graph cis/trans edge-count reversal, maximum marginal site response | >0.03 | 1.000000 | 1.000000 | Excludes fixed answers that ignore molecule-link evidence; a strong graph response is observed from both exact-enumeration methods. |


**Additional diagnostic, not a verifier gate:** The same hap1 junction reassignment also changes the hap0 fraction by 0.021352 in the reference method and 0.012826 in the independent method. These measurements can help assess cross-haplotype leakage, but `tests/test_workflow.py::test_phased_joint_allele_splice_read_response` asserts only the per-case hap1 increase (>0.03); it does **not** enforce a hap0 drift bound such as <0.04. This observation therefore does not affect grading and must not be counted among the seven enforced inequalities.

For every *positive* response gate, a deliberately frozen-output method has measured response **zero** and fails; the independently implemented calibrated methods both pass. Upper drift gates are not sufficient alone—an unresponsive method could also have zero drift—so the verifier simultaneously applies planted-truth accuracy, uncertainty and the positive-response checks. These measured fixtures demonstrate separation from that no-response defect; they do not establish universal optimal cutoffs. The tight 0.030 batch robustness bound is explicitly recorded rather than presented as a generous safety margin. The calibration report's `input_sha256` entries pin the exact assays used. The seven verifier-enforced perturbation checks are distinct from the separate casewise baseline-discrimination checks below; the measured responses are retained and the revised datasets independently rechecked.

### Phase-graph score calibration and independent methods

The agent-visible phase-graph protocol explicitly states four withheld-data gates: absolute marginal posterior error <0.03 at each site, absolute complete-MAP posterior error <0.025 per case, average absolute marginal error <0.01 per dataset, and average Brier score <0.08 per dataset. With planted binary phase bit b and reported marginal probability p, Brier = mean((p-b)^2) over all loci of each separately evaluated withheld dataset. The first three are computational-agreement requirements against a separately implemented exact Bayes posterior, whereas Brier directly measures probabilistic accuracy against generated truth. The planted truth is isolated from the agent. The six-decimal tuple-enumeration reference and separately written integer-bitmask implementation agreed exactly on all 560 marginals and 51 MAP posterior values in the 23 sealed and 28 hidden cases (maximum observed discrepancies 0.000000). The 0.03 and 0.025 allowances are deliberately conservative ceilings rather than claimed empirical optimizer errors; they accommodate different exact numeric summation strategies without requiring identical floating-point implementations. The <0.01 dataset mean error protects against a consistent small bias.

In the reproducible seeded fixtures the independent method had Brier scores 0.000999 (sealed) and 0.014686 (hidden), while reporting p=0 for every locus would score approximately 0.455285 and 0.417197, respectively; p=0.5 everywhere scores 0.25. In a separate 32-seed simulation (seeds 1200 through 1231, 16 cases per seed), both independent inference methods agreed to six decimals; the largest dataset-level mean Brier was 0.043993. Thus the <0.08 gate allows variation above the worst observed calibration dataset while excluding fixed probabilities, but it is a benchmark-specific simulation threshold, not evidence of clinical calibration. The source, seed-wise scores, sample counts, hashes and recalculation command are in `authoring/provenance/calibrate_phase_graph.py` and `authoring/evidence/phase_graph_calibration_2026-11-09.json`.

### Likelihood-key generation and numerical tolerance

The synthetic sealed and hidden likelihood answer keys are reproducible without any hidden manual values. Run `python3 authoring/provenance/validate_likelihood_keys.py` from the task root. The script rebuilds UMI evidence from `tests/{sealed,hidden}_observations.json` using `tests/independent_evidence.py:independently_verify`, confirms the pinned evidence files, and regenerates likelihood outputs from those counts and `tests/assets/spikein_calibration.json` using `tests/independent_calibration.py:calculate` (derivative-root bisection). These results reproduce all 23 sealed and 28 hidden likelihood-key records exactly, including interval endpoints, likelihoods and deltas. The separate reference implementation `solution/inference.py:infer` uses golden-section search and independently agrees to within 0.000001 across the original fixtures. Across 408 further seeded, count-perturbed case/control panels, its maximum disagreement with the derivative implementation is also 0.000001. The sensitivity seed, source-file SHA-256 hashes and comparisons are in `authoring/evidence/likelihood_key_provenance_2026-11-09.json`.

The verifier's absolute `3e-6` tolerance for displayed estimates, log-likelihoods, interval endpoints and deltas corresponds to three units of the specified six-decimal output scale. It allows a small margin for numerical convergence and independent optimizers beyond their observed 0.000001 disagreement; it is not a relative or biological accuracy tolerance. The comparison also requires six-decimal rounding, integer molecule counts, `[0,1]` fractions/endpoints, interval ordering and internal delta consistency, preventing the tolerance from admitting structurally invalid answers. This calibration is computational and applies only to the supplied synthetic assays.

### CTRF reporting on failed verifier runs

`tests/test.sh` initializes reward 0 and invokes the bounded pytest runner to emit JUnit. The trusted `tests/run_isolated.py` executes the complete 72-name manifest in nine isolated pytest subprocess groups (20, 13, 1, 1, 4, 9, 11, 8, 5), merges their JUnit results and prevents sporadic long-running single-process pytest stalls. An unconditional shell EXIT trap calls `tests/ctrf_report.py` to create `/logs/verifier/ctrf.json` on normal shell exit, including when pytest fails, is interrupted by a normal exception, JUnit is absent or malformed, or final consistency/process checks fail. When JUnit has failed test cases, each one is recorded with failed status and its diagnostic; a missing/invalid report generates an explicit failed execution record. Reward 1 is still awarded only after all 72 named tests pass and process-hygiene checks succeed. The reporting implementation handles passing, failing, missing, malformed and failed-gate JUnit results; these are reporting-path behaviors rather than evidence of scientific correctness. `tests/ctrf_report.py` and `tests/test.sh` are the executable sources of truth. The current 72-test run evidence is `authoring/evidence/verifier_budget_alignment_2026-11-10.json`. Forced process termination outside shell control cannot be guaranteed to run a trap.



### Verifier rejection, CTRF and perturbation semantics (current revision)

Invalid-input tests require an explicit nonzero process exit (`WorkflowNonzeroExit`), not an arbitrary `AssertionError`: agent timeouts, zero-exit missing JSON, stale outputs and cleanup failures cannot satisfy rejection controls. The 35-second per-invocation limit remains enforced separately. Across all nine groups, a new root-only baseline cache prevents repeated default sealed/hidden runs, and a trusted ledger enforces no more than 100 actual solver invocations and at most 900 seconds of total solver wall time. The verifier outer timeout is 1100 seconds to allow Python setup, fixture copying, cleanup and CTRF reporting; isolated subprocess fail-safe limits are 950 seconds rather than independent 130-second budgets. These are not additional solver allowances. The lane-B mixed-standard longitudinal perturbation requires a nonzero change exceeding one six-decimal unit in at least one case. Separately, the treated-followup junction-exchange test requires a change greater than 0.03 in each sealed case, as specified in `environment/data/longitudinal_protocol.txt`; both independent estimators exceed this casewise bound. Unaffected cases are graded against independently planted truth using the published mean-error, interval-coverage and width criteria. Exact unchanged JSON is not required: shared-variance estimation may legitimately alter unmodified intervals.

The CTRF converter runs from the verifier shell EXIT trap on successful and unsuccessful runs. Its `results.summary` includes epoch-millisecond `start` and `stop` fields and each `results.tests` item includes `duration` in milliseconds, parsed from the JUnit testcase time attribute. These values are also present for explicit failure/diagnostic records (duration 0). Missing/invalid JUnit, failed tests and failed final reward checks therefore remain visible, without turning a failed run into reward 1. The numerical calibration eligibility thresholds are specified in `environment/data/batch_assay_protocol.txt` and match `solution/inference.py:reconcile_batches`: standards need at least 100 unique reads per accepted replicate, two accepted replicates at four fractions spanning 0.8, a median response range of at least 0.25 across the extreme standards, and 120 panel-specific unique reads to retain each measured batch.

## Provenance and reproduction

* Generators: `authoring/provenance/generate_fixture.py`, `generate_observations.py`, `generate_batch_assay.py`, `generate_phased_assay.py`, `generate_longitudinal_assay.py`, `generate_phase_graph.py`. The input generator uses pinned random seeds; avoid replacing the supplied sealed/hidden fixtures during verifier runs.
* Science-validation scripts: `authoring/provenance/phased_model_audit.py`, `validate_longitudinal_two_methods.py`. An alternative solver is retained under `authoring/provenance/independent_submission/` for genuinely separate method checks.
* Expected coding keys: the reference oracle generated sealed and hidden expected JSON; `authoring/provenance/validate_coding_keys.py` independently checks 92 resolved CDS/protein records by exon reconstruction and codon-table translation, and the manifest pins their digests.
* Expected likelihood keys: `authoring/provenance/validate_likelihood_keys.py` reproduces them from seeded observations, independently grouped UMI evidence and spike-in controls with the derivative-bisection estimator; the distinct golden-section reference supplies the numerical comparison and tolerance sensitivity evidence.
* Batch calibration: `authoring/provenance/validate_batch_assay.py` reproduces twelve fixture/seed panels per estimator, with five reused numeric seed streams across both fixture populations and the measured lane-B perturbation response. `authoring/provenance/validate_submission.py` checks leakage, donor-level shortcuts, planted-truth hashes, declared manifest hashes and the 72-test manifest.
* Source provenance: `authoring/provenance/source_manifest.json` and `authoring/provenance/public_variant_snapshot.json`; all biological accession claims are pinned to the source record versions in these files. Optional regeneration dependencies are documented in `authoring/provenance/requirements.txt`. Agent deliverables depend on standard-library Python only.
* The verifier itself is in `tests/test_workflow.py`; `tests/test.sh` checks completion and process hygiene before assigning reward. The task's project metadata is in `task.toml`, following the mandatory bundle-structure schema.

## Public scientific references

- [NCBI EDA Gene ID 1896](https://www.ncbi.nlm.nih.gov/gene/1896)
- [EDA-A1 protein NP_001390.1](https://www.ncbi.nlm.nih.gov/protein/NP_001390.1)
- [EDA-A1 transcript NM_001399.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_001399.5)
- [EDA-A2 transcript NM_001005609.2](https://www.ncbi.nlm.nih.gov/nuccore/NM_001005609.2)
- [ClinVar RCV006905138.1](https://www.ncbi.nlm.nih.gov/clinvar/RCV006905138.1/)

These source links are contextual annotations; the artificial transcript coordinates, sequence edits and RNA measurements used for this exercise do not establish patient diagnoses or genuine biological findings.

### Current verifier-budget and anchor review (10 November 2026)

The changed FASTA/genome mutation probe checks integrity of the canonical EDA-A1/A2 protein anchor as an input-wide consistency condition. A changed anchor is rejected nonzero; a reference-allele disagreement in one edit remains a per-variant `reference_mismatch` record. The batch perturbation requires movement in the biologically specified direction but no undisclosed >0.07 fraction jump. The lane-B mixed-standard perturbation must change at least one rounded longitudinal case effect by more than one six-decimal unit. The distinct treated-followup junction-exchange perturbation must change each sealed case effect by >0.03; its current calibrated minimum response across two methods and withheld fixtures is 0.268244. No additional ten-case >0.02 threshold is enforced. The separate `invocation_accounting.json` report gives the observed invocation count and solver time for each verifier attempt. The synthetic batch-calibration evidence comprises two fixture populations and five numeric seed streams reused across them; it is not twelve independent RNG streams.

### Independent verifier-budget and per-haplotype validation (10 November 2026)

The current reference and independently implemented solvers each passed the 72-test manifest from the regenerated casewise fixtures with 100 actual agent invocations. The current measured solver times were 50.743094 seconds for the reference and 34.181586 seconds for the independent method, both within the 900-second cumulative solver budget. The nine bounded group sizes are 20, 13, 1, 1, 4, 9, 11, 8 and 5; current results are pinned in `authoring/evidence/verifier_budget_alignment_2026-11-10.json`. The separate verifier timeout is 1100 seconds to cover trusted setup and reporting; its subprocess fail-safe timeout does not reduce the agent budget. Baseline result reuse is scoped to a fresh, root-owned directory per verifier attempt, never accessible to the solver UID.

The refreshed `authoring/evidence/phased_multivariate_audit.json` now computes the phased-RNA MAE, maximum error, number of marginal errors ≥0.08, interval coverage and width limits **separately for hap0 and hap1** in both sealed and hidden fixtures for each of two different estimators. All eight method/fixture/haplotype combinations satisfy the published gates. The exact current observations, hashes, group sizes and commands are recorded in `authoring/evidence/verifier_budget_alignment_2026-11-10.json`.

### Current casewise discrimination and input-boundary validation (10 November 2026)

The revised nonlinear batch-capture simulation was regenerated from pinned seeds. Per-case raw-pooling controls demonstrate that a method-free A1-unique ratio with fixed uncertainty fails on every sealed and hidden case; independent calibration methods still satisfy all original statistical limits across their twelve calibration datasets. No scientific scoring boundary was tightened. The consequence regression now distinguishes an intact native terminal stop shifted by a same-frame stop insertion from a genuinely deleted native termination site. The input-validation tests explicitly reject duplicated DNA case identifiers, negative batch evidence counts, phase panels with fewer than ten loci, and phase edge replicate identifiers reused across cases. The reference and independent implementations share these input-contract requirements but retain distinct RNA estimation algorithms.

### Current individual-case baseline discrimination and counterfactual calibration — 10 November 2026

The revised synthetic generators now include non-neutral case/control batch shifts, detectable phased-allele case changes, and stronger alternating-sign longitudinal effects, all sampled deterministically without reference-estimator output. Public scientific accuracy, coverage and width thresholds are preserved. The planted truths and assays are independently regenerated using `authoring/provenance/generate_batch_assay.py`, `generate_phased_assay.py`, and `generate_longitudinal_assay.py`; the new parameter ranges are simulation design, not clinical priors.

In each fixture, copying control estimates into case estimates fails the published individual-case batch or phased truth limits (0 cases locally plausible for either shortcut). The fixed longitudinal +0.25, −0.25 and 0.0 case-effect guesses also fail the published case-point-error threshold for every withheld case. A global median or first-case value may be close to an individual longitudinal mean by chance; such constant guesses nevertheless fail the *joint* case-specific response test. Thus each case is checked with a controlled RNA-evidence counterfactual: the same estimate cannot be emitted for the original and changed evidence. Frozen outputs have zero response by definition. The verifier checks the response of the submitted method on the altered inputs, not an expected reference coefficient.

The smallest measured case-level changes, across both independently coded methods and both withheld synthetic fixtures, are:

| Case-level counterfactual | Minimum required response | Smallest measured response | Headroom to cutoff |
|---|---:|---:|---:|
| Batch case-group A1 shift | >0.020 | 0.051487 | 0.031487 |
| Phased hap1 junction reassignment | >0.030 | 0.073524 | 0.043524 |
| Paired-donor followup junction exchange | >0.030 | 0.268244 | 0.238244 |

These measured lower bounds exceed the thresholds without reference-output equality. Constant, copied-control and frozen answers respond by zero. The hidden-case responses shown here were computed by the independent authoring checks; the test manifest performs its stated sealed/hidden accuracy checks and its specified casewise perturbations within the 100-invocation budget.

Three additional **explicitly documented** scientific response gates are contained inside existing named tests, with no extra solver invocations: every batch case must increase its case A1 estimate by >0.02 after moving case-group A2-unique evidence to A1-unique; every sealed phased case must increase hap1 A1 fraction by >0.03 after reassigning up to 80 hap1_A2 junction reads per case library; and every sealed longitudinal case must change the case treatment-by-time effect by >0.03 after exchanging treated-followup hap1 and unphased A1/A2 evidence. Each cutoff is lower than the smallest measured change under both independent correct estimators, while a fixed, copied or no-change method has zero response. Per-method, per-fixture minimum responses and baseline counts are generated into `authoring/evidence/current_submission_validation_2026-11-10.json` by `python3 authoring/provenance/validate_submission.py`. These response checks are synthetic causal sanity checks, not estimates of clinical measurement sensitivity.

## Phase-graph cardinality and RNA rounding audit

The phase-graph verifier tests 9 loci with all removed-site edges deleted, so
only cardinality can explain rejection. A separate well-formed 13-locus panel
probes the upper bound. The reference also rejects 16 cleanly connected loci.

The simple RNA evidence delta accepts two deterministic, biologically
equivalent six-decimal rounding orders: subtract the displayed rounded beta
means, or subtract the unrounded beta means before rounding the difference.
The verifier recomputes the permitted values from checked independent UMI
counts in every evidence comparison, including molecule-order, single-case,
and independent-calculation tests; it does not enforce one convention in a
second test. No broad tolerance is used. Other JSON values, intervals and
classification remain checked exactly. Both conventions are also tested
against a deliberately incorrect three-millionth delta.
