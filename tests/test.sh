#!/bin/sh
set -u
umask 077
mkdir -p /logs/verifier
chmod 700 /logs/verifier
printf '0\n' > /logs/verifier/reward.txt
# Report every verifier outcome, including pytest failures, malformed JUnit and
# process-cleanup/reward-gate failures. The trap runs on every shell exit.
rm -f /logs/verifier/junit.xml /logs/verifier/ctrf.json /logs/verifier/verification_failure.txt
date +%s%3N > /logs/verifier/started_at_ms.txt
trap 'python3 /tests/ctrf_report.py /logs/verifier/junit.xml /logs/verifier/ctrf.json /logs/verifier/verification_failure.txt /logs/verifier/started_at_ms.txt || true' EXIT
cd /app || exit 0
# Assert the dedicated deliverable was transferred; no starter workbench should be visible here.
# The shell test report is authoritative for both Oracle and NOP.
for f in cli.py coding.py evidence.py audit.py inference.py; do
    if [ ! -s "/app/submission/$f" ]; then
        printf 'ARTIFACT HANDOFF ERROR: missing /app/submission/%s\n' "$f" >&2
    fi
done
if [ -f /app/submission/coding.py ] && grep -q 'Coding reconstruction is unfinished' /app/submission/coding.py; then
    printf 'ARTIFACT HANDOFF ERROR: unfinished coding.py received instead of submitted code\n' >&2
fi
if python3 /tests/run_isolated.py /logs/verifier/junit.xml; then
    # A session-changing agent descendant could escape the process-group kill.
    # Sweep the entire demoted uid before the trusted reward gate executes.
    if ! python3 /tests/process_hygiene.py; then
        printf 'Untrusted processes survived verifier cleanup\n' >&2
        printf 'Untrusted processes survived verifier cleanup\n' > /logs/verifier/verification_failure.txt
        exit 0
    fi
    python3 - <<'PY'
from pathlib import Path
import xml.etree.ElementTree as ET
p=Path('/logs/verifier/junit.xml')
root=ET.parse(p).getroot()
suites=[root] if root.tag=='testsuite' else list(root.findall('testsuite'))
assert len(suites)==1
s=suites[0]
cases=list(s.findall('testcase'))
expected={'test_authentic_protein_audit','test_sealed_coding_evidence','test_holdout_coding_evidence','test_output_schemas_and_order','test_known_negative_strand_consequences','test_molecule_coverage_and_uncertainty','test_case_order_independence','test_molecule_order_independence','test_single_case_reconciliation','test_no_unexplained_identity_baseline','test_calibrated_mixture_against_independent_expected','test_likelihood_molecule_counts_require_integer_json_type','test_likelihood_numeric_output_must_be_six_decimal_rounded','test_each_graded_likelihood_instance_requires_estimation','test_calibration_change_is_reflected','test_mixture_low_coverage_is_explicit','test_variants_order_independence','test_low_quality_duplicate_molecules_do_not_change_results','test_ignoring_all_variants_fails_every_coding_record','test_null_and_boundary_cases_are_scientifically_checked','test_deliverables_are_standard_library_only','test_full_order_invariance_all_outputs','test_umi_ranking_with_conflicting_signal_and_qualities','test_shared_molecules_are_not_unique_isoform_support','test_calibration_changes_pseudocount_coverage_and_mle_endpoints','test_invalid_allele_pair_and_coordinate_types','test_real_eda_protein_variant_concordance','test_curated_protein_panel_reorder_and_isoform_deletion','test_curated_variant_mismatch_and_assertion_conflicts','test_authentic_human_eda_grch38_crosswalk','test_authentic_crosswalk_is_evidence_dependent','test_variant_rna_case_id_mismatches_are_rejected','test_curated_accession_assertion_schema_and_provenance','test_curated_crosswalk_exact_schema_uniqueness_and_deduplication','test_reference_protein_anchors_are_mutation_sensitive','test_transcript_hgvs_coordinate_and_allele_reconciliation','test_multi_exon_reason_codes_are_exhaustive','test_projection_mismatch_requires_independently_consistent_cds_anchor','test_stop_lost_even_when_later_inframe_stop_remains','test_replicated_batch_assay_ground_truth','test_batch_assay_input_order_is_invariant','test_batch_specific_calibration_is_used','test_batch_corrupted_calibration_replicate_is_robust','test_batch_two_informative_lanes_are_sufficient','test_batch_rounding_counts_and_unidentifiable_schema','test_batch_neutral_controls_require_measurement_response','test_batch_missing_ids_and_single_lane_are_unidentifiable','test_batch_unidentifiable_and_cross_sample_id_validation','test_phased_multivariate_mixture_truth_and_intervals','test_phased_mixture_input_permutations','test_phased_two_valid_batches_suffice','test_phased_calibration_mix_response','test_phased_joint_allele_splice_read_response','test_phased_single_lane_and_missing_ids_are_not_guessed','test_phased_single_contaminated_standard_robustness','test_phased_pooled_and_fixed_answer_baselines_do_not_solve','test_phased_duplicate_calibration_or_bad_molecule_types_rejected','test_phased_absent_mixed_standards_is_unidentifiable','test_phased_two_informative_batches_can_still_be_unidentifiable','test_longitudinal_paired_truth_and_calibrated_intervals','test_longitudinal_measurement_order_invariance','test_longitudinal_reacts_to_shifted_treatment_junction_calls','test_longitudinal_missing_donor_and_insufficient_case','test_longitudinal_invalid_library_identifiers_and_counts','test_longitudinal_raw_batch_pooling_is_insufficient','test_longitudinal_case_ids_match_dna','test_longitudinal_depends_on_independent_mixed_standard_calibration','test_phase_graph_scientific_posterior_and_schema','test_phase_graph_requires_real_global_inference','test_phase_graph_input_order_invariance','test_phase_graph_responds_to_calibration_and_edge_data','test_phase_graph_unidentifiable_panels_and_bad_inputs'}
assert len(cases)==len(expected)
assert {r.get('name') for r in cases}==expected
assert all(not any(c.find(t) is not None for t in ('skipped','failure','error')) for c in cases)
assert int(s.get('tests',-1))==len(expected)
assert all(int(s.get(k,-1))==0 for k in ('failures','errors','skipped'))
# CTRF is emitted unconditionally by the EXIT trap, never just on success.

PY
    if [ $? -eq 0 ]; then
        # Verify a second time immediately before writing the reward.
        if python3 /tests/process_hygiene.py; then
            printf '1\n' > /logs/verifier/reward.txt
        else
            printf 'Untrusted processes survived final cleanup\n' > /logs/verifier/verification_failure.txt
        fi
    else
        printf 'JUnit reward-gate consistency assertions failed\n' > /logs/verifier/verification_failure.txt
    fi
else
    # Failed runs must also release detached children before returning.
    python3 /tests/process_hygiene.py || true
fi
