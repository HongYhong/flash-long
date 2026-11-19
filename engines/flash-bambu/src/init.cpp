#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" {
    SEXP bambu_next_probe_abi(SEXP ir_test_sexp, SEXP n_threads_sexp);
    SEXP bambu_next_em_theta(SEXP XSEXP, SEXP YSEXP, SEXP maxiterSEXP, SEXP minvalueSEXP, SEXP convSEXP);
    SEXP bambu_next_emWithL1(SEXP ASEXP, SEXP A_fullSEXP, SEXP A_uniqueSEXP,
                             SEXP YSEXP, SEXP KSEXP,
                             SEXP maxiterSEXP, SEXP minvalueSEXP, SEXP convSEXP);
    SEXP bambu_next_csr_is_valid(SEXP ptr_sexp);
    SEXP bambu_next_csr_cache_lookup(
        SEXP n_tx_sexp, SEXP ex_starts_sexp, SEXP ex_ends_sexp,
        SEXP eq_lens_sexp, SEXP eq_vals_sexp, SEXP first_tx_sexp, SEXP last_tx_sexp);
    SEXP bambu_next_csr_cache_put(
        SEXP xptr_sexp, SEXP n_tx_sexp, SEXP ex_starts_sexp, SEXP ex_ends_sexp,
        SEXP eq_lens_sexp, SEXP eq_vals_sexp, SEXP first_tx_sexp, SEXP last_tx_sexp);
    SEXP bambu_next_csr_cache_clear();
    SEXP bambu_next_prepare_gtf_cpp(SEXP gtf_path_sexp, SEXP n_threads_sexp);
    SEXP bambu_next_extract_bam_and_motifs_cpp(
        SEXP bam_path_sexp, SEXP bai_path_sexp, SEXP run_name_sexp,
        SEXP genome_seqnames_sexp, SEXP genome_seqlengths_sexp, SEXP ann_seqnames_sexp,
        SEXP fa_path_sexp, SEXP fai_path_sexp, SEXP track_reads_sexp, SEXP n_threads_sexp);
    SEXP bambu_next_fill_motifs_from_r_cpp(
        SEXP bsd_ptr_sexp, SEXP seq_names_sexp, SEXP seq_chars_sexp);
    SEXP bambu_next_export_unique_junctions_cpp(SEXP bsd_ptr_sexp);
    SEXP bambu_next_prepare_junc_ec_cpp(
        SEXP bsd_ptr_sexp, SEXP csr_ptr_sexp, SEXP stranded_sexp);
    SEXP bambu_next_build_read_classes_cpp(
        SEXP bsd_ptr_sexp, SEXP csr_ptr_sexp, SEXP stranded_sexp, SEXP can_run_hc_sexp,
        SEXP plus_ss_sexp, SEXP plus_se_sexp, SEXP plus_es_sexp, SEXP plus_ee_sexp,
        SEXP minus_ss_sexp, SEXP minus_se_sexp, SEXP minus_es_sexp, SEXP minus_ee_sexp);
    SEXP bambu_next_get_bsd_active_seqs_cpp(SEXP bsd_ptr_sexp);
    SEXP bambu_next_subset_bsd_seq_cpp(SEXP bsd_ptr_sexp, SEXP seq_1b_sexp);
    SEXP bambu_next_get_bsd_tracked_cpp(SEXP bsd_ptr_sexp);
    SEXP bambu_next_build_csr_from_flat(
        SEXP seq_names_sexp, SEXP seq_lengths_sexp,
        SEXP tx_names_sexp, SEXP tx_gene_names_sexp,
        SEXP tx_seq_1b_sexp, SEXP tx_strand_code_sexp, SEXP tx_ex_counts_sexp,
        SEXP ex_start_sexp, SEXP ex_end_sexp, SEXP ex_rank_sexp, SEXP ex_end_rank_sexp,
        SEXP eq_vals_1b_sexp, SEXP eq_lens_sexp,
        SEXP compute_eq_if_null_sexp, SEXP n_threads_sexp);
    SEXP bambu_next_assign_gene_ids_cpp(
        SEXP q_seq_names_sexp,
        SEXP q_tx_ex_counts_sexp,
        SEXP q_ex_seq_1b_sexp,
        SEXP q_ex_start_sexp,
        SEXP q_ex_end_sexp,
        SEXP q_ex_strand_sexp,
        SEXP csr_ptr_sexp,
        SEXP ann_ex_counts_sexp,
        SEXP ann_gene_names_sexp,
        SEXP ann_ex_seq_1b_sexp,
        SEXP ann_ex_start_sexp,
        SEXP ann_ex_end_sexp,
        SEXP ann_ex_strand_sexp,
        SEXP min_exon_overlap_sexp,
        SEXP fusion_mode_sexp,
        SEXP max_chain_iter_sexp);
    SEXP bambu_next_is_rc_compatible_cpp(
        SEXP q_seq_names_sexp,
        SEXP q_tx_ex_counts_sexp,
        SEXP q_ex_seq_1b_sexp,
        SEXP q_ex_start_sexp,
        SEXP q_ex_end_sexp,
        SEXP q_ex_strand_sexp,
        SEXP q_ex_rank_sexp,
        SEXP q_ex_end_rank_sexp,
        SEXP csr_ptr_sexp,
        SEXP n_threads_sexp);
    SEXP bambu_next_count_polya_faidx_cpp(
        SEXP fa_path_sexp,
        SEXP fai_path_sexp,
        SEXP seq_names_sexp,
        SEXP seq_lengths_sexp,
        SEXP start_seq_1b_sexp,
        SEXP start_s_sexp,
        SEXP start_e_sexp,
        SEXP end_seq_1b_sexp,
        SEXP end_s_sexp,
        SEXP end_e_sexp,
        SEXP strand_code_sexp,
        SEXP n_threads_sexp);
    SEXP bambu_next_kway_merge_spliced_cpp(
        SEXP samples_list_sexp,
        SEXP min_read_count_sexp,
        SEXP min_read_frac_sexp,
        SEXP min_tx_me_sexp,
        SEXP min_tx_se_sexp);
    SEXP bambu_next_merge_unspliced_cpp(
        SEXP seq_levels_sexp,
        SEXP seq_idx_1b_sexp,
        SEXP strand_code_sexp,
        SEXP start_sexp,
        SEXP end_sexp,
        SEXP read_count_sexp,
        SEXP gene_read_prop_sexp,
        SEXP tx_score_sexp,
        SEXP tx_score_nofit_sexp,
        SEXP stranded_sexp,
        SEXP min_read_count_sexp,
        SEXP min_read_frac_sexp,
        SEXP min_tx_single_exon_sexp);
    SEXP bambu_next_build_exons_from_introns_cpp(
        SEXP q_tx_seq_1b_sexp,
        SEXP q_start_sexp,
        SEXP q_end_sexp,
        SEXP q_tx_strand_code_sexp,
        SEXP q_is_vec_sexp,
        SEXP q_ie_vec_sexp);
    SEXP bambu_next_add_new_spliced_cpp(
        SEXP csr_ptr_sexp,
        SEXP q_seq_names_sexp,
        SEXP q_tx_ex_counts_sexp,
        SEXP q_ex_seq_1b_sexp,
        SEXP q_ex_start_sexp,
        SEXP q_ex_end_sexp,
        SEXP q_ex_strand_sexp,
        SEXP q_ex_rank_sexp,
        SEXP q_ex_end_rank_sexp,
        SEXP q_in_counts_sexp,
        SEXP q_in_start_sexp,
        SEXP q_in_end_sexp,
        SEXP min_exon_distance_sexp,
        SEXP min_ps_dist_sexp,
        SEXP min_ps_dist_se_sexp,
        SEXP n_threads_sexp);
    SEXP bambu_next_filter_unspliced_overlap_cpp(
        SEXP csr_ptr_sexp,
        SEXP q_seq_names_sexp,
        SEXP q_seq_1b_sexp,
        SEXP q_start_sexp,
        SEXP q_end_sexp,
        SEXP q_strand_code_sexp,
        SEXP min_overlap_sexp);
    SEXP bambu_next_extend_min_eq_csr_cpp(
        SEXP ref_csr_ptr_sexp,
        SEXP n_novel_sexp,
        SEXP seq_names_sexp,
        SEXP seq_lengths_sexp,
        SEXP tx_names_sexp,
        SEXP tx_gene_names_sexp,
        SEXP tx_seq_1b_sexp,
        SEXP tx_strand_code_sexp,
        SEXP tx_ex_counts_sexp,
        SEXP ex_start_sexp,
        SEXP ex_end_sexp,
        SEXP ex_rank_sexp,
        SEXP ex_end_rank_sexp,
        SEXP n_threads_sexp);
    SEXP bambu_next_calc_dist_to_annotation_cpp(
        SEXP q_seq_names_sexp,
        SEXP q_tx_ex_counts_sexp,
        SEXP q_ex_seq_1b_sexp,
        SEXP q_ex_start_sexp,
        SEXP q_ex_end_sexp,
        SEXP q_ex_strand_sexp,
        SEXP q_ex_rank_sexp,
        SEXP q_ex_end_rank_sexp,
        SEXP csr_ptr_sexp,
        SEXP max_dist_sexp,
        SEXP ps_dist_sexp,
        SEXP ps_dist_se_sexp,
        SEXP ignore_strand_sexp,
        SEXP n_threads_sexp);
    SEXP bambu_next_quantify_cpp(
        SEXP csr_ptr_sexp,
        SEXP q_seq_names_sexp,
        SEXP q_tx_ex_counts_sexp,
        SEXP q_ex_seq_1b_sexp,
        SEXP q_ex_start_sexp,
        SEXP q_ex_end_sexp,
        SEXP q_ex_strand_sexp,
        SEXP q_ex_rank_sexp,
        SEXP q_ex_end_rank_sexp,
        SEXP q_read_count_sexp,
        SEXP q_rc_names_sexp,
        SEXP isore_params_sexp,
        SEXP em_params_sexp,
        SEXP return_dist_table_sexp,
        SEXP n_threads_sexp);
    SEXP bambu_next_write_gtf_cpp(
        SEXP file_path_sexp,
        SEXP seq_names_sexp,
        SEXP tx_names_sexp,
        SEXP tx_gene_names_sexp,
        SEXP tx_ex_counts_sexp,
        SEXP ex_seq_1b_sexp,
        SEXP ex_start_sexp,
        SEXP ex_end_sexp,
        SEXP ex_strand_code_sexp,
        SEXP ex_rank_sexp);
}

static const R_CallMethodDef CallEntries[] = {
    {"bambu_next_probe_abi", (DL_FUNC) &bambu_next_probe_abi, 2},
    {"bambu_next_em_theta", (DL_FUNC) &bambu_next_em_theta, 5},
    {"bambu_next_emWithL1", (DL_FUNC) &bambu_next_emWithL1, 8},
    {"bambu_next_csr_is_valid", (DL_FUNC) &bambu_next_csr_is_valid, 1},
    {"bambu_next_csr_cache_lookup", (DL_FUNC) &bambu_next_csr_cache_lookup, 7},
    {"bambu_next_csr_cache_put", (DL_FUNC) &bambu_next_csr_cache_put, 8},
    {"bambu_next_csr_cache_clear", (DL_FUNC) &bambu_next_csr_cache_clear, 0},
    {"bambu_next_prepare_gtf_cpp", (DL_FUNC) &bambu_next_prepare_gtf_cpp, 2},
    {"bambu_next_extract_bam_and_motifs_cpp", (DL_FUNC) &bambu_next_extract_bam_and_motifs_cpp, 10},
    {"bambu_next_fill_motifs_from_r_cpp", (DL_FUNC) &bambu_next_fill_motifs_from_r_cpp, 3},
    {"bambu_next_export_unique_junctions_cpp", (DL_FUNC) &bambu_next_export_unique_junctions_cpp, 1},
    {"bambu_next_prepare_junc_ec_cpp", (DL_FUNC) &bambu_next_prepare_junc_ec_cpp, 3},
    {"bambu_next_build_read_classes_cpp", (DL_FUNC) &bambu_next_build_read_classes_cpp, 12},
    {"bambu_next_get_bsd_active_seqs_cpp", (DL_FUNC) &bambu_next_get_bsd_active_seqs_cpp, 1},
    {"bambu_next_subset_bsd_seq_cpp", (DL_FUNC) &bambu_next_subset_bsd_seq_cpp, 2},
    {"bambu_next_get_bsd_tracked_cpp", (DL_FUNC) &bambu_next_get_bsd_tracked_cpp, 1},
    {"bambu_next_build_csr_from_flat", (DL_FUNC) &bambu_next_build_csr_from_flat, 15},
    {"bambu_next_assign_gene_ids_cpp", (DL_FUNC) &bambu_next_assign_gene_ids_cpp, 16},
    {"bambu_next_is_rc_compatible_cpp", (DL_FUNC) &bambu_next_is_rc_compatible_cpp, 10},
    {"bambu_next_count_polya_faidx_cpp", (DL_FUNC) &bambu_next_count_polya_faidx_cpp, 12},
    {"bambu_next_kway_merge_spliced_cpp", (DL_FUNC) &bambu_next_kway_merge_spliced_cpp, 5},
    {"bambu_next_merge_unspliced_cpp", (DL_FUNC) &bambu_next_merge_unspliced_cpp, 13},
    {"bambu_next_build_exons_from_introns_cpp", (DL_FUNC) &bambu_next_build_exons_from_introns_cpp, 6},
    {"bambu_next_add_new_spliced_cpp", (DL_FUNC) &bambu_next_add_new_spliced_cpp, 16},
    {"bambu_next_filter_unspliced_overlap_cpp", (DL_FUNC) &bambu_next_filter_unspliced_overlap_cpp, 7},
    {"bambu_next_extend_min_eq_csr_cpp", (DL_FUNC) &bambu_next_extend_min_eq_csr_cpp, 14},
    {"bambu_next_calc_dist_to_annotation_cpp", (DL_FUNC) &bambu_next_calc_dist_to_annotation_cpp, 14},
    {"bambu_next_quantify_cpp", (DL_FUNC) &bambu_next_quantify_cpp, 15},
    {"bambu_next_write_gtf_cpp", (DL_FUNC) &bambu_next_write_gtf_cpp, 10},
    {NULL, NULL, 0}
};

extern "C" void R_init_bambu(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
