#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: monitor_process_tree.sh ROOT_PID PHASE_FILE OUTPUT_TSV [INTERVAL_SEC]" >&2
    exit 2
fi

root_pid=$1
phase_file=$2
output_tsv=$3
interval=${4:-1}

start_ns=$(date +%s%N)
printf 'elapsed_seconds\tphase\ttree_rss_kib\ttree_pss_kib\ttree_pss_anon_kib\ttree_pss_file_kib\tprocess_count\tsample_span_ms\n' > "${output_tsv}"

while kill -0 "${root_pid}" 2>/dev/null; do
    t0_ns=$(date +%s%N)
    elapsed=$(awk -v s="${start_ns}" -v n="${t0_ns}" 'BEGIN { printf "%.3f", (n - s) / 1e9 }')
    phase="run"
    if [[ "${phase_file}" != "-" && -s "${phase_file}" ]]; then
        phase=$(tr -d '\r\n\t ' < "${phase_file}")
    fi

    mapfile -t tree_snapshot < <(
        ps -eo pid=,ppid=,rss= | awk -v root="${root_pid}" '
            {
                parent[$1] = $2
                rss[$1] = $3
            }
            END {
                keep[root] = 1
                changed = 1
                while (changed) {
                    changed = 0
                    for (pid in parent) {
                        if (!keep[pid] && keep[parent[pid]]) {
                            keep[pid] = 1
                            changed = 1
                        }
                    }
                }
                total = 0
                count = 0
                for (pid in keep) {
                    if (keep[pid] && pid in rss) {
                        total += rss[pid]
                        count += 1
                    }
                }
                print total + 0, count + 0
                for (pid in keep) {
                    if (keep[pid] && pid in rss) print pid, rss[pid]
                }
            }
        '
    )
    if [[ ${#tree_snapshot[@]} -eq 0 ]]; then
        break
    fi

    read -r ps_rss_kib process_count <<< "${tree_snapshot[0]}"
    tree_rss_kib=0
    tree_pss_kib=0
    tree_pss_anon_kib=0
    tree_pss_file_kib=0

    for ((i = 1; i < ${#tree_snapshot[@]}; i++)); do
        read -r pid pid_ps_rss <<< "${tree_snapshot[$i]}"
        if [[ -r "/proc/${pid}/smaps_rollup" ]]; then
            read -r r_rss pss pss_anon pss_file < <(
                awk '
                    $1 == "Rss:" { r = $2 }
                    $1 == "Pss:" { p = $2 }
                    $1 == "Pss_Anon:" { a = $2 }
                    $1 == "Pss_File:" { f = $2 }
                    END { printf "%d %d %d %d\n", r + 0, p + 0, a + 0, f + 0 }
                ' "/proc/${pid}/smaps_rollup" 2>/dev/null || echo "0 0 0 0"
            )
            if (( r_rss < pss )); then
                r_rss=${pss}
            fi
            tree_rss_kib=$((tree_rss_kib + r_rss))
            tree_pss_kib=$((tree_pss_kib + pss))
            tree_pss_anon_kib=$((tree_pss_anon_kib + pss_anon))
            tree_pss_file_kib=$((tree_pss_file_kib + pss_file))
        else
            tree_rss_kib=$((tree_rss_kib + pid_ps_rss))
        fi
    done
    if (( tree_rss_kib < ps_rss_kib && tree_rss_kib == 0 )); then
        tree_rss_kib=${ps_rss_kib}
    fi

    t1_ns=$(date +%s%N)
    span_ms=$(awk -v a="${t0_ns}" -v b="${t1_ns}" 'BEGIN { printf "%.2f", (b - a) / 1e6 }')

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${elapsed}" "${phase}" "${tree_rss_kib}" "${tree_pss_kib}" \
        "${tree_pss_anon_kib}" "${tree_pss_file_kib}" "${process_count}" "${span_ms}" \
        >> "${output_tsv}"
    sleep "${interval}"
done
