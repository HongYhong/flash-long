"""Streaming IntronGraph, IntronCollector, and IntronPath threading with deque BFS."""

from __future__ import annotations

import logging
from collections import defaultdict, deque
from typing import Dict, List, Optional, Tuple

from isoquant_next.core.utils import find_closest, overlaps

logger = logging.getLogger("IsoQuant")


class TerminalVertex:
    polya = -10
    read_end = -11
    tss_right = -12
    polyt = -20
    read_start = -21
    tss_left = -22

    right = (polya, read_end, tss_right)
    left = (polyt, read_start, tss_left)


def is_terminal_vertex(v: Tuple[int, int]) -> bool:
    return v[0] in TerminalVertex.right


def is_starting_vertex(v: Tuple[int, int]) -> bool:
    return v[0] in TerminalVertex.left


class IntronCollector:
    """Collects, clusters, and tracks discarded/corrected splice junctions per locus."""

    def __init__(self, gene_info, delta: int = 0) -> None:
        self.gene_info = gene_info
        self.known_introns = set(self.gene_info.intron_profiles.features)
        self.delta = delta
        self.clustered_introns: Dict[Tuple[int, int], int] = defaultdict(int)
        self.intron_correction_map: Dict[Tuple[int, int], Tuple[int, int]] = {}
        self.discarded_introns: set[Tuple[int, int]] = set()

    def collect_introns(self, read_assignments) -> Dict[Tuple[int, int], int]:
        all_introns: Dict[Tuple[int, int], int] = defaultdict(int)
        for assignment in read_assignments:
            if not assignment.corrected_introns or assignment.multimapper:
                continue
            for intron in assignment.corrected_introns:
                all_introns[intron] += 1
        return all_introns

    def construct_similar_intron_map(self, all_introns: Dict[Tuple[int, int], int]):
        ordered = sorted(all_introns.keys())
        similar_map: Dict[Tuple[int, int], list] = defaultdict(list)
        for i, intron in enumerate(ordered):
            j = i + 1
            while j < len(ordered) and abs(ordered[j][0] - intron[0]) <= self.delta:
                if abs(ordered[j][1] - intron[1]) <= self.delta:
                    similar_map[intron].append(ordered[j])
                    similar_map[ordered[j]].append(intron)
                j += 1
        return similar_map

    def cluster_introns(self, all_introns: Dict[Tuple[int, int], int], min_count: int) -> None:
        similar_map = self.construct_similar_intron_map(all_introns)
        for count, intron in sorted(((v, k) for k, v in all_introns.items()), reverse=True):
            if intron in self.known_introns:
                self.clustered_introns[intron] = count
            elif intron in similar_map:
                similar = [
                    (count, si) for si in similar_map[intron] if si in self.clustered_introns
                ]
                if similar:
                    sub = sorted(similar, reverse=True)[0]
                    self.clustered_introns[sub[1]] += sub[0]
                    self.intron_correction_map[intron] = sub[1]
                else:
                    self.clustered_introns[intron] = count
            elif count < min_count:
                self.discarded_introns.add(intron)
            else:
                self.clustered_introns[intron] = count

    def process(self, read_assignments, min_count: int) -> None:
        all_introns = self.collect_introns(read_assignments)
        self.cluster_introns(all_introns, min_count)

    def add_substitute(self, original_intron: Tuple[int, int], substitute_intron: Tuple[int, int]) -> None:
        self.clustered_introns[substitute_intron] += self.clustered_introns[original_intron]
        del self.clustered_introns[original_intron]
        self.intron_correction_map[original_intron] = substitute_intron

    def discard(self, intron: Tuple[int, int]) -> None:
        self.discarded_introns.add(intron)
        if intron in self.clustered_introns:
            del self.clustered_introns[intron]

    def substitute(self, v: Tuple[int, int]) -> Tuple[int, int]:
        return self.intron_correction_map.get(v, v)

    def simplify_correction_map(self) -> None:
        all_introns = sorted(self.intron_correction_map.keys())
        to_remove = set()
        for intron in all_introns:
            subs = self.intron_correction_map[intron]
            if subs in self.discarded_introns:
                to_remove.add(intron)
                continue
            if subs not in self.intron_correction_map:
                continue
            while subs in self.intron_correction_map:
                subs = self.intron_correction_map[subs]
            if subs in self.discarded_introns:
                to_remove.add(intron)
                continue
            self.intron_correction_map[intron] = subs

        for intron in to_remove:
            self.discard(intron)
            del self.intron_correction_map[intron]


class IntronGraph:
    """Directed splice-junction graph with tip/bulge collapse and polyA/TSS terminal vertex attachment."""

    def __init__(
        self,
        params,
        gene_info,
        read_assignments,
        polya_predictions: Optional[List[int]] = None,
        tss_predictions: Optional[List[int]] = None,
    ) -> None:
        self.params = params
        self.gene_info = gene_info
        self.read_assignments = read_assignments
        self.polya_predictions = polya_predictions
        self.tss_predictions = tss_predictions

        self.incoming_edges: Dict[Tuple[int, int], set] = defaultdict(set)
        self.outgoing_edges: Dict[Tuple[int, int], set] = defaultdict(set)
        self.intron_collector = IntronCollector(gene_info, params.delta)
        self.max_coverage = 0
        self.edge_weights: Dict[Tuple[Tuple[int, int], Tuple[int, int]], int] = defaultdict(int)

        self.starting_known_positions: Dict[Tuple[int, int], List[int]] = defaultdict(list)
        self.terminal_known_positions: Dict[Tuple[int, int], List[int]] = defaultdict(list)
        for t, introns in self.gene_info.all_isoforms_introns.items():
            if not introns:
                continue
            self.starting_known_positions[introns[0]].append(self.gene_info.all_isoforms_exons[t][0][0])
            self.terminal_known_positions[introns[-1]].append(self.gene_info.all_isoforms_exons[t][-1][1])

        self.intron_collector.process(read_assignments, self.params.min_novel_intron_count)
        self.construct()
        if getattr(self.params, "debug", False):
            self.print_graph()
        self.simplify()
        self.attach_terminal_positions()
        if getattr(self.params, "debug", False):
            self.print_graph()

    def add_edge(self, v1: Tuple[int, int], v2: Tuple[int, int]) -> None:
        v1 = self.intron_collector.substitute(v1)
        v2 = self.intron_collector.substitute(v2)
        self.outgoing_edges[v1].add(v2)
        self.incoming_edges[v2].add(v1)
        self.edge_weights[(v1, v2)] += 1

    def is_isolated(self, v: Tuple[int, int]) -> bool:
        return not self.outgoing_edges[v] and not self.incoming_edges[v]

    def is_monointron(self, v: Tuple[int, int]) -> bool:
        no_out = not self.outgoing_edges[v] or all(is_terminal_vertex(x) for x in self.outgoing_edges[v])
        no_in = not self.incoming_edges[v] or all(is_starting_vertex(x) for x in self.incoming_edges[v])
        return no_in and no_out

    def signleton_dead_end(self, v: Tuple[int, int]) -> set:
        path = set()
        while self.intron_collector.clustered_introns[v] == 1:
            path.add(v)
            if len(self.outgoing_edges[v]) != 1:
                break
            v = next(iter(self.outgoing_edges[v]))
        return path if len(self.outgoing_edges[v]) == 0 else set()

    def signleton_dead_start(self, v: Tuple[int, int]) -> set:
        path = set()
        while self.intron_collector.clustered_introns[v] == 1:
            path.add(v)
            if len(self.incoming_edges[v]) != 1:
                break
            v = next(iter(self.incoming_edges[v]))
        return path if len(self.incoming_edges[v]) == 0 else set()

    def get_outgoing(self, intron: Tuple[int, int], v_type: Optional[int] = None) -> List[Tuple[int, int]]:
        if intron not in self.outgoing_edges:
            return []
        if v_type is None:
            return sorted(v for v in self.outgoing_edges[intron] if v[0] >= 0)
        return sorted(v for v in self.outgoing_edges[intron] if v[0] == v_type)

    def get_incoming(self, intron: Tuple[int, int], v_type: Optional[int] = None) -> List[Tuple[int, int]]:
        if intron not in self.incoming_edges:
            return []
        if v_type is None:
            return sorted(v for v in self.incoming_edges[intron] if v[0] >= 0)
        return sorted(v for v in self.incoming_edges[intron] if v[0] == v_type)

    def collapse_vertex(self, to_collapse: Tuple[int, int], substitute_vertex: Tuple[int, int]) -> None:
        self.outgoing_edges[substitute_vertex].update(self.outgoing_edges[to_collapse])
        for i in self.outgoing_edges[to_collapse]:
            self.incoming_edges[i].remove(to_collapse)
            self.incoming_edges[i].add(substitute_vertex)
        self.incoming_edges[substitute_vertex].update(self.incoming_edges[to_collapse])
        for i in self.incoming_edges[to_collapse]:
            self.outgoing_edges[i].remove(to_collapse)
            self.outgoing_edges[i].add(substitute_vertex)
        self.intron_collector.add_substitute(to_collapse, substitute_vertex)

    def construct(self) -> None:
        for assignment in self.read_assignments:
            if assignment.multimapper or any(
                intron in self.intron_collector.discarded_introns for intron in assignment.corrected_introns
            ):
                continue
            introns = assignment.corrected_introns
            for i in range(len(introns) - 1):
                self.add_edge(introns[i], introns[i + 1])
        self.max_coverage = (
            max(self.intron_collector.clustered_introns.values())
            if self.intron_collector.clustered_introns
            else 0
        )

    def simplify(self) -> None:
        self.clean_tips_and_bulges()
        self.remove_singleton_dead_ends()
        self.remove_isolates()
        self.intron_collector.simplify_correction_map()

    def clean_tips_and_bulges(self) -> None:
        to_remove = set()
        for current in sorted(self.outgoing_edges.keys()):
            sub = self.collapse_vertex_set(self.outgoing_edges[current])
            for i in sorted(sub.keys()):
                if i not in to_remove:
                    to_remove.add(i)
                    self.collapse_vertex(i, sub[i])
        for i in to_remove:
            del self.outgoing_edges[i]
            del self.incoming_edges[i]

        to_remove.clear()
        for current in sorted(self.incoming_edges.keys()):
            sub = self.collapse_vertex_set(self.incoming_edges[current])
            for i in sorted(sub.keys()):
                if i not in to_remove:
                    to_remove.add(i)
                    self.collapse_vertex(i, sub[i])
        for i in to_remove:
            del self.outgoing_edges[i]
            del self.incoming_edges[i]

    def remove_singleton_dead_ends(self) -> None:
        to_clean: Dict[Tuple[int, int], set] = {}
        for current in sorted(self.outgoing_edges.keys()):
            if self.intron_collector.clustered_introns[current] < self.params.singleton_adjacent_cov:
                continue
            out_paths = [self.signleton_dead_end(i) for i in self.outgoing_edges[current]]
            if any(len(p) == 0 for p in out_paths):
                continue
            to_clean[current] = {v for p in out_paths for v in p}

        for intron, verts in to_clean.items():
            self.outgoing_edges[intron] = set()
            for i in verts:
                self.outgoing_edges.pop(i, None)
                self.incoming_edges.pop(i, None)

        to_clean.clear()
        for current in sorted(self.incoming_edges.keys()):
            if self.intron_collector.clustered_introns[current] < self.params.singleton_adjacent_cov:
                continue
            in_paths = [self.signleton_dead_start(i) for i in self.incoming_edges[current]]
            if any(len(p) == 0 for p in in_paths):
                continue
            to_clean[current] = {v for p in in_paths for v in p}

        for intron, verts in to_clean.items():
            self.incoming_edges[intron] = set()
            for i in verts:
                self.outgoing_edges.pop(i, None)
                self.incoming_edges.pop(i, None)

    def remove_isolates(self) -> None:
        to_remove = set()
        isolated = {intron for intron in self.intron_collector.clustered_introns if self.is_isolated(intron)}
        sub = self.collapse_vertex_set(isolated)
        for i in sorted(sub.keys()):
            to_remove.add(i)
            self.collapse_vertex(i, sub[i])
        for i in to_remove:
            del self.outgoing_edges[i]
            del self.incoming_edges[i]

        to_remove.clear()
        cutoff = self.params.min_novel_isolated_intron_abs
        for intron in isolated:
            if intron not in self.intron_collector.clustered_introns or intron in self.intron_collector.known_introns:
                continue
            if self.intron_collector.clustered_introns[intron] < cutoff:
                self.intron_collector.discard(intron)
                to_remove.add(intron)
        for i in to_remove:
            del self.outgoing_edges[i]
            del self.incoming_edges[i]

    def collapse_vertex_set(self, vertex_set) -> Dict[Tuple[int, int], Tuple[int, int]]:
        if len(vertex_set) <= 1:
            return {}
        approved = set()
        substitute_dict: Dict[Tuple[int, int], Tuple[int, int]] = {}
        for count, vertex in sorted(((self.intron_collector.clustered_introns[i], i) for i in vertex_set), reverse=True):
            similar = []
            for i in sorted(approved):
                sd = abs(i[0] - vertex[0])
                ed = abs(i[1] - vertex[1])
                if (
                    sd < self.params.graph_clustering_distance
                    and ed < self.params.graph_clustering_distance
                    and count < self.intron_collector.clustered_introns[i] * self.params.graph_clustering_ratio
                ):
                    similar.append((sd + ed, i))
            if not similar:
                approved.add(vertex)
            else:
                substitute_dict[vertex] = sorted(similar)[0][1]
        return substitute_dict

    def attach_terminal_positions(self) -> None:
        polya_ends, read_ends, polyt_starts, read_starts = self.collect_terminal_positions()
        introns = sorted(self.intron_collector.clustered_introns)
        self._attach_side(introns, polya_ends, read_ends, read_end=True)
        self._attach_side(introns, polyt_starts, read_starts, read_end=False)

    def _attach_side(self, introns: List, polya_confirmed_positions: dict, read_terminal_positions: dict, read_end: bool) -> None:
        polya_positions = {
            intron: self._refine_positions(
                self.cluster_polya_positions(polya_confirmed_positions[intron], intron, read_end),
                self.polya_predictions,
            )
            for intron in introns
        }

        tss_positions = {}
        remaining_terminal_positions = {}
        for intron in introns:
            confirmed: Dict[int, int] = {}
            remaining: Dict[int, int] = {}
            for pos, cnt in read_terminal_positions[intron].items():
                nearest, diff = find_closest(pos, self.tss_predictions) if self.tss_predictions else (None, None)
                if nearest is not None and diff <= self.params.apa_delta:
                    confirmed[nearest] = confirmed.get(nearest, 0) + cnt
                else:
                    remaining[pos] = cnt
            tss_positions[intron] = confirmed
            remaining_terminal_positions[intron] = remaining

        cutoffs = {}
        extra_positions = {}
        for intron in introns:
            clustered_polyas = polya_positions[intron]
            cutoff = self.params.terminal_position_abs
            if clustered_polyas:
                cutoff = max(cutoff, max(clustered_polyas.values()) * self.params.terminal_position_rel)
                extra: Dict[int, int] = {}
                furthest = max(clustered_polyas.keys()) if read_end else min(clustered_polyas.keys())
                for pos, cnt in remaining_terminal_positions[intron].items():
                    if (read_end and pos >= furthest + self.params.apa_delta) or (
                        not read_end and pos <= furthest - self.params.apa_delta
                    ):
                        extra[pos] = cnt
            else:
                extra = remaining_terminal_positions[intron]

            neighbors = self.outgoing_edges if read_end else self.incoming_edges
            if intron in neighbors and neighbors[intron]:
                neighbor_cov = max(self.intron_collector.clustered_introns[i] for i in neighbors[intron])
                cutoff = max(cutoff, neighbor_cov * self.params.terminal_internal_position_rel)
            cutoffs[intron] = cutoff
            extra_positions[intron] = extra

        terminal_predictions = self.polya_predictions if read_end else self.tss_predictions
        terminal_positions = {
            intron: self._refine_positions(
                self.cluster_terminal_positions(extra_positions[intron], read_end=read_end, cutoff=cutoffs[intron]),
                terminal_predictions,
            )
            for intron in introns
        }

        polya_v = TerminalVertex.polya if read_end else TerminalVertex.polyt
        read_v = TerminalVertex.read_end if read_end else TerminalVertex.read_start
        tss_v = TerminalVertex.tss_right if read_end else TerminalVertex.tss_left
        edges = self.outgoing_edges if read_end else self.incoming_edges
        for intron in introns:
            for pos in polya_positions[intron]:
                edges[intron].add((polya_v, pos))
            for pos in tss_positions[intron]:
                edges[intron].add((tss_v, pos))
            for pos in terminal_positions[intron]:
                edges[intron].add((read_v, pos))

    def _refine_positions(
        self, clustered: Dict[int, int], predicted_positions: Optional[List[int]]
    ) -> Dict[int, int]:
        if not clustered or not predicted_positions:
            return clustered
        refined: Dict[int, int] = {}
        for pos, cnt in clustered.items():
            nearest, diff = find_closest(pos, predicted_positions)
            new_pos = nearest if (nearest is not None and diff <= self.params.apa_delta) else pos
            refined[new_pos] = refined.get(new_pos, 0) + cnt
        return refined

    def collect_terminal_positions(self):
        polya_ends = defaultdict(lambda: defaultdict(int))
        read_ends = defaultdict(lambda: defaultdict(int))
        polyt_starts = defaultdict(lambda: defaultdict(int))
        read_starts = defaultdict(lambda: defaultdict(int))

        for assignment in self.read_assignments:
            if assignment.multimapper or not assignment.corrected_introns:
                continue
            if any(i in self.intron_collector.discarded_introns for i in assignment.corrected_introns):
                continue

            starting_intron = self.intron_collector.substitute(assignment.corrected_introns[0])
            read_start = assignment.corrected_exons[0][0]
            if read_start < starting_intron[0]:
                polyt_det = assignment.strand == "-" and (
                    assignment.polya_info.external_polyt_pos != -1 or assignment.polya_info.internal_polyt_pos != -1
                )
                if polyt_det:
                    polyt_starts[starting_intron][read_start] += 1
                elif not self.is_start_internal(starting_intron, read_start):
                    read_starts[starting_intron][read_start] += 1

            term_intron = self.intron_collector.substitute(assignment.corrected_introns[-1])
            read_end = assignment.corrected_exons[-1][1]
            if read_end > term_intron[1]:
                polya_det = assignment.strand == "+" and (
                    assignment.polya_info.external_polya_pos != -1 or assignment.polya_info.internal_polya_pos != -1
                )
                if polya_det:
                    polya_ends[term_intron][read_end] += 1
                elif not self.is_end_internal(term_intron, read_end):
                    read_ends[term_intron][read_end] += 1

        return polya_ends, read_ends, polyt_starts, read_starts

    def cluster_polya_positions(self, position_dict, intron: Tuple[int, int], read_end: bool) -> Dict[int, int]:
        clustered_counts: Dict[int, int] = {}
        if not position_dict:
            return clustered_counts

        known_positions = self.terminal_known_positions[intron] if read_end else self.starting_known_positions[intron]
        while position_dict:
            top_pos = max(position_dict.items(), key=lambda x: x[1])[0]
            nearest, diff = find_closest(top_pos, known_positions)
            if nearest and diff <= self.params.apa_delta:
                top_pos = nearest
            if read_end:
                assert top_pos > intron[1]
            else:
                assert top_pos < intron[0]
            total_cnt = 0
            for pos in range(top_pos - self.params.apa_delta, top_pos + self.params.apa_delta + 1):
                if pos in position_dict:
                    total_cnt += position_dict.pop(pos)
            clustered_counts[top_pos] = total_cnt

        max_count = max(clustered_counts.values())
        if max_count == self.params.terminal_position_abs:
            return clustered_counts
        cutoff = max(max_count * self.params.terminal_position_rel, self.params.terminal_position_abs)
        return {k: v for k, v in clustered_counts.items() if v >= cutoff}

    def cluster_terminal_positions(self, position_dict, read_end: bool, cutoff: float = 0) -> Dict[int, int]:
        if not position_dict:
            return {}
        total_cnt = sum(position_dict.values())
        if total_cnt < cutoff:
            return {}
        pos = max(position_dict.keys()) if read_end else min(position_dict.keys())
        return {pos: total_cnt}

    def is_start_internal(self, intron: Tuple[int, int], read_start: int) -> bool:
        return any(inc[1] - self.params.delta <= read_start for inc in self.incoming_edges[intron])

    def is_end_internal(self, intron: Tuple[int, int], read_end: int) -> bool:
        return any(out[0] + self.params.delta >= read_end for out in self.outgoing_edges[intron])

    def get_connected_component(self, path) -> set:
        q = deque()
        processed = set()
        for intron in path:
            if intron[0] < 0:
                continue
            q.append(intron)
            processed.add(intron)
        while q:
            intron = q.popleft()
            for i in self.outgoing_edges[intron]:
                if i[0] >= 0 and i not in processed:
                    processed.add(i)
                    q.append(i)
            for i in self.incoming_edges[intron]:
                if i[0] >= 0 and i not in processed:
                    processed.add(i)
                    q.append(i)
        return processed

    def get_max_component_coverage(self, path) -> int:
        comp = self.get_connected_component(path)
        if not comp:
            return 0
        return max(self.intron_collector.clustered_introns[i] for i in comp)

    def get_overlapping_component_max_coverage(self, coordinates: Tuple[int, int]) -> int:
        processed = set()
        for intron in set(self.outgoing_edges.keys()) | set(self.incoming_edges.keys()):
            if intron in processed:
                continue
            max_r = intron[1]
            for v in self.outgoing_edges[intron]:
                max_r = max(max_r, v[1] if is_terminal_vertex(v) else v[0])
            min_l = intron[0]
            for v in self.incoming_edges[intron]:
                min_l = min(min_l, v[1])
            if overlaps(coordinates, (min_l, max_r)):
                processed.update(self.get_connected_component([intron]))
        if not processed:
            return 0
        return max(self.intron_collector.clustered_introns[i] for i in processed)

    def print_graph(self) -> None:
        logger.debug("Vertices: %d, substituted: %d", len(self.intron_collector.clustered_introns), len(self.intron_collector.intron_correction_map))


class IntronPathProcessor:
    def __init__(self, params, intron_graph: IntronGraph) -> None:
        self.params = params
        self.intron_graph = intron_graph
        self.all_vertices = set(self.intron_graph.intron_collector.clustered_introns.keys())
        for edge_set in self.intron_graph.outgoing_edges.values():
            self.all_vertices.update(edge_set)
        for edge_set in self.intron_graph.incoming_edges.values():
            self.all_vertices.update(edge_set)
        self.visited = set()

    def visit_vertex(self, v) -> None:
        if v in self.all_vertices:
            self.visited.add(v)

    def visit_path(self, p) -> None:
        for v in p:
            self.visit_vertex(v)

    def thread_introns(self, introns):
        path = []
        for intron in introns:
            if intron in self.intron_graph.intron_collector.discarded_introns:
                return None
            path.append(self.intron_graph.intron_collector.substitute(intron))
        return path

    def thread_ends(self, intron, end: int, trusted: bool = False):
        possible_polyas = self.intron_graph.get_outgoing(intron, TerminalVertex.polya)
        if trusted:
            for v in possible_polyas:
                if abs(v[1] - end) <= self.params.apa_delta:
                    return v

        outgoing_introns = self.intron_graph.get_outgoing(intron)
        if outgoing_introns:
            rightmost_exon_end = max(i[0] for i in outgoing_introns) - 1
            if not trusted and end <= rightmost_exon_end + self.params.delta:
                return None

        all_ends = sorted(
            self.intron_graph.get_outgoing(intron, TerminalVertex.read_end)
            + self.intron_graph.get_outgoing(intron, TerminalVertex.tss_right)
            + possible_polyas,
            key=lambda x: x[1],
        )
        if not all_ends:
            return None
        rightmost_end = all_ends[-1]
        if trusted and end >= rightmost_end[1] and rightmost_end[0] == TerminalVertex.read_end:
            return rightmost_end
        if not trusted and end <= rightmost_end[1] + self.params.apa_delta and (
            len(all_ends) <= 1 or end > all_ends[-2][1]
        ):
            return rightmost_end
        return None

    def thread_starts(self, intron, start: int, trusted: bool = False):
        possible_polyas = self.intron_graph.get_incoming(intron, TerminalVertex.polyt)
        if trusted:
            for v in possible_polyas:
                if abs(v[1] - start) <= self.params.apa_delta:
                    return v

        incoming_introns = self.intron_graph.get_incoming(intron)
        if incoming_introns:
            leftmost_exon_start = min(i[1] for i in incoming_introns) + 1
            if not trusted and start >= leftmost_exon_start - self.params.delta:
                return None

        all_starts = sorted(
            self.intron_graph.get_incoming(intron, TerminalVertex.read_start)
            + self.intron_graph.get_incoming(intron, TerminalVertex.tss_left)
            + possible_polyas,
            key=lambda x: x[1],
        )
        if not all_starts:
            return None
        leftmost_start = all_starts[0]
        if trusted and start <= leftmost_start[1] and leftmost_start[0] == TerminalVertex.read_start:
            return leftmost_start
        if not trusted and start >= leftmost_start[1] and (len(all_starts) <= 1 or start < all_starts[1][1]):
            return leftmost_start
        return None


class IntronPathStorage:
    def __init__(self, params, path_processor: IntronPathProcessor) -> None:
        self.params = params
        self.path_processor = path_processor
        self.intron_graph = path_processor.intron_graph
        self.paths: Dict[tuple, int] = defaultdict(int)
        self.fl_paths: set[tuple] = set()
        self.paths_to_reads: Dict[tuple, list] = defaultdict(list)

    def fill(self, read_assignments) -> None:
        for a in read_assignments:
            if a.multimapper:
                continue
            intron_path = self.path_processor.thread_introns(a.corrected_introns)
            if not intron_path:
                continue
            read_end = a.corrected_exons[-1][1]
            is_end_trusted = a.strand == "+" and (
                a.polya_info.external_polya_pos != -1 or a.polya_info.internal_polya_pos != -1
            )
            term_v = self.path_processor.thread_ends(intron_path[-1], read_end, is_end_trusted)
            if term_v:
                intron_path.append(term_v)

            read_start = a.corrected_exons[0][0]
            is_start_trusted = a.strand == "-" and (
                a.polya_info.external_polyt_pos != -1 or a.polya_info.internal_polyt_pos != -1
            )
            start_v = self.path_processor.thread_starts(intron_path[0], read_start, is_start_trusted)
            if start_v:
                intron_path = [start_v] + intron_path

            path_tuple = tuple(intron_path)
            self.paths[path_tuple] += 1
            if term_v and start_v:
                if not self.params.requires_polya_for_construction or (
                    term_v[0] == TerminalVertex.polya or start_v[0] == TerminalVertex.polyt
                ):
                    self.fl_paths.add(path_tuple)
            self.paths_to_reads[path_tuple].append(a)
