# -*- coding: utf-8 -*-
"""
General support-edge obstruction classifier for low-cardinality QLS(5).

Input data:
    supports = tuple of coordinate-support tuples, one for each new label.
    edges    = tuple of pairs of new labels that must be orthogonal.

The classifier uses only elementary Hilbert-space facts.  It is deterministic
and independent of the label-pattern enumerators.
"""

from __future__ import annotations

from itertools import combinations


def cliques(vertices, edges):
    edge_set = {tuple(sorted(e)) for e in edges}
    for size in range(2, len(vertices) + 1):
        for clique in combinations(vertices, size):
            if all(tuple(sorted((a, b))) in edge_set for a, b in combinations(clique, 2)):
                yield clique


def classify(supports, edges, first_new_label=5):
    """Return a deterministic obstruction record, or ``unclassified``.

    The three rules are the ones used in the c=9 certificate:

    * clique_union_dim
    * two_2d_overlap_one_edge
    * same_2d_two_step
    * same_2d_basis_common_neighbor
    * same_3d_two_independent_common_neighbors
    """
    vertices = tuple(range(first_new_label, first_new_label + len(supports)))
    edge_set = {tuple(sorted(e)) for e in edges}
    support = {v: set(supports[v - first_new_label]) for v in vertices}

    reductions = []
    changed = True
    while changed:
        changed = False
        for u, v in edge_set:
            for a, b in ((u, v), (v, u)):
                Sa = support[a]
                Sb = support[b]
                inter = Sa & Sb
                if len(Sa) == 2 and len(inter) == 1:
                    coord = next(iter(inter))
                    if coord in Sb:
                        Sb.remove(coord)
                        reductions.append({
                            "two_dimensional_neighbor": a,
                            "reduced_vertex": b,
                            "removed_coordinate": coord,
                            "neighbor_support": sorted(Sa),
                        })
                        changed = True
                        if len(Sb) < 2:
                            return {
                                "reason": "forced_support_reduction_degeneracy",
                                "witness": {
                                    "reductions": reductions,
                                    "degenerate_vertex": b,
                                    "remaining_support": sorted(Sb),
                                },
                            }
        for u, w in edge_set:
            if support[u] != support[w] or len(support[u]) != 3:
                continue
            S = set(support[u])
            planes = []
            common_neighbors = []
            for v in vertices:
                if v in (u, w):
                    continue
                if tuple(sorted((u, v))) not in edge_set or tuple(sorted((w, v))) not in edge_set:
                    continue
                P = support[v] & S
                outside = support[v] - S
                if len(P) == 2 and len(outside) == 1:
                    planes.append(P)
                    common_neighbors.append(v)
            if len(planes) < 2:
                continue
            line = set.intersection(*planes)
            if len(line) != 1:
                continue
            coord = next(iter(line))
            for a in (u, w):
                if coord in support[a]:
                    support[a].remove(coord)
                    reductions.append({
                        "same_3d_edge": [u, w],
                        "reduced_vertex": a,
                        "removed_coordinate": coord,
                        "reason": "common_neighbor_planes_force_complement_line",
                    })
                    changed = True
                    if len(support[a]) < 2:
                        return {
                            "reason": "forced_support_reduction_degeneracy",
                            "witness": {
                                "reductions": reductions,
                                "degenerate_vertex": a,
                                "remaining_support": sorted(support[a]),
                            },
                        }
            for v, P in zip(common_neighbors, planes):
                for remove_coord in sorted(P - {coord}):
                    if remove_coord in support[v]:
                        support[v].remove(remove_coord)
                        reductions.append({
                            "same_3d_edge": [u, w],
                            "reduced_vertex": v,
                            "removed_coordinate": remove_coord,
                            "forced_line": coord,
                            "reason": "common_neighbor_projection_line",
                        })
                        changed = True
                        if len(support[v]) < 2:
                            return {
                                "reason": "forced_support_reduction_degeneracy",
                                "witness": {
                                    "reductions": reductions,
                                    "degenerate_vertex": v,
                                    "remaining_support": sorted(support[v]),
                                },
                            }

    for clique in cliques(vertices, edges):
        union_dim = len(set().union(*(support[v] for v in clique)))
        if union_dim < len(clique):
            return {
                "reason": "clique_union_dim",
                "witness": {
                    "clique": list(clique),
                    "union_dim": union_dim,
                    "clique_size": len(clique),
                    "reductions": reductions,
                },
            }

    for u, v in combinations(vertices, 2):
        if tuple(sorted((u, v))) not in edge_set:
            continue
        Su = support[u]
        Sv = support[v]
        if len(Su) == 2 and len(Sv) == 2 and len(Su & Sv) == 1:
            return {
                "reason": "two_2d_overlap_one_edge",
                "witness": {
                    "edge": [u, v],
                    "support_u": sorted(Su),
                    "support_v": sorted(Sv),
                    "intersection": sorted(Su & Sv),
                    "reductions": reductions,
                },
            }

    for u, w in combinations(vertices, 2):
        S = support[u]
        if len(S) != 2 or support[w] != S:
            continue
        for v in vertices:
            if v in (u, w) or support[v] != S:
                continue
            if tuple(sorted((u, v))) in edge_set and tuple(sorted((v, w))) in edge_set:
                return {
                    "reason": "same_2d_two_step",
                    "witness": {
                        "path": [u, v, w],
                        "support": sorted(S),
                        "reductions": reductions,
                    },
                }

    for u, w in combinations(vertices, 2):
        S = support[u]
        if len(S) != 3 or support[w] != S:
            continue
        common = [
            v for v in vertices
            if v not in (u, w)
            and tuple(sorted((u, v))) in edge_set
            and tuple(sorted((w, v))) in edge_set
        ]
        for a, b in combinations(common, 2):
            Sa = support[a]
            Sb = support[b]
            if (
                len(Sa) == 2
                and len(Sb) == 2
                and Sa <= S
                and Sb <= S
                and len(Sa & Sb) == 1
                and Sa | Sb == S
            ):
                return {
                    "reason": "same_3d_two_independent_common_neighbors",
                    "witness": {
                        "same_support_vertices": [u, w],
                        "same_support": sorted(S),
                        "common_neighbors": [a, b],
                        "neighbor_supports": [sorted(Sa), sorted(Sb)],
                        "reductions": reductions,
                    },
                }

    for u, v in combinations(vertices, 2):
        S = support[u]
        if len(S) != 2 or support[v] != S or tuple(sorted((u, v))) not in edge_set:
            continue
        for w in vertices:
            if w in (u, v):
                continue
            if tuple(sorted((u, w))) not in edge_set or tuple(sorted((v, w))) not in edge_set:
                continue
            outside = support[w] - S
            if len(outside) < 2:
                return {
                    "reason": "same_2d_basis_common_neighbor",
                    "witness": {
                        "basis_edge": [u, v],
                        "common_neighbor": w,
                        "basis_support": sorted(S),
                        "neighbor_support": sorted(support[w]),
                        "outside_support_size": len(outside),
                        "reductions": reductions,
                    },
                }

    return {
        "reason": "unclassified",
        "witness": {},
    }
