# -*- coding: utf-8 -*-
"""
Deterministic certificate for nonexistence of QLS(5) with cardinality 9.

This script is the proof-oriented version of qls5_c9_local_graph_search.py:
it does not use random vector witnesses and has no time limit.  It assumes the
active graph sieve has already produced qls5_c9_active_graph_general.json.

Pipeline:
1. Read the 14 active-cell graph types surviving the local support sieve.
2. Enumerate all normalized old/new label fillings for each graph.
3. Convert each filling to its support-edge type.
4. Classify the type by the deterministic obstructions in
   qls5_c9_support_obstruction.py.
"""

from __future__ import annotations

from collections import Counter
import ast
import json
import os

from qls5_c9_support_obstruction import classify


N = 5
T = 4
OLD = tuple(range(N))
NEW = tuple(range(N, N + T))
OUTDIR = os.path.dirname(os.path.abspath(__file__))


def canonicalize_new_labels(grid):
    mapping = {}
    nxt = N
    out = []
    for row in grid:
        new_row = []
        for x in row:
            if x >= N:
                if x not in mapping:
                    mapping[x] = nxt
                    nxt += 1
                x = mapping[x]
            new_row.append(x)
        out.append(tuple(new_row))
    return tuple(out)


def enumerate_graph(mat):
    grid = [[None] * N for _ in range(N)]
    grid[0] = list(OLD)
    row_old = [set() for _ in range(N)]
    col_old = [{c} for c in OLD]
    row_new = [set() for _ in range(N)]
    col_new = [set() for _ in range(N)]
    row_all = [set(grid[0])] + [set() for _ in range(1, N)]
    col_all = [{c} for c in OLD]
    zero_sets = {k: set() for k in NEW}
    new_edges = set()

    cells = [(r, c) for r in range(1, N) for c in range(N)]
    row_deg = {r: sum(mat[r - 1][c] for c in range(N)) for r in range(1, N)}
    col_deg = {c: sum(mat[r - 1][c] for r in range(1, N)) for c in range(N)}
    cells.sort(key=lambda rc: (mat[rc[0] - 1][rc[1]], row_deg[rc[0]] + col_deg[rc[1]]), reverse=True)

    seen_patterns = set()
    support_type_counter = Counter()

    def add_new_constraints_for_old(r, c, old_label):
        affected = row_new[r] | col_new[c]
        old_zero = {k: set(zero_sets[k]) for k in affected}
        ok = True
        for k in affected:
            zero_sets[k].add(old_label)
            if len(zero_sets[k]) > 3:
                ok = False
        return ok, old_zero

    def restore_zero(old_zero):
        for k, z in old_zero.items():
            zero_sets[k] = z

    def add_constraints_for_new(r, c, new_label):
        add = row_old[r] | col_old[c]
        if len(zero_sets[new_label] | add) > 3:
            return None
        old_z = set(zero_sets[new_label])
        old_edges = set(new_edges)
        zero_sets[new_label] |= add
        for k in row_new[r] | col_new[c]:
            if k != new_label:
                new_edges.add(tuple(sorted((k, new_label))))
        return old_z, old_edges

    def support_edge_key():
        supports = tuple(tuple(i for i in OLD if i not in zero_sets[k]) for k in NEW)
        return repr((supports, tuple(sorted(new_edges))))

    def bt(i):
        if i == len(cells):
            if not all(any(k in row_new[r] for r in range(1, N)) for k in NEW):
                return
            pattern = canonicalize_new_labels(grid)
            if pattern in seen_patterns:
                return
            seen_patterns.add(pattern)
            support_type_counter[support_edge_key()] += 1
            return

        r, c = cells[i]
        active = mat[r - 1][c] == 1
        labels = NEW if active else OLD
        for label in labels:
            if label in row_all[r] or label in col_all[c]:
                continue
            if not active and label == c:
                continue
            grid[r][c] = label
            row_all[r].add(label)
            col_all[c].add(label)
            if active:
                upd = add_constraints_for_new(r, c, label)
                if upd is not None:
                    old_z, old_edges = upd
                    row_new[r].add(label)
                    col_new[c].add(label)
                    bt(i + 1)
                    col_new[c].remove(label)
                    row_new[r].remove(label)
                    zero_sets[label] = old_z
                    new_edges.clear()
                    new_edges.update(old_edges)
            else:
                ok, old_zero = add_new_constraints_for_old(r, c, label)
                if ok:
                    row_old[r].add(label)
                    col_old[c].add(label)
                    bt(i + 1)
                    col_old[c].remove(label)
                    row_old[r].remove(label)
                restore_zero(old_zero)
            col_all[c].remove(label)
            row_all[r].remove(label)
            grid[r][c] = None

    bt(0)
    return seen_patterns, support_type_counter


def main():
    graph_data = json.load(open(os.path.join(OUTDIR, "qls5_c9_active_graph_general.json"), encoding="utf-8"))
    graphs = graph_data["local_support_ok"]["graph_types"]

    global_support_types = Counter()
    graph_summaries = []
    for index, item in enumerate(graphs, 1):
        mat = tuple(tuple(row) for row in item["matrix"])
        patterns, support_types = enumerate_graph(mat)
        global_support_types.update(support_types)
        graph_summaries.append({
            "graph_index": index,
            "active_edges": item["active_edges"],
            "row_degrees": item["row_degrees"],
            "col_degrees": item["col_degrees"],
            "label_patterns": len(patterns),
            "support_edge_types": len(support_types),
        })

    reason_type_counts = Counter()
    reason_pattern_counts = Counter()
    records = []
    for type_key, count in sorted(global_support_types.items()):
        supports, edges = ast.literal_eval(type_key)
        result = classify(supports, edges)
        reason_type_counts[result["reason"]] += 1
        reason_pattern_counts[result["reason"]] += count
        records.append({
            "type": type_key,
            "pattern_count": count,
            "reason": result["reason"],
            "witness": result["witness"],
        })

    out = {
        "active_graph_types": len(graphs),
        "graph_summaries": graph_summaries,
        "label_pattern_count": sum(global_support_types.values()),
        "support_edge_type_count": len(global_support_types),
        "unclassified_count": reason_type_counts.get("unclassified", 0),
        "reason_type_counts": dict(sorted(reason_type_counts.items())),
        "reason_label_pattern_counts": dict(sorted(reason_pattern_counts.items())),
        "records": records,
    }
    path = os.path.join(OUTDIR, "qls5_c9_deterministic_certificate.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print("active_graph_types", out["active_graph_types"])
    print("label_pattern_count", out["label_pattern_count"])
    print("support_edge_type_count", out["support_edge_type_count"])
    print("unclassified_count", out["unclassified_count"])
    print("reason_type_counts", out["reason_type_counts"])
    print("reason_label_pattern_counts", out["reason_label_pattern_counts"])
    print("wrote", path)


if __name__ == "__main__":
    main()
