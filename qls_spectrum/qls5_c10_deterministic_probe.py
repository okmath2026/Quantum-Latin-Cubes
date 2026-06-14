# -*- coding: utf-8 -*-
"""
Deterministic support-edge probe for QLS(5) with cardinality 10.

This is the c=10 analogue of qls5_c9_deterministic_certificate.py.  It reads
the 40 active graph types surviving the local support sieve in
qls5_c10_active_graph_general.json, enumerates normalized old/new label
fillings, collapses them to support-edge types, and classifies those types by
the elementary obstruction rules in qls_support_obstruction_general.py.

If unclassified_count is zero, the script is a nonexistence certificate.
If not, the unclassified records identify the next algebraic cases to study.
"""

from __future__ import annotations

from collections import Counter
import argparse
import ast
import json
import os
import time

from qls_support_obstruction_general import classify


N = 5
T = 5
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


def enumerate_graph(mat, deadline=None, max_patterns=None):
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
    cells.sort(
        key=lambda rc: (mat[rc[0] - 1][rc[1]], row_deg[rc[0]] + col_deg[rc[1]]),
        reverse=True,
    )

    seen_patterns = set()
    support_type_counter = Counter()
    timed_out = False

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
        nonlocal timed_out
        if deadline is not None and time.time() > deadline:
            timed_out = True
            return
        if max_patterns is not None and len(seen_patterns) >= max_patterns:
            timed_out = True
            return
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
            if timed_out:
                return

    bt(0)
    return seen_patterns, support_type_counter, timed_out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=None)
    parser.add_argument("--max-graphs", type=int, default=None)
    parser.add_argument("--max-patterns-per-graph", type=int, default=None)
    args = parser.parse_args()

    deadline = time.time() + args.seconds if args.seconds is not None else None
    graph_data = json.load(open(os.path.join(OUTDIR, "qls5_c10_active_graph_general.json"), encoding="utf-8"))
    graphs = graph_data["local_support_ok"]["graph_types"]
    if args.max_graphs is not None:
        graphs = graphs[:args.max_graphs]

    global_support_types = Counter()
    graph_summaries = []
    stopped_early = False
    for index, item in enumerate(graphs, 1):
        mat = tuple(tuple(row) for row in item["matrix"])
        patterns, support_types, timed_out = enumerate_graph(
            mat,
            deadline=deadline,
            max_patterns=args.max_patterns_per_graph,
        )
        global_support_types.update(support_types)
        graph_summaries.append({
            "graph_index": index,
            "active_edges": item["active_edges"],
            "row_degrees": item["row_degrees"],
            "col_degrees": item["col_degrees"],
            "label_patterns": len(patterns),
            "support_edge_types": len(support_types),
            "incomplete": timed_out,
        })
        print(
            "graph", index,
            "edges", item["active_edges"],
            "patterns", len(patterns),
            "support_types", len(support_types),
            "incomplete", timed_out,
        )
        if timed_out:
            stopped_early = True
            break

    reason_type_counts = Counter()
    reason_pattern_counts = Counter()
    records = []
    for type_key, count in sorted(global_support_types.items()):
        supports, edges = ast.literal_eval(type_key)
        result = classify(supports, edges, first_new_label=N)
        reason_type_counts[result["reason"]] += 1
        reason_pattern_counts[result["reason"]] += count
        records.append({
            "type": type_key,
            "pattern_count": count,
            "support_sizes": [len(s) for s in supports],
            "edge_count": len(edges),
            "reason": result["reason"],
            "witness": result["witness"],
        })

    out = {
        "complete": not stopped_early and len(graph_summaries) == len(graph_data["local_support_ok"]["graph_types"]),
        "active_graph_types_available": len(graph_data["local_support_ok"]["graph_types"]),
        "active_graph_types_processed": len(graph_summaries),
        "graph_summaries": graph_summaries,
        "label_pattern_count": sum(global_support_types.values()),
        "support_edge_type_count": len(global_support_types),
        "unclassified_count": reason_type_counts.get("unclassified", 0),
        "reason_type_counts": dict(sorted(reason_type_counts.items())),
        "reason_label_pattern_counts": dict(sorted(reason_pattern_counts.items())),
        "records": records,
    }
    suffix = "partial" if not out["complete"] else "complete"
    path = os.path.join(OUTDIR, f"qls5_c10_deterministic_probe_{suffix}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print("complete", out["complete"])
    print("active_graph_types_processed", out["active_graph_types_processed"])
    print("label_pattern_count", out["label_pattern_count"])
    print("support_edge_type_count", out["support_edge_type_count"])
    print("unclassified_count", out["unclassified_count"])
    print("reason_type_counts", out["reason_type_counts"])
    print("reason_label_pattern_counts", out["reason_label_pattern_counts"])
    print("wrote", path)


if __name__ == "__main__":
    main()
