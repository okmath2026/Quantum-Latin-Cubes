# -*- coding: utf-8 -*-
"""
Incremental low-cardinality probe for QLS(5).

This is a parameterized and resumable version of the c=10 deterministic probe.
It is designed for the c=11 frontier, where the old generic local-support sieve
does not finish quickly.

Important implementation detail:
    New labels are introduced in canonical first-occurrence order.  This avoids
    enumerating all permutations of the new labels and then canonicalizing only
    at the leaves.

Typical commands:
    python qls5_incremental_probe.py --t 6 --seconds 600
    python qls5_incremental_probe.py --t 6 --start 41 --end 80
    python qls5_incremental_probe.py --t 5 --graph-source local_support_ok
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
OLD = tuple(range(N))
OUTDIR = os.path.dirname(os.path.abspath(__file__))


def graph_sort_key(item):
    return (
        item["active_edges"],
        tuple(item["row_degrees"]),
        tuple(item["col_degrees"]),
        item["graph_key"],
    )


def enumerate_graph(mat, t, deadline=None, max_patterns=None):
    new_labels = tuple(range(N, N + t))
    grid = [[None] * N for _ in range(N)]
    grid[0] = list(OLD)
    row_old = [set() for _ in range(N)]
    col_old = [{c} for c in OLD]
    row_new = [set() for _ in range(N)]
    col_new = [set() for _ in range(N)]
    row_all = [set(grid[0])] + [set() for _ in range(1, N)]
    col_all = [{c} for c in OLD]
    zero_sets = {k: set() for k in new_labels}
    new_edges = set()

    cells = [(r, c) for r in range(1, N) for c in range(N)]
    row_deg = {r: sum(mat[r - 1][c] for c in range(N)) for r in range(1, N)}
    col_deg = {c: sum(mat[r - 1][c] for r in range(1, N)) for c in range(N)}
    cells.sort(
        key=lambda rc: (mat[rc[0] - 1][rc[1]], row_deg[rc[0]] + col_deg[rc[1]]),
        reverse=True,
    )
    active_suffix = [0] * (len(cells) + 1)
    for i in range(len(cells) - 1, -1, -1):
        r, c = cells[i]
        active_suffix[i] = active_suffix[i + 1] + (1 if mat[r - 1][c] == 1 else 0)

    support_type_counter = Counter()
    pattern_count = 0
    timed_out = False
    max_label_used = N - 1

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
        supports = tuple(tuple(i for i in OLD if i not in zero_sets[k]) for k in new_labels)
        return repr((supports, tuple(sorted(new_edges))))

    def active_label_options(current_max):
        # Canonical first-occurrence rule: already introduced labels may be
        # reused; the next unseen label may be introduced.  Later unseen labels
        # cannot appear before the next one.
        upper = min(current_max + 1, N + t - 1)
        return range(N, upper + 1)

    def bt(i, current_max):
        nonlocal timed_out, pattern_count
        if deadline is not None and time.time() > deadline:
            timed_out = True
            return
        if max_patterns is not None and pattern_count >= max_patterns:
            timed_out = True
            return

        introduced = max(0, current_max - N + 1)
        if introduced + active_suffix[i] < t:
            return

        if i == len(cells):
            if introduced != t:
                return
            pattern_count += 1
            support_type_counter[support_edge_key()] += 1
            return

        r, c = cells[i]
        active = mat[r - 1][c] == 1
        labels = active_label_options(current_max) if active else OLD
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
                    bt(i + 1, max(current_max, label))
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
                    bt(i + 1, current_max)
                    col_old[c].remove(label)
                    row_old[r].remove(label)
                restore_zero(old_zero)
            col_all[c].remove(label)
            row_all[r].remove(label)
            grid[r][c] = None
            if timed_out:
                return

    bt(0, max_label_used)
    return pattern_count, support_type_counter, timed_out


def classify_records(global_support_types, t):
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
    return reason_type_counts, reason_pattern_counts, records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--t", type=int, required=True, help="number of new labels; cardinality is 5+t")
    parser.add_argument("--graph-source", default="old_and_new_ok", choices=["old_and_new_ok", "local_support_ok"])
    parser.add_argument("--start", type=int, default=1, help="1-based graph index after sorting")
    parser.add_argument("--end", type=int, default=None, help="inclusive 1-based graph index after sorting")
    parser.add_argument("--seconds", type=float, default=None)
    parser.add_argument("--max-patterns-per-graph", type=int, default=None)
    parser.add_argument("--out-prefix", default=None)
    parser.add_argument("--per-graph-only", action="store_true", help="skip global classification/record output")
    args = parser.parse_args()

    c = N + args.t
    graph_path = os.path.join(OUTDIR, f"qls5_c{c}_active_graph_general.json")
    graph_data = json.load(open(graph_path, encoding="utf-8"))
    source = graph_data[args.graph_source]
    if source is None:
        raise SystemExit(f"Graph source {args.graph_source!r} is null in {graph_path}")
    graphs = sorted(source["graph_types"], key=graph_sort_key)

    start = max(1, args.start)
    end = args.end if args.end is not None else len(graphs)
    selected = list(enumerate(graphs, 1))[start - 1:end]
    deadline = time.time() + args.seconds if args.seconds is not None else None
    out_prefix = args.out_prefix or f"qls5_c{c}_incremental_{args.graph_source}"

    global_support_types = Counter()
    graph_summaries = []
    stopped_early = False

    for graph_index, item in selected:
        mat = tuple(tuple(row) for row in item["matrix"])
        pattern_count, support_types, timed_out = enumerate_graph(
            mat,
            args.t,
            deadline=deadline,
            max_patterns=args.max_patterns_per_graph,
        )
        global_support_types.update(support_types)
        summary = {
            "graph_index": graph_index,
            "active_edges": item["active_edges"],
            "row_degrees": item["row_degrees"],
            "col_degrees": item["col_degrees"],
            "label_patterns": pattern_count,
            "support_edge_types": len(support_types),
            "incomplete": timed_out,
        }
        graph_summaries.append(summary)

        graph_out = {
            "n": N,
            "t": args.t,
            "cardinality": c,
            "graph_source": args.graph_source,
            "summary": summary,
            "support_edge_types": [
                {"type": k, "count": v} for k, v in sorted(support_types.items())
            ],
        }
        graph_path_out = os.path.join(OUTDIR, f"{out_prefix}_graph_{graph_index:03d}.json")
        with open(graph_path_out, "w", encoding="utf-8") as f:
            json.dump(graph_out, f, ensure_ascii=False, indent=2)

        print(
            "graph", graph_index,
            "edges", item["active_edges"],
            "patterns", pattern_count,
            "support_types", len(support_types),
            "incomplete", timed_out,
        )
        if timed_out:
            stopped_early = True
            break

    if args.per_graph_only:
        out = {
            "complete": (
                not stopped_early
                and len(graph_summaries) == len(selected)
                and start == 1
                and end == len(graphs)
            ),
            "n": N,
            "t": args.t,
            "cardinality": c,
            "graph_source": args.graph_source,
            "available_graph_types": len(graphs),
            "range_start": start,
            "range_end": end,
            "processed_graph_types": len(graph_summaries),
            "graph_summaries": graph_summaries,
            "label_pattern_count": sum(s["label_patterns"] for s in graph_summaries),
            "support_edge_type_count": None,
            "unclassified_count": None,
            "reason_type_counts": None,
            "reason_label_pattern_counts": None,
            "records": [],
        }
        suffix = f"{start:03d}_{start + len(graph_summaries) - 1:03d}_per_graph_only"
        out_path = os.path.join(OUTDIR, f"{out_prefix}_{suffix}.json")
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(out, f, ensure_ascii=False, indent=2)
        print("complete", out["complete"])
        print("processed_graph_types", out["processed_graph_types"])
        print("label_pattern_count", out["label_pattern_count"])
        print("support_edge_type_count", out["support_edge_type_count"])
        print("unclassified_count", out["unclassified_count"])
        print("wrote", out_path)
        return

    reason_type_counts, reason_pattern_counts, records = classify_records(global_support_types, args.t)
    complete_range = (
        not stopped_early
        and len(graph_summaries) == len(selected)
        and start == 1
        and end == len(graphs)
    )
    out = {
        "complete": complete_range,
        "n": N,
        "t": args.t,
        "cardinality": c,
        "graph_source": args.graph_source,
        "available_graph_types": len(graphs),
        "range_start": start,
        "range_end": end,
        "processed_graph_types": len(graph_summaries),
        "graph_summaries": graph_summaries,
        "label_pattern_count": sum(global_support_types.values()),
        "support_edge_type_count": len(global_support_types),
        "unclassified_count": reason_type_counts.get("unclassified", 0),
        "reason_type_counts": dict(sorted(reason_type_counts.items())),
        "reason_label_pattern_counts": dict(sorted(reason_pattern_counts.items())),
        "records": records,
    }
    suffix = "complete" if complete_range else f"{start:03d}_{start + len(graph_summaries) - 1:03d}_partial"
    out_path = os.path.join(OUTDIR, f"{out_prefix}_{suffix}.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print("complete", out["complete"])
    print("processed_graph_types", out["processed_graph_types"])
    print("label_pattern_count", out["label_pattern_count"])
    print("support_edge_type_count", out["support_edge_type_count"])
    print("unclassified_count", out["unclassified_count"])
    print("reason_type_counts", out["reason_type_counts"])
    print("reason_label_pattern_counts", out["reason_label_pattern_counts"])
    print("wrote", out_path)


if __name__ == "__main__":
    main()
