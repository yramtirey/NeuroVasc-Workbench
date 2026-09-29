import pandas as pd

from neurovasc.graph.branch_graph import BranchGraph


def branch_metrics_table(
    branch_graph: BranchGraph,
) -> pd.DataFrame:
    """
    Convert branch-level measurements into a tabular structure.
    """

    records = []

    for branch in branch_graph.branches:

        records.append(
            {
                "branch_id": branch.branch_id,
                "start_node": branch.start_node,
                "end_node": branch.end_node,
                "length_mm": branch.length_mm,
                "chord_length_mm": branch.chord_length_mm,
                "tortuosity": branch.tortuosity,
                "centerline_points": len(
                    branch.world_path
                ),
            }
        )

    return pd.DataFrame(
        records
    )