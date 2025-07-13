from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
from typing import Any
import matplotlib.pyplot as plt
import numpy as np
import polars as pl


def get_slurm_df(file: Path) -> pl.DataFrame:
    rows = []

    with open(file, "r") as f:
        for line in f:
            fields = dict(re.findall(r"(\S+?)=(\S+)", line))
            rows.append(fields)

    df = pl.DataFrame(rows)
    fine_df = df.select(
        pl.col("JobId"),
        pl.col("JobState"),
        pl.col("StartTime"),
        pl.col("EndTime"),
        pl.col("Tres"),
    )

    return fine_df


@dataclass
class UsageTimeline:
    """
    x: Timestamps
    y: Number of processes being used at corresponding x.
       The last value of a correctly terminated job should be 0.
    """

    x: list[float]
    y: list[int]


class UsageTimelineProvider:
    def __init__(self, *args, **kwargs) -> None:
        self.timelines: list[UsageTimeline] = []

    def get_timelines(self) -> list[UsageTimeline]:
        return self.timelines

    def get_unified_timelines(self) -> tuple[list[float, list[list[int]]]]:
        xs = [i.x for i in self.get_timelines()]
        ys = [i.y for i in self.get_timelines()]
        shared_x = list(np.unique(np.concatenate(xs)))

        new_ys = [[] for _ in range(len(ys))]
        ys_idx = [0 for _ in range(len(ys))]

        for x in shared_x:
            for idx, x_list in enumerate(xs):
                if x in x_list:
                    orig_y_list = ys[idx]
                    current_orig_y_idx = ys_idx[idx]
                    new_ys[idx].append(orig_y_list[current_orig_y_idx])
                    # Pick updated orig_y value the next time
                    ys_idx[idx] += 1
                else:
                    new_ys[idx].append(None)

        for y in new_ys:
            self.fill_step_timeline(y)

        return shared_x, new_ys

    """
    Replace None values from a step timeline with appropriate values
    """

    def fill_step_timeline(self, lst: list[Any], boundary_val: Any = 0):
        # Propagate the last value to the left
        for i in range(len(lst) - 1, 0, -1):
            if lst[i] is None:
                lst[i] = boundary_val
            else:
                break

        # Propagate values to the right following a step post logic
        last_seen = boundary_val

        for i in range(len(lst)):
            if lst[i] is not None:
                last_seen = lst[i]
            lst[i] = last_seen

    def normalize_by_global_start(self):
        timelines = self.timelines
        global_start = min([min(t.x) for t in timelines])
        self.timelines = [
            UsageTimeline([x - global_start for x in t.x], t.y) for t in timelines
        ]


class DMRProvider(UsageTimelineProvider):
    def __init__(self, log_files: list[str | Path]):
        super().__init__()

        data: list[list[tuple[float, int]]] = []
        for path in log_files:
            # Load datasets
            cur_data = self.load_dmr_data(path)
            data.append(cur_data)

        # Flatten to find the earliest start time
        all_timestamps = [timestamp for timeline in data for timestamp, _ in timeline]
        if not all_timestamps:
            raise ValueError("No data found in either file.")

        self.timelines = [
            UsageTimeline(*[list(tup) for tup in zip(*execution)]) for execution in data
        ]

        self.normalize_by_global_start()
        self.fix_last_ys()

    def get_timelines(self):
        return self.timelines

    def load_dmr_data(self, filepath: str | Path) -> list[tuple[float, int]]:
        init_data: list[tuple[float, int]] = []
        last_fini = None

        with open(filepath, "r") as file:
            for line in file:
                stripped = line.strip()

                if stripped.startswith("DMR INIT"):
                    timestamp_match = re.search(r"Timestamp:\s*(\d+\.\d+)", stripped)
                    nodes_match = re.search(r"Nodes:\s*(\d+)", stripped)
                    if timestamp_match and nodes_match:
                        timestamp = float(timestamp_match.group(1))
                        nodes = int(nodes_match.group(1))
                        init_data.append((timestamp, nodes))

                elif stripped.startswith("DMR FINI"):
                    timestamp_match = re.search(r"Timestamp:\s*(\d+\.\d+)", stripped)
                    nodes_match = re.search(r"Nodes:\s*(\d+)", stripped)
                    if timestamp_match and nodes_match:
                        last_fini = (
                            float(timestamp_match.group(1)),
                            int(nodes_match.group(1)),
                        )

        if last_fini:
            init_data.append(last_fini)

        return init_data

    """
    Enforce last y value to be 0 for all timelines
    """

    def fix_last_ys(self):
        for timeline in self.timelines:
            timeline.y[-1] = 0


class BaselineProvider(UsageTimelineProvider):
    def __init__(self, job_ids: list[str], jc_df: pl.DataFrame):
        super().__init__()

        for job_id in job_ids:
            workload = jc_df.row(by_predicate=pl.col("JobId") == job_id)
            _, _, start, end, tres = workload
            start_ts = datetime.fromisoformat(start).timestamp()
            end_ts = datetime.fromisoformat(end).timestamp()
            tres_fields = dict(re.findall(r"(\w+)=(\w+)", tres))

            usage_timeline = UsageTimeline(
                [start_ts, end_ts], [int(tres_fields["cpu"]), 0]
            )
            self.timelines.append(usage_timeline)

        self.normalize_by_global_start()


class WorkflowProvider(UsageTimelineProvider):
    def __init__(self, log_files: list[str | Path], jc_df: pl.DataFrame):
        super().__init__()

        for log_file in log_files:
            steps = self.walk_workflow_log(log_file)

            times = []
            nodes = []
            for step in steps:
                name, jobid = step

                job_row = jc_df.row(by_predicate=pl.col("JobId") == jobid)
                _, _, start, end, tres = job_row
                start_ts = datetime.fromisoformat(start).timestamp()
                end_ts = datetime.fromisoformat(end).timestamp()
                tres_fields = dict(re.findall(r"(\w+)=(\w+)", tres))
                if name == "clustering": # or name == "silhouette":
                    # FIXME workflow should be rerun using the computeSlim
                    # service instead of the compute one to obtain more
                    # reliable results. Current data suggest that the job
                    # duration would not change, but it should be verified.
                    coeff = int(tres_fields["cpu"])
                elif name == "silhouette":
                    coeff = 1
                else:
                    coeff = 0

                times.extend([start_ts, end_ts])
                nodes.extend([coeff, 0])

            self.timelines.append(UsageTimeline(times, nodes))

        self.normalize_by_global_start()

    def walk_workflow_log(self, file: Path) -> list[tuple[str, str]]:
        with open(file, "r") as fp:
            data = fp.read()

        pattern = r"Scheduled job /loop/([^/]+)/\d+\.\d+ with job id (\d+)"
        res = [(i.group(1), i.group(2)) for i in re.finditer(pattern, data)]
        return res


def stacked_step_area(ax, x: list[float], ys: list[list[int]]):
    """
    Combine many step plots summing them vertically.
    x values are assumed to be increase monotonically
    """
    offsets = np.zeros(len(x))
    for idx, y in enumerate(ys):
        ax.fill_between(
            [x[0], *x, x[-1]],
            [0, *(offsets + y), 0],
            [0, *offsets, 0],
            step="post",
            alpha=1.0,
            label=f"Workload {idx+1}"
        )
        offsets += y


current_dir = Path(__file__).parent.resolve()


def main():
    mall = DMRProvider(
        [
            current_dir / "two-dmr-exec-newest/output_668.txt",
            current_dir / "two-dmr-exec-newest/output_670.txt",
        ]
    )

    jc_df = get_slurm_df(current_dir / "slurm_jobcomp.log")
    
    baseline = BaselineProvider(["976", "977"], jc_df=jc_df)

    # TODO workflow slurm logs report 3 cpus used even during silhouette
    workflow = WorkflowProvider(
        [
            current_dir.parent / "perf" / "workflow_sleep_duo1_2" / "log.txt",
            current_dir.parent / "perf" / "workflow_sleep_duo1_1" / "log.txt",
        ],
        jc_df=jc_df,
    )

    scaling = 0.85
    fig, axs = plt.subplots(3, sharex=True, figsize=(6.4 * scaling, 4.8 * scaling))

    axs[0].title.set_text("Baseline compute nodes usage")
    stacked_step_area(axs[0], *baseline.get_unified_timelines())

    axs[1].title.set_text("Workflow compute nodes usage")
    stacked_step_area(axs[1], *workflow.get_unified_timelines())

    axs[2].title.set_text("Malleability compute nodes usage")
    stacked_step_area(axs[2], *mall.get_unified_timelines())
    axs[2].legend()

    axs[0].set_yticks(np.arange(0, 3 + 1, 1))
    axs[1].set_yticks(np.arange(0, 3 + 1, 1))
    axs[2].set_yticks(np.arange(0, 3 + 1, 1))

    fig.supxlabel("Time since earliest workload start (s)")
    fig.supylabel("Number of allocated compute nodes")

    fig.tight_layout()

    fig.savefig(current_dir / "comparison_plot.pdf")


if __name__ == "__main__":
    main()
