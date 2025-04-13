#!/usr/bin/env python3
from dataclasses import dataclass
import polars as pl
from pathlib import Path

# Beware:
#   This script sums up PENDING and RUNNING times for each SLURM job.
#   Other approaches from this repository consider only the RUNNING part.
#   Use get_workflow_metrics.py with a SLURM jobcomp file to get comparable results.


@dataclass
class Execution:
    node_seconds: float
    timespan: float
    start_ns: int = 0
    end_ns: int = 0


def process_log(file: Path) -> Execution:
    df = pl.read_csv(file)
    df_filtered = (
        df.filter(pl.col("name").str.starts_with("/build") == False)
        .with_columns((pl.col("end_time") - pl.col("start_time")).alias("delta"))
        .with_columns(
            pl.when(pl.col("name").str.contains("/loop/annealing"))
            .then(0)
            .when(pl.col("name").str.contains("/loop/silhouette"))
            .then(1)
            .when(pl.col("name").str.contains("/loop/clustering"))
            .then(3)
            .otherwise(0)
            .alias("coefficient")
        )
        .with_columns(
            (pl.col("delta") * pl.col("coefficient")).alias("nanoseconds-node")
        )
    )

    seconds_node_sum = df_filtered["nanoseconds-node"].sum() / 1000000000
    total_timespan = (
        df_filtered["end_time"].max() - df_filtered["start_time"].min()
    ) / 1000000000

    return Execution(
        timespan=total_timespan,
        node_seconds=seconds_node_sum,
        start_ns=df_filtered["start_time"].min(),
        end_ns=df_filtered["end_time"].max(),
    )


def aggregate_results_single_exec(prefix: str):
    result_dir = Path(__file__).parent.parent.resolve() / "perf"
    dirnames = [f"{prefix}{i}" for i in range(1, 6)]

    executions = []
    for name in dirnames:
        file = result_dir / name / "report.csv"
        execution = process_log(file=file)
        executions.append(execution)

    df = pl.DataFrame(executions)

    mean_df = df.select(
        pl.col("node_seconds").mean(),
        pl.col("node_seconds").std().alias("node_seconds_std"),
        pl.col("timespan").mean(),
        pl.col("timespan").std().alias("timespan_std"),
    )
    mean_df.write_csv(result_dir / f"{prefix}_metrics.csv")


def aggregate_results_dual_exec(prefix: str):
    result_dir = Path(__file__).parent.parent.resolve() / "perf"
    dirnames = [(f"{prefix}{i}_1", f"{prefix}{i}_2") for i in range(1, 6)]

    executions = []
    for first, second in dirnames:
        file1 = result_dir / first / "report.csv"
        execution1 = process_log(file=file1)
        file2 = result_dir / second / "report.csv"
        execution2 = process_log(file=file2)
        tot = Execution(
            node_seconds=execution1.node_seconds + execution2.node_seconds,
            timespan=(
                max(execution1.end_ns, execution2.end_ns)
                - min(execution1.start_ns, execution2.start_ns)
            )
            / 1000000000,
        )

        executions.append(tot)

    df = pl.DataFrame(executions)

    mean_df = df.select(
        pl.col("node_seconds").mean(),
        pl.col("node_seconds").std().alias("node_seconds_std"),
        pl.col("timespan").mean(),
        pl.col("timespan").std().alias("timespan_std"),
    )
    mean_df.write_csv(result_dir / f"{prefix}_metrics.csv")


def main():
    aggregate_results_single_exec("workflow")
    aggregate_results_single_exec("workflow_sleep")

    aggregate_results_dual_exec("workflow_duo")
    aggregate_results_dual_exec("workflow_sleep_duo")


if __name__ == "__main__":
    main()
