from dataclasses import dataclass
import polars as pl
from pathlib import Path


@dataclass
class Execution:
    node_seconds: float
    timespan: float


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

    return Execution(timespan=total_timespan, node_seconds=seconds_node_sum)


def aggregate_results(prefix: str):
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

def main():
    aggregate_results("workflow")
    aggregate_results("workflow_sleep")


if __name__ == "__main__":
    main()
