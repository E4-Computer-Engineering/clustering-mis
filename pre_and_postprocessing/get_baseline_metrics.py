from dataclasses import dataclass
import polars as pl
from pathlib import Path
import re


@dataclass
class Execution:
    node_seconds: float
    timespan: float


def time_to_seconds(orig: str) -> float:
    minutes, seconds = orig[:-1].split("m")
    minutes = int(minutes)
    seconds = float(seconds)
    return minutes * 60 + seconds


def process_log(file: Path) -> Execution:
    with open(file, "r") as fp:
        data = fp.read()
    time_re = r"^real\t(.*)\n"
    orig_time_re = re.search(time_re, data, re.MULTILINE)
    orig_time = orig_time_re.group(1)
    time_seconds = time_to_seconds(orig_time)
    return Execution(timespan=time_seconds, node_seconds=3 * time_seconds)


def aggregate_results_single_exec(prefix: str):
    result_dir = Path(__file__).parent.parent.resolve() / "perf"
    dirnames = [f"{prefix}{i}" for i in range(1, 6)]

    executions = []
    for name in dirnames:
        file = result_dir / name / "log.txt"
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
        file1 = result_dir / first / "log.txt"
        execution1 = process_log(file=file1)
        file2 = result_dir / second / "log.txt"
        execution2 = process_log(file=file2)
        tot = Execution(
            node_seconds=execution1.node_seconds + execution2.node_seconds,
            timespan=execution1.timespan + execution2.timespan,
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
    aggregate_results_single_exec("baseline")
    aggregate_results_single_exec("baseline_sleep")

    aggregate_results_dual_exec("baseline_duo")
    aggregate_results_dual_exec("baseline_sleep_duo")


if __name__ == "__main__":
    main()
