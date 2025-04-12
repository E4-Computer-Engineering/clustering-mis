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


def aggregate_results(prefix: str):
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


def main():
    aggregate_results("baseline")
    aggregate_results("baseline_sleep")


if __name__ == "__main__":
    main()
