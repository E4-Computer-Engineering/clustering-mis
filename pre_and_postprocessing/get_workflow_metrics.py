from datetime import datetime
import re
from dataclasses import dataclass
import polars as pl
from pathlib import Path


@dataclass
class Execution:
    node_seconds: float
    timespan: float
    start_s: float = 0.0
    end_s: float = 0.0


def get_slurm_df(file: Path) -> pl.DataFrame:
    rows = []

    with open(file, "r") as f:
        for line in f:
            fields = dict(re.findall(r"(\S+?)=(\S+)", line))
            rows.append(fields)

    df = pl.DataFrame(rows)
    fine_df = df.select(
        pl.col("JobId"), pl.col("JobState"), pl.col("StartTime"), pl.col("EndTime")
    ).filter(pl.col("JobState") == "COMPLETED")

    return fine_df


def walk_log(file: Path) -> list[tuple[str, str]]:
    with open(file, "r") as fp:
        data = fp.read()

    pattern = r"Scheduled job /loop/([^/]+)/\d+\.\d+ with job id (\d+)"
    res = [(i.group(1), i.group(2)) for i in re.finditer(pattern, data)]
    return res


def process_log_with_jobcomp(file: Path, jc_df: pl.DataFrame) -> Execution:
    steps = walk_log(file)

    starts = []
    ends = []
    load = 0.0

    for step in steps:
        name, jobid = step
        if name == "clustering":
            coeff = 3
        elif name == "silhouette":
            coeff = 1
        elif name == "annealing":
            coeff = 0
        else:
            print("Unexpected input!")

        job_row = jc_df.row(by_predicate=pl.col("JobId") == jobid)
        _, _, start, end = job_row
        start_ts = datetime.fromisoformat(start).timestamp()
        end_ts = datetime.fromisoformat(end).timestamp()

        load += coeff * (end_ts - start_ts)
        starts.append(start_ts)
        ends.append(end_ts)

    global_start = min(starts)
    global_end = max(ends)

    return Execution(
        node_seconds=load,
        timespan=(global_end - global_start),
        start_s=global_start,
        end_s=global_end,
    )


def aggregate_results_single_exec(prefix: str, jc_df: pl.DataFrame):
    result_dir = Path(__file__).parent.parent.resolve() / "perf"
    dirnames = [f"{prefix}{i}" for i in range(1, 6)]

    executions = []
    for name in dirnames:
        file = result_dir / name / "log.txt"
        execution = process_log_with_jobcomp(file=file, jc_df=jc_df)
        executions.append(execution)

    df = pl.DataFrame(executions)

    mean_df = df.select(
        pl.col("node_seconds").mean(),
        pl.col("node_seconds").std().alias("node_seconds_std"),
        pl.col("timespan").mean(),
        pl.col("timespan").std().alias("timespan_std"),
    )
    mean_df.write_csv(result_dir / f"{prefix}_metrics.csv")


def aggregate_results_dual_exec(prefix: str, jc_df: pl.DataFrame):
    result_dir = Path(__file__).parent.parent.resolve() / "perf"
    dirnames = [(f"{prefix}{i}_1", f"{prefix}{i}_2") for i in range(1, 6)]

    executions = []
    for first, second in dirnames:
        file1 = result_dir / first / "log.txt"
        execution1 = process_log_with_jobcomp(file=file1, jc_df=jc_df)
        file2 = result_dir / second / "log.txt"
        execution2 = process_log_with_jobcomp(file=file2, jc_df=jc_df)
        tot = Execution(
            node_seconds=execution1.node_seconds + execution2.node_seconds,
            timespan=(
                max(execution1.end_s, execution2.end_s)
                - min(execution1.start_s, execution2.start_s)
            ),
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
    jc_path = Path("slurm_jobcomp.log")
    if not jc_path.exists():
        print(
            "This script needs a slurm_jobcomp.log file (look under /var/log/ on the controller node)"
        )
        exit(1)

    jc_df = get_slurm_df(file=jc_df)

    aggregate_results_single_exec("workflow", jc_df)
    aggregate_results_single_exec("workflow_sleep", jc_df)

    aggregate_results_dual_exec("workflow_duo", jc_df)
    aggregate_results_dual_exec("workflow_sleep_duo", jc_df)


if __name__ == "__main__":
    main()
