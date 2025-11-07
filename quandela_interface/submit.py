import os
from typing import Any
from dotenv import load_dotenv
import requests
import sys
import time


# TODO can we use floats?
def submit(token: str, QUBOmatrix: list[list[float]]) -> dict[str, str]:
    url = "https://api.cloud.quandela.com/qt/cvarvqe"

    payload = {
        "max_iterations": 10,
        "platform_name": "sim:ascella",
        "qubo_matrix": QUBOmatrix,
    }
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "Authorization": f"Bearer {token}",
    }

    # response = requests.post(url, json=payload, headers=headers).json()
    # Placeholder job_id until everything is ready
    return {"job_id": "af1cc335-a4e4-44f0-9180-e27f4daabdd8"}


def retrieve(token: str, job_id: str) -> dict[str, Any]:
    url = f"https://api.cloud.quandela.com/qt/cvarvqe/{job_id}/results"

    headers = {"Accept": "application/json", "Authorization": f"Bearer {token}"}

    return requests.get(url, headers=headers).json()


def check_status(token: str, job_id: str) -> dict[str, Any]:
    url = f"https://api.cloud.quandela.com/qt/cvarvqe/{job_id}/status"

    headers = {"Accept": "application/json", "Authorization": f"Bearer {token}"}

    return requests.get(url, headers=headers).json()


def read_qubo(file: str) -> list[list[float]]:
    res = []

    with open(file, "r") as f:
        for line in f:
            row = [float(x) for x in line.split()]
            res.append(row)

    assert(len(res) == len(res[0]))

    return res


def write_results(file: str, results_dict: dict[str, Any]):
    bitstring = results_dict["bitstring"]
    loss = results_dict["loss"]

    with open(file, "w") as f:
        f.write(",".join(bitstring))
        f.write(f"\n{loss}\n")


def wait_for_completion(token: str, job_id: str, retry_interval: int):
    timestamp_start = time.time()
    while True:
        status_dict = check_status(token, job_id)

        if "error" in status_dict:
            print(status_dict["error"])
            sys.exit(1)
        else:
            status = status_dict["status"]

        if status in ("canceled", "error"):
            print(f"Job {status}. Exiting.")
            sys.exit(1)
        elif status == "completed":
            break
        elif status == "running":
            continue

        if time.time() - timestamp_start > 1800:
            print("Timeout waiting for job to complete. Exiting.")
            sys.exit(1)

        time.sleep(retry_interval)


def main():
    load_dotenv()
    QUANDELA_TOKEN = os.getenv("QUANDELA_TOKEN")
    RETRY_INTERVAL = int(os.getenv("RETRY_INTERVAL", "10"))

    QUBO = read_qubo(sys.argv[1])

    submit_dict = submit(QUANDELA_TOKEN, QUBO)

    if "job_id" in submit_dict:
        job_id = submit_dict["job_id"]
    else:
        msg = submit_dict.get("error", submit_dict.get("message", "No message"))
        print(msg)
        sys.exit(1)

    wait_for_completion(QUANDELA_TOKEN, job_id, RETRY_INTERVAL)

    retrieve_dict = retrieve(QUANDELA_TOKEN, job_id)

    if "error" in retrieve_dict:
        print(retrieve_dict["error"])
        sys.exit(1)
    else:
        write_results(sys.argv[2], retrieve_dict["results"])


if __name__ == "__main__":
    main()
