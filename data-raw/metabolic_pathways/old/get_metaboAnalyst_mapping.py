import requests
import pandas as pd
import os
import json


script_dir = os.path.dirname("/mnt/c/code/playdata/data-raw/metabolic_pathways/")
os.chdir(script_dir)

file_path = "chebi_compounds_20240801_0501.tsv"

data_chebi = pd.read_csv(file_path, sep="\t", encoding="latin1")

# Extract the ID column
ids = data_chebi["ID"].tolist()


def chunk_list(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i : i + chunk_size]


output_folder = "./mapping"
os.makedirs(output_folder, exist_ok=True)

# URL and headers for the POST request
url = "https://rest.xialab.ca/api/mapcompounds"
headers = {
    "Content-Type": "application/json",
}

for i, chunk in enumerate(chunk_list(ids, 5000)):

    chunk_str = [str(item) for item in chunk]

    query_list = ";".join(chunk_str)
    payload = json.dumps({"queryList": query_list, "inputType": "chebi"})

    response = requests.post(url, data=payload, headers=headers)

    output_file = os.path.join(output_folder, f"mapped_chunk_{i + 1}.json")
    with open(output_file, "w") as f:
        f.write(response.text)
