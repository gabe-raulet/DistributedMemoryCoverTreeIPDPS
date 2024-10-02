import sys
import json
from pathlib import Path
import subprocess as sp

test = dict()

test["points_filename"] = "points"
test["thread_counts"] = [8,16,32,64,128,256]
test["trials"] = 3

case1, case2, case3 = dict(), dict(), dict()

case1["radius"] = 0.1
case1["split_ratio"] = 0.5
case1["min_hub_size"] = 10
case1["tag"] = "case1"

test["cases"] = [case1]

with open("input.json", "w") as f:
    f.write(json.dumps(test, indent=4) + "\n")
