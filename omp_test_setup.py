import sys
import json
from pathlib import Path
import subprocess as sp

test = dict()

test["points_filename"] = "points"
test["thread_counts"] = [1,2,3,4,5,6,7,8,9,10,11,12]
test["trials"] = 3

case1, case2, case3 = dict(), dict(), dict()

case1["radius"] = 0.5
case1["split_ratio"] = 0.5
case1["min_hub_size"] = 10
case1["tag"] = "case1"

case2["radius"] = 0.5
case2["split_ratio"] = 0.7
case2["min_hub_size"] = 10
case2["tag"] = "case2"

case3["radius"] = 0.5
case3["split_ratio"] = 0.9
case3["min_hub_size"] = 10
case3["tag"] = "case3"

test["cases"] = [case1, case2, case3]

with open("input.json", "w") as f:
    f.write(json.dumps(test, indent=4) + "\n")
