import sys
import json
from pathlib import Path
import subprocess as sp

test = dict()

test["points_filename"] = "points"
test["thread_counts"] = list(reversed([1,2,4,8,16,32,64,128,256]))
test["trials"] = 1

case1, case2, case3, case4 = dict(), dict(), dict(), dict()

case1["radius"] = 2.0
case1["split_ratio"] = 0.80
case1["min_hub_size"] = 20
case1["tag"] = "case1"

#case2["radius"] = 0.25
#case2["split_ratio"] = 0.75
#case2["min_hub_size"] = 20
#case2["tag"] = "case2"
#
#case3["radius"] = 0.25
#case3["split_ratio"] = 0.70
#case3["min_hub_size"] = 20
#case3["tag"] = "case3"
#
#case4["radius"] = 0.25
#case4["split_ratio"] = 0.65
#case4["min_hub_size"] = 20
#case4["tag"] = "case4"

test["cases"] = [case1]

with open("input.json", "w") as f:
    f.write(json.dumps(test, indent=4) + "\n")
