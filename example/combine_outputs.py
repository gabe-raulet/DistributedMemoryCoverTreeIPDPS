import json
import sys
from pathlib import Path

paths = list(Path(".").glob("out.*.json"))
json_dicts = [json.load(open(str(p), "r")) for p in paths]
json_dict = dict({"experiments" : json_dicts})
json.dump(json_dict, sys.stdout, indent=4)
for p in paths: p.unlink()
