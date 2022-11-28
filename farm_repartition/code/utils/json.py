import json
from json import loads


def import_json(path: str) -> dict:
    with open(path, 'r') as f:
        data = json.loads(f.read())
    return data 

def dump_json(data:dict, path:str):
    with open(path, 'w') as f:
        json.dump(data, f)
