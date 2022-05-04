import json

def import_json(path: str) -> dict:
    with open(path, 'r') as f:
        data = json.loads(f.read())
    return data 


    