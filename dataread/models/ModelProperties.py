import json
import os
from pprint import pprint
__all__ = ['ModelProperties']

def get_file_content(path:str) -> str:
    """
    raises FileNotFoundError exception
    """

    with open(path, 'r') as file:
        data = file.read()
        return data

class ModelProperties():
    attrs = {}
    task_id:str = None
    dataset_id:str = None
    year:int = None

    def __init__(self, dataset_id, task_id) -> None:
        self.dataset_id = dataset_id
        self.task_id = task_id
    
    def parse(self, parameters_json_value:str):
        self.attrs = json.loads(parameters_json_value)
        self.parse_metadata()  

    def parse_metadata(self):
        print('[CALL parse_metadata]')
        pprint(self.attrs)
        if 'metadata' not in self.attrs:
            raise RuntimeError('Invalid input metadata schema 1')

        metadata = self.attrs['metadata']
        
        if 'year' not in metadata:
            raise RuntimeError('Invalid input metadata schema 2')

        self.year = int(self.attrs['metadata']['year'])

        self.dataset_id = '-'.join([
            metadata['zone'],
            str(metadata['year']),
            str(metadata['depth_min']),
            str(metadata['depth_max']),
        ])

    @property
    def file_template(self) -> str:
        return f'/media/share/data/{self.dataset_id}/{{param}}/{{param}}{{zone}}modelNetCDF{self.year}-01to{self.year + 1}-01.nc'
        

    @property
    def parameters(self) -> dict:
        return self.attrs['parameters']


    def isDataDownloadTaskCompleted(self) -> bool:
        try:
            _ = get_file_content(f'/media/share/data/{self.dataset_id}/task.mark')
            return True
        except FileNotFoundError:
            return False

        pass

    def getMonthlySimulationsPath(self, i:int) -> str:
        return f'{self.results_dir_path}/monthly_simulations_{i:03d}.nc'

    @property
    def results_dir_path(self) -> str:
        return f'/media/share/results/{self.dataset_id}/{self.task_id}'
    