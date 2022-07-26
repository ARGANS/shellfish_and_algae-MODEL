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
    source_path:str = None
    destination_path:str = None
    year:int = None

    def __init__(self, source_path, destination_path) -> None:
        self.source_path = source_path
        self.destination_path = destination_path
        
    def parse(self, parameters_json_value:str):
        self.attrs = json.loads(parameters_json_value)

        if 'metadata' not in self.attrs:
            raise RuntimeError('Invalid input metadata schema 1')

        metadata = self.attrs['metadata']
        
        if 'year' not in metadata:
            raise RuntimeError('Invalid input metadata schema 2')

        self.year = int(self.attrs['metadata']['year'])

    @property
    def file_template(self) -> str:
        return f'{self.source_path}/_pretreated/{{param}}/{{param}}{{zone}}modelNetCDF{self.year}-01to{self.year + 1}-01.nc'
        

    @property
    def parameters(self) -> dict:
        return self.attrs['parameters']


    def isDataDownloadTaskCompleted(self) -> bool:
        try:
            _ = get_file_content(f'{self.source_path}/task.mark')
            return True
        except FileNotFoundError:
            return False

    def getMonthlySimulationsPath(self, i:int) -> str:
        return f'{self.results_dir_path}/monthly_simulations_{i:03d}.nc'

    @property
    def results_dir_path(self) -> str:
        return self.destination_path
    