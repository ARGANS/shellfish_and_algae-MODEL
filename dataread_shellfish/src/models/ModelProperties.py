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
        print('INPUT JSON')
        print(self.attrs)

        if 'metadata' not in self.attrs:
            raise RuntimeError('Invalid input dataset schema 1')

        datasets = self.attrs['dataset_parameters']
        
        if 'year' not in datasets:
            raise RuntimeError('Invalid input dataset schema 2')

        self.year = int(datasets['year'])

    @property
    def file_template(self) -> str:
        return f'{self.source_path}/{{Parameter}}/{{Parameter}}{{Place}}modelNetCDF{self.year}-01to{self.year + 1}-01.nc'
        

    @property
    def parameters(self) -> dict:
        return self.attrs['parameters']


    def isDataDownloadTaskCompleted(self) -> bool:
        try:
            _ = get_file_content(f'{self.source_path}/end.mark')
            return True
        except FileNotFoundError:
            return False

    def getMonthlySimulationsPath(self, i:int) -> str:
        return f'{self.destination_path}/monthly_simulations_{i:03d}.nc'
    