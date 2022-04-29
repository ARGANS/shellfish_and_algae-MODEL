import json
import os
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
    task_id = None
    
    def parse(self, parameters_json_value:str):
        self.attrs = json.loads(parameters_json_value)
        self.parse_metadata()  

    def parse_metadata(self):
        if 'metadata' not in self.attrs:
            raise RuntimeError('Invalid input metadata schema')

        metadata = self.attrs['metadata']
        
        if 'year' not in metadata or \
            not isinstance(metadata['year'], int):
            raise RuntimeError('Invalid input metadata schema')


        self.task_id = '-'.join([
            metadata['zone'],
            str(metadata['year']),
            str(metadata['depth_min']),
            str(metadata['depth_max']),
        ])

    @property
    def file_template(self) -> str:
        year:int = self.attrs['metadata']['year']
        # '/media/share/data/{zone}/{param}/{param}{zone}modelNetCDF2021-01to2022-01.nc',
        return f'/media/share/data/{self.task_id}/{{param}}/{{param}}{{zone}}modelNetCDF{year}-01to{year+1}-01.nc'
        

    @property
    def parameters(self) -> dict:
        return self.attrs['parameters']


    def isDataDownloadTaskCompleted(self) -> bool:
        try:
            _ = get_file_content(f'/media/share/data/{self.task_id}/task.mark')
            return True
        except FileNotFoundError:
            return False

        pass

    def getMonthlySimulationsPath(self, i:int) -> str:
        return f'{self.results_dir_path}/monthly_simulations_{i:03d}.nc'

    @property
    def results_dir_path(self) -> str:
        return f'/media/share/results/{self.task_id}'
    