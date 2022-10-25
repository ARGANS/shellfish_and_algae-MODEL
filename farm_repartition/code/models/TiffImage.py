from dataclasses import dataclass

@dataclass
class TiffImage:
    lonmin:float
    lonmax:float
    latmin:float
    latmax:float
    pasx:float
    pasy:float