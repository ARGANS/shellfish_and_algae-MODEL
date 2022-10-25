from datetime import datetime

def ncdf_date():
    # GDL> ncdf_date()
    # 2022-09-15T19:37:08Z  
    return  datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")

def hours_since_1900(year):
    """
    hours since 1900-1-1 0:0:
    """
    return int((datetime(year,1,1).timestamp() - datetime(1900,1,1).timestamp()) / 3600)
