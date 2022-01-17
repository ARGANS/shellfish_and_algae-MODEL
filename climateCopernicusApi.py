import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'grib',
        'variable': '10m_u_component_of_wind',
        'year': '2021',
        'month': '12',
        'day': '16',
        'time': '12:00',
    },
    'download.grib')