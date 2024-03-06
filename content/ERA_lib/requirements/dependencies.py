import piplite

async def install():

    await piplite.install([
        "widgetsnbextension==4.0.10", # included (see: jupyter-lite.json)
        "nbformat==5.9.2",            # included (see: jupyter-lite.json)
        "ipywidgets==8.1.2",          # included (see: jupyter-lite.json)
        "plotly==5.19.0",             # included (see: jupyter-lite.json)
        
        "numpy==1.26.1",              # download at runtime (platform/os specific)
        "matplotlib==3.5.2",          # download at runtime (platform/os specific)
        "pandas==1.5.3",              # download at runtime (platform/os specific)
        "astropy==5.3.2"              # download at runtime (platform/os specific)
    ])
