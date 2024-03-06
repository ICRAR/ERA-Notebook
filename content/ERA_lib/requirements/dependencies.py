import piplite

async def install():

    await piplite.install([
        "widgetsnbextension==4.0.10",
        "nbformat==5.9.2",
        "ipywidgets==8.1.2",
        "numpy==1.26.4",
        "matplotlib==3.8.3",
        "pandas==2.2.1",
        "plotly==5.19.0",
        "astropy==6.0.0"
    ])
