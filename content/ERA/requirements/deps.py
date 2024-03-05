import piplite

async def install():
    await piplite.install("https://files.pythonhosted.org/packages/99/bc/82a8c3985209ca7c0a61b383c80e015fd92e74f8ba0ec1af98f9d6ca8dce/widgetsnbextension-4.0.10-py3-none-any.whl")
    await piplite.install(["nbformat","ipywidgets","numpy","matplotlib","pandas","plotly","astropy"])
