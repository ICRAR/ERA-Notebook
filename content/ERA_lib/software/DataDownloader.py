class DataDownloader:
    import nbformat
    import ipywidgets as widgets
    import datetime
    import asyncio
    import pyodide
    
    def __init__(self):

        self.base_url = "https://pub-c6077d5d916b42829baa6590f95bbc34.r2.dev/"

        self.DAY = None
        self.TIME = None
        self.CAL = None
        
        self.dayPicker01 = self.widgets.Dropdown(
            options=[
                ('2023/04/20 (Solar Eclipse)', "2023_04_20/"),
                ('2023/04/21', "2023_04_21/")
            ],
            value="2023_04_21/",
            description='Day',
        )

        self.timePicker01 = self.widgets.TimePicker(
            description='Time',
            value=self.datetime.time(12, 0),
            disabled=False,
            min=self.datetime.time(8, 0),
            max=self.datetime.time(17, 0)
        )

        self.calPicker01 = self.widgets.Checkbox(
            value=False,
            description='Calibrate',
            disabled=False,
            indent=True
        )

        def on_button_click(b):
            self.DAY  = self.dayPicker01.value
            self.TIME = self.timePicker01.value
            self.CAL = self.calPicker01.value
            self.determineURL(self.DAY,self.TIME,self.CAL)
            self.downloadFile(self.url,self.filename)

        self.button = self.widgets.Button(description="Download Data")
        self.button.on_click(on_button_click)

        display(self.dayPicker01)
        display(self.timePicker01)
        display(self.calPicker01)
        display(self.button)

    
    def determineURL(self,day,time,cal):

        if day == "2023_04_20/":
            timestamp = 1681920000
        elif day == "2023_04_21/":
            timestamp = 1682006400
        else:
            raise "Invalid Day"
        
        timestamp += time.hour * 3600 + time.minute * 60
        timestamp = int( ( int(timestamp/300 - 0.5) + 1 ) * 300 )
        
        if cal:
            cal_type = "fullcal"
        else:
            cal_type = "bpcal"

        self.filename = str(timestamp) + "_" + cal_type + ".uvfits"
        self.url = self.base_url + day + self.filename

    
    def downloadFile(self,url,filename):

        async def do_download():
            
            res = await self.pyodide.http.pyfetch(url)
            data = await res.bytes()
            
            with open( filename ,"wb") as f:
                f.write( data )

        loop = self.asyncio.get_event_loop()
        loop.run_until_complete( do_download() )


