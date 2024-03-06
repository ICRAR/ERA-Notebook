class FileUploader:
    import nbformat
    import ipywidgets as widgets
    
    def __init__(self):
        self.uploader = self.widgets.FileUpload(
            accept='',
            multiple=False,
            description='Select File'
        )
        display(self.uploader)
    
    def write_file(self,file_path="telescope-data.uvfits"):
        uploaded_file = self.uploader.value[0]
        with open(file_path, 'wb') as f: 
            f.write(uploaded_file.content)