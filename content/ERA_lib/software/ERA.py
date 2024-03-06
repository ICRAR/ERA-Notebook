import os
import nbformat
import ipywidgets as widgets
import numpy as np
import astropy.io.fits as fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import pandas
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from .drawWCS import drawWCSGrid

# --- Plotly Colorscale (sqrt scale) ---
viridis = px.colors.sequential.Viridis
colorscale = [
    [0.0**2, viridis[0]],
    [0.2**2, viridis[2]],
    [0.4**2, viridis[4]],
    [0.6**2, viridis[6]],
    [0.8**2, viridis[8]],
    [1.0**2, viridis[9]],
]

# -- Constants --
c = 299792458.0

class DataFile:

    def __init__(self,file_path,print_metadata=False):
        
        # -- Load FITS File --
        hdu=fits.open(file_path, memmap=False)
        
        # -- FITS Header Convenience Functions --
        NAXIS  = hdu[0].header["NAXIS"]
        PCOUNT = hdu[0].header["PCOUNT"]
        CTYPE  = [""]*NAXIS
        for i in range(len(CTYPE)):
            try:
                CTYPE[i] = hdu[0].header["CTYPE"+str(i+1)]
            except:
                pass
            
        def CRVAL(name):
            for i in range(len(CTYPE)):
                if name == CTYPE[i]:
                    return hdu[0].header["CRVAL"+str(i+1)]
        
        def CDELT(name):
            for i in range(len(CTYPE)):
                if name == CTYPE[i]:
                    return hdu[0].header["CDELT"+str(i+1)]
        
        def AxisSize(name):
            for i in range(len(CTYPE)):
                if name == CTYPE[i]:
                    return hdu[0].header["NAXIS"+str(i+1)]
        
        def AxisIndex(name):
            for i in range(len(CTYPE)):
                if name == CTYPE[i]:
                    return int( hdu[0].header["NAXIS"] - (i+1) )

        # -- Telescope Info --
        self.Nant = int( hdu[1].header["NAXIS2"] )                       # Number of antenna's
        self.Nbls = int( self.Nant * (self.Nant+1) / 2 )                 # Number of baselines
        
        # -- Observation Info & Data --
        self.Npol  = AxisSize("STOKES")                                  # Number of polarisations
        self.Nfreq = AxisSize("FREQ")                                    # Number of frequency bins
        self.freq  = np.arange(self.Nfreq)*CDELT("FREQ") + CRVAL("FREQ") # Frequency Axis
        self.freq_center = self.freq[int(self.Nfreq/2)]                  # Center freq to use after averaging
        self.Ntimes = len(np.unique(hdu[0].data[:]['DATE']))             # Number of unique timesteps in the data
        self.Nblts  = len(hdu[0].data)                                   # Total number of baselines for all timesteps
        if self.Nblts / self.Nbls != self.Ntimes:
            print("Beware: Not all baselines exist for every timestep")
        
        self.OBSRA   = CRVAL("RA")                                       # Phase Center (Right Ascension)
        self.OBSDEC  = CRVAL("DEC")                                      # Phase Center (Declination)
        
        self.UVWs = np.zeros([self.Nblts,3])
        self.UVWs[:,0] = hdu[0].data[:].par('UU') * self.freq_center     # Mult by freq to get UVW's in (Number of wavelengths)
        self.UVWs[:,1] = hdu[0].data[:].par('VV') * self.freq_center     # Mult by freq to get UVW's in (Number of wavelengths)
        self.UVdist = np.sqrt(np.sum( self.UVWs[:,0:2]**2 , axis=1))     # Compute baseline lengths (Used in figure & to compute image resolution)

        # -- Extract the visibilities as complex numbers --
        self.vis = np.zeros(self.Nblts,dtype=np.complex64)
        for i in range(self.Nblts):
            real = np.average( hdu[0].data[i][PCOUNT].take(indices=0, axis=NAXIS-1-1).take(indices=0, axis=NAXIS-1-2) ) # Take Real Part & first index of polarisation, average over all other axes
            imag = np.average( hdu[0].data[i][PCOUNT].take(indices=1, axis=NAXIS-1-1).take(indices=0, axis=NAXIS-1-2) ) # Take Imag Part & first index of polarisation, average over all other axes
            self.vis[i] = real + 1j*imag
        
        # -- Remove Autos --
        cross = np.argwhere( self.UVdist > 0 )[:,0]
        self.UVdist = self.UVdist[cross]
        self.UVWs = self.UVWs[cross,:]
        self.vis = self.vis[cross]
        self.Nblts = len(cross)
        self.Nbls = int(len(cross)/self.Ntimes)
        #Nbls = Nbls - Nant
        #Nblts = Nbls*Ntimes
        
        # -- Generate Antenna Names --
        self.ant_name = hdu[1].data.field(0)    # Antenna Names
        self.ants_in_baseline = []                   # Antenna's that make up each baseline/visibility (One could derive this from BASELINES, but the convention differs between software)
        k = 0
        for p in range(self.Nant):
            for q in range(p+1,self.Nant,1):
                self.ants_in_baseline.append( [self.ant_name[p],self.ant_name[q]] )
                k += 1
        
        # Compute image resolution & npix
        self.uvmax = np.max(self.UVdist) *(c/self.freq_center) # Max UV dist in metres
        self.lam = c / self.freq_center
        self.resolution = 1.2 * self.lam / self.uvmax
        pixres = self.resolution / 3
        self.npix = int(2/pixres)
        
        # Define image coordinate system
        self.l_axis = np.linspace(-1,1,num=self.npix)
        self.m_axis = np.linspace(-1,1,num=self.npix)
        self.l = np.tile( self.l_axis , (len(self.l_axis),1) )
        self.m = np.tile( self.m_axis , (len(self.m_axis),1) ).T
        r2 = self.l**2 + self.m**2; r2[ r2 >= 1 ] = 1
        self.n = np.sqrt( 1 - r2 )

        # Define RA,DEC coordinate system
        def gen_radec(npix,obsra,obsdec):
            # Convert to radians
            obsra  = obsra*np.pi/180.0
            obsdec = obsdec*np.pi/180.0
        
            # Define direction cosine coords - (l,m,n) grid
            grid = np.moveaxis( np.array([self.n,self.l,self.m]) , 0, -1)
        
            # -- Rotate for Dec --
            roty = np.array([
              [ np.cos(obsdec) ,0  ,np.sin(obsdec)],
              [0                ,1  ,0              ],
              [-np.sin(obsdec) ,0  ,np.cos(obsdec)]
            ])
            # -- Rotate for RA --
            rotz = np.array([
              [np.cos(-obsra) ,-np.sin(-obsra) ,0 ],
              [np.sin(-obsra) , np.cos(-obsra) ,0 ],
              [0            ,0             ,1 ]
            ])
        
            rot  = np.matmul( roty , rotz )
            grid = np.matmul( grid , rot )
        
            ra  = np.arctan2(grid[:,:,1],grid[:,:,0]) * 180/np.pi
            dec = np.arctan2(grid[:,:,2], np.sqrt( grid[:,:,0]**2 + grid[:,:,1]**2 ) ) * 180/np.pi
        
            ra[ r2 >= 1 ] = 0
            dec[ r2 >= 1 ] = 0
        
            return ra, dec
        
        self.ra_coords, self.dec_coords = gen_radec(self.npix,self.OBSRA,self.OBSDEC)

        if print_metadata:
            print("Num Antennas:\t\t",       self.Nant)
            print("Num Baselines:\t\t",      self.Nbls)
            print("Num Time Steps:\t\t",     self.Ntimes)
            print("Num Frequency Bins:\t",   self.Nfreq)
            print("Center Frequency:\t",     self.freq_center / 1e6, "MHz" )
            print("Array of frequnecies:\t", self.freq / 1e6, "MHz")
            print("Phase Center (RA):\t",    "%.3f" % self.OBSRA , "degrees")
            print("Phase Center (DEC):\t",   "%.3f" % self.OBSDEC, "degrees")
    
    
    # ------ Imaging Functions ------
    def grid_nearest(self,select=None):
        if select is None:
            select = np.arange(self.Nblts,dtype=int)
        Nvis = len(select) # The number of visibilities selected to include in the data
        VISgrid = np.zeros([self.npix,self.npix],dtype=np.complex64)
        WEIGHT  = np.zeros([self.npix,self.npix],dtype=int)
        for i in range(Nvis):
            x = int(np.round( self.UVWs[select[i],0]*2 ))
            y = int(np.round( self.UVWs[select[i],1]*2 ))
            VISgrid[y,x]   += self.vis[select[i]]
            VISgrid[-y,-x] += np.conj(self.vis[select[i]])
            WEIGHT[y,x]    += 1
            WEIGHT[-y,-x]  += 1
        WEIGHT[ WEIGHT == 0 ] = 1 # avoid dividing by 0
        VISgrid = VISgrid / WEIGHT  # Is this natural weighting?? IDK, maybe its uniform. Which one is best? I cant remember
        return VISgrid

    
    def gen_image(self,select=None):
        if select is None:
            np.arange(self.Nblts,dtype=int)
        VISgrid = self.grid_nearest(select=select)
        img = np.fft.fftshift(np.fft.fftshift(np.fft.fft2( VISgrid ),axes=0),axes=1) / self.npix
        return np.abs( img )

    
    # ------ Figures ------
    def fig_uv_coverage(self):
        """
        plot u,v coverage including all timesteps
        """
        fig = go.Figure(
            data = [
                go.Scattergl(
                  name = "",
                  x = self.UVWs[:,0]*(c/self.freq_center),
                  y = self.UVWs[:,1]*(c/self.freq_center),
                  text = [ self.ants_in_baseline[k][0] + " x " + self.ants_in_baseline[k][1] for k in range(self.Nbls) ]*self.Ntimes,
                  mode = "markers",
                  marker = dict(
                      showscale = False, # disable colorbar
                      size = 6
                  )
                ),
                go.Scattergl(
                  name = "conjugate",
                  x = -self.UVWs[:,0]*(c/self.freq_center),
                  y = -self.UVWs[:,1]*(c/self.freq_center),
                  text = [ self.ants_in_baseline[k][1] + " x " + self.ants_in_baseline[k][0] for k in range(self.Nbls) ]*self.Ntimes,
                  mode = "markers",
                  marker = dict(
                      showscale = False, # disable colorbar
                      size = 6
                  )
                )
            ],
            layout = {
                "template": "simple_white",
                "autosize": False,
                "width": 600,
                "height": 600,
                "xaxis": {
                    "title_text": "East baseline [m]",
                    "title_font": {"size": 20}
                  },
                "yaxis": {
                    "title_text": "North baseline [m]",
                    "title_font": {"size": 20}
                  },
                "showlegend": False
            }
        )
        return fig

    
    def fig_vis_v_baseline(self):
        fig = go.Figure(
            data = go.Scattergl(
                    name = "",
                    x = self.UVdist*(c/self.freq_center),
                    y = np.abs(self.vis),
                    text = [ self.ants_in_baseline[k][0] + " x " + self.ants_in_baseline[k][1] for k in range(self.Nbls) ]*self.Ntimes,
                    mode = "markers",
                    marker = dict(
                        showscale = False, # disable colorbar
                        size = 6
                    )
            ),
            layout = {
                "template": "simple_white",
                "autosize": False,
                "width": 800,
                "height": 600,
                "xaxis": {
                    "title_text": "Baseline length [m]",
                    "title_font": {"size": 20}
                  },
                "yaxis": {
                    "title_text": "Visibility amplitude [Jy]",
                    "title_font": {"size": 20}
                  },
                "showlegend": False
            }
        )
        return fig

    
    def fig_all_sky_image(self):
        # Generate Image
        img = self.gen_image()
        
        # Generate Figure
        fig = go.Figure(
            data = [
              go.Heatmap(
                name = "",
                z = img,
                x = self.l_axis,
                y = self.m_axis,
                hovertemplate="RA: %{text[0]:,.2f} <br />DEC: %{text[1]:,.2f} <br />Flux Density: %{text[2]:,.3f} Jy",
                text = [[  [ self.ra_coords[y,x] , self.dec_coords[y,x] , img[y,x] ]  for x in range(self.npix) ] for y in range(self.npix) ],
                colorscale = colorscale
              )
            ],
            layout = {
                "template": "simple_white",
                "autosize": False,
                "width": 700,
                "height": 700,
                "xaxis": {
                    "title_text": "",
                    "title_font": {"size": 20}
                  },
                "yaxis": {
                    "title_text": "",
                    "title_font": {"size": 20}
                  }
            }
        )
        
        # --- Manualy Draw WCS Grid ---
        drawWCSGrid(fig, self.OBSRA, self.OBSDEC)

        return fig

    
    def fig_img_v_baseline(self):

        # Generate data for each animation frame
        Nframes = 16
        data = np.zeros([self.npix,self.npix,Nframes])
        uvmax = np.ceil(np.max(self.UVdist) * (c/self.freq_center) ) / (c/self.freq_center)
        
        for n in range(Nframes):
            uv_cutoff = (n+1)/Nframes * uvmax
            select = np.argwhere( self.UVdist <= uv_cutoff )[:,0]
            data[:,:,n] = self.gen_image(select=select)
        
        # --- Define animation frames ---
        frames = []
        slider_steps = []
        for n in range(Nframes):
            uv_cutoff = (n+1)/Nframes * uvmax *(c/self.freq_center) # (metres)
            
            # --- define data for this frame ---
            frames.append(
                go.Heatmap(
                    name = "",
                    z = data[:,:,n],
                    x = self.l_axis,
                    y = self.m_axis,
                    colorscale = colorscale,
                    visible = False #All frames not visible by default (We will manually enable frame 0 later)
                )
            )
            
            # --- Define plot title and slider related stuff for this frame ---
            slider_step = {
              "method": "restyle", #"update",
              "label": np.round(uv_cutoff,2),
              "args": [
                        {"visible": [False] * Nframes}
                      ]
            }
            slider_step["args"][0]["visible"][n] = True
            slider_steps.append(slider_step)
        
        # --- Create Figure ---
        fig = go.Figure(
            data=frames,
            layout = {
                "template": "simple_white",
                "autosize": False,
                "width": 700,
                "height": 750,
                "xaxis": {
                    "title_text": "",
                    "title_font": {"size": 20}
                  },
                "yaxis": {
                    "title_text": "",
                    "title_font": {"size": 20}
                  },
                "sliders": [{
                    "active": 0,
                    "currentvalue": {"prefix": "Max Baseline Length (m): "},
                    "pad": {"t": 50},
                    "steps": slider_steps
                  }]
            }
        )
        
        # --- Manualy Draw WCS Grid ---
        drawWCSGrid(fig, self.OBSRA, self.OBSDEC)
        
        fig.data[0].visible = True
        
        return fig


    def fig_imgvis_v_baseline(self):
        # Generate data for each animation frame
        Nframes = 16
        fig = make_subplots(rows=1, cols=2)
        
        uvmax  = np.ceil(np.max(self.UVdist) * (c/self.freq_center) ) / (c/self.freq_center)
        vismax = np.max(np.abs(self.vis)) * 1.05
        vismin = np.min(np.abs(self.vis)) / 1.05
        
        steps = []
        for n in range(Nframes):
            uv_cutoff = (n+1)/Nframes * uvmax
            select = np.argwhere( self.UVdist <= uv_cutoff )[:,0]
            
            # --- Generate Data ---
            uv_data  = self.UVdist[select] * (c/self.freq_center)
            vis_data = np.abs(self.vis[select])
            img_data = self.gen_image(select=select)
            
            # --- Add UVplot Trace ---
            fig.add_trace(
              go.Scattergl(
                  name = "",
                  x = uv_data,
                  y = vis_data,
                  mode = "markers",
                  marker = dict(
                      color = "black",
                      showscale = False, # disable colorbar
                      size = 4
                  ),
                  visible = False #All frames not visible by default (We will manually enable frame 0 later)
              ),
              row=1, col=1
            )
            
            fig.add_trace(
              go.Heatmap(
                name = "",
                z = img_data,
                x = self.l_axis,
                y = self.m_axis,
                colorscale = colorscale,
                visible = False #All frames not visible by default (We will manually enable frame 0 later)
              ),
              row=1, col=2
            )
            
            # Define slide data
            step = {
              "method": 'restyle',
              "args": ['visible', ['legendonly'] * (2*Nframes) ],
              "label": np.round(uv_cutoff*(c/self.freq_center),2)
            }
            step['args'][1][2*n  ] = True
            step['args'][1][2*n+1] = True
            steps.append(step)
        
        fig.update_layout(width=1100, height=600,autosize=False)
        fig.layout["template"] = "simple_white"
        fig.layout["sliders"] = [{
            "steps": steps,
            "currentvalue": {"prefix": "Max Baseline Length (m): "}
        }]
        
        for n in range(Nframes):
            fig.layout["xaxis"+str(n*2+1)] = {
              "range": [0,uvmax],
              "domain": [0,0.5]
            }
            fig.layout["yaxis"+str(n*2+1)] = {
              "range": [vismin,vismax],
              "domain": [0,1]
            }
        
        fig.data[0].visible = True
        fig.data[1].visible = True
        
        return fig


    def fig_img_v_uv(self):
        # Generate data for each animation frame
        Nframes = 12
        fig = make_subplots(rows=1, cols=2)
        
        uvmax  = np.ceil(np.max(self.UVdist) * (c/self.freq_center) ) / (c/self.freq_center)
        
        steps = []
        for n in range(Nframes):
            uv_cutoff = (n+1)/Nframes * uvmax
            select = np.argwhere( self.UVdist <= uv_cutoff )[:,0]
            
            # --- Generate Data ---
            UVs = self.UVWs[select,0:2] * (c/self.freq_center)
            UVs = np.tile( UVs , (2,1) )
            UVs[int(UVs.shape[0]/2):,:] *= -1
            
            img_data = self.gen_image(select=select)
            
            # --- Add UVplot Trace ---
            fig.add_trace(
              go.Scattergl(
                  name = "",
                  x = UVs[:,0],
                  y = UVs[:,1],
                  mode = "markers",
                  marker = dict(
                      color = "black",
                      showscale = False, # disable colorbar
                      size = 4
                  ),
                  visible = False #All frames not visible by default (We will manually enable frame 0 later)
              ),
              row=1, col=1
            )
            
            fig.add_trace(
              go.Heatmap(
                name = "",
                z = img_data,
                x = self.l_axis,
                y = self.m_axis,
                colorscale = colorscale,
                visible = False #All frames not visible by default (We will manually enable frame 0 later)
              ),
              row=1, col=2
            )
            
            # Define slide data
            step = {
              "method": 'restyle',
              "args": ['visible', ['legendonly'] * (2*Nframes) ],
              "label": np.round(uv_cutoff*(c/self.freq_center),2)
            }
            step['args'][1][2*n  ] = True
            step['args'][1][2*n+1] = True
            steps.append(step)
        
        fig.update_layout(width=1100, height=650,autosize=False)
        fig.layout["template"] = "simple_white"
        fig.layout["sliders"] = [{
            "steps": steps,
            "currentvalue": {"prefix": "Max Baseline Length (m): "}
        }]
        
        for n in range(Nframes):
          fig.layout["xaxis"+str(n*2+1)] = {
              "range": [-uvmax,uvmax],
              "domain": [0,0.5]
          }
          fig.layout["yaxis"+str(n*2+1)] = {
              "range": [-uvmax,uvmax],
              "domain": [0,1]
          }
        
        fig.data[0].visible = True
        fig.data[1].visible = True
        
        return fig
        


