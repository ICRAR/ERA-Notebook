import numpy as np

# ------ Plotly Convenience Functions ------
def latitude_line(lat,dec=0,N=100):
    lat = lat/180.0*np.pi
    dec = dec/180.0*np.pi
    t = np.linspace(-np.pi,np.pi,N)
    line = np.array([
        np.cos(t)*np.cos(np.ones(N)*lat),
        np.sin(t)*np.cos(np.ones(N)*lat),
        np.sin(np.ones(N)*lat)
    ]).T
    rot = np.array([
        [ np.cos(-dec) ,0  ,np.sin(-dec)],
        [0            ,1  ,0          ],
        [-np.sin(-dec) ,0  ,np.cos(-dec)]
    ])
    line = np.matmul( line , rot )

    # --- remove points that are not visible ---
    idx = np.argwhere( line[:,0] > 0 )[:,0]
    line = line[idx,:]

    # --- Generate SVG path ---
    if len(line) == 0:
        return ""
    path = f'M {line[0,1]}, {line[0,2]}'
    for k in range(1,len(line)):
        path += f'L{line[k,1]}, {line[k,2]}'
    return path

def longitude_line(lon,dec=0,N=100):
    lon = lon/180.0*np.pi
    dec = dec/180.0*np.pi
    t = np.linspace(-np.pi,np.pi,N)
    line = np.array([
        np.cos(t),
        np.zeros(N),
        np.sin(t)
    ]).T
    rotz = np.array([
        [np.cos(-lon) ,-np.sin(-lon) ,0 ],
        [np.sin(-lon) , np.cos(-lon) ,0 ],
        [0            ,0             ,1 ]
    ])
    line = np.matmul( line , rotz )
    roty = np.array([
        [ np.cos(-dec) ,0  ,np.sin(-dec)],
        [0             ,1  ,0           ],
        [-np.sin(-dec) ,0  ,np.cos(-dec)]
    ])
    line = np.matmul( line , roty )

    # --- remove points that are not visible ---
    idx = np.argwhere( line[:,0] > 0 )[:,0]
    line = line[idx,:]

    # --- Generate SVG path ---
    if len(line) == 0:
        return ""
    path = f'M {line[0,1]}, {line[0,2]}'
    for k in range(1,len(line)):
        path += f'L{line[k,1]}, {line[k,2]}'
    return path

# Manualy Draw WCS Grid
def drawWCSGrid(fig,OBSRA,OBSDEC):
    fig.add_shape(type="circle",
      xref="x", yref="y",
      x0=-1, y0=-1,
      x1=1, y1=1,
      opacity=1,
      fillcolor="rgba(0,0,0,0)",
      line_color="white",
      line_width=1
    )

    lat_resolution = 30
    for lat in np.arange(-90,90,lat_resolution):
        fig.add_shape(
            type="path",
            path=latitude_line(lat=lat,dec=OBSDEC),
            fillcolor="rgba(0,0,0,0)",
            line_color="white",
            line_width=1
        )

    lon_resolution = 30
    lon_offset = OBSRA % lon_resolution
    for lon in np.arange(-90,90,lon_resolution):
        fig.add_shape(
            type="path",
            path=longitude_line(lon=lon+lon_offset,dec=OBSDEC),
            fillcolor="rgba(0,0,0,0)",
            line_color="white",
            line_width=1
        )

