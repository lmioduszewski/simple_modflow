
import random
import plotly.graph_objects as go
import numpy as np
from skimage import measure

#AVERAGE FUNCTION
def average(lst=None):
    """return the average value of a list of numbers"""
    for i in range(len(lst)):
        lst[i] = int(lst[i])
        
    return (sum(lst) / len(lst))

#DRAW LINE OF A POLYGON IN DASH FIGURE
def dpl(fig=None, x0=None, y0=None, x1=None, y1=None, color='RoyalBlue', width=3):
    """Draw a line that is part of a polygon onto a plotly plot, i.e. DrawPolyLine = dpl.
    If the line drawn closes the polygon, the polygon object is created and returned.
    Defaults to Royal Blue with a width of 3.
    The function requires x and y coordinates of both the beginning and end of the line."""
    
    polyline = fig.add_shape(type="line", x0=x0, y0=y0, x1=x1, y1=y1, 
                                line=dict(color=color,width=width))
    return polyline

#CREATE AN SVG CLOSED PATH FOR X,Y COORDINATES
def to_svg_string(coordinates):
    """Generate an SVG string from a list of (x,y) coordinate tuples

    Args:
        coordinates (list): a list of (x,y) coordinate tuples

    Returns:
        string: the returned SVG string
    """
    path_data = "".join([f" {coord[0]},{coord[1]} L" for coord in coordinates])
    path_data = "M" + path_data[:-1] + " Z"
    return path_data

#GENRATE A RANDOM RBGA COLOR
def random_rgba(a=1.0):
    """Generate a random RGBA color in the format rgba(r, g, b, a)
    with an alpha value argument. The alpha value defaults to 1.0.
    
    Args:
        a - an integer between 0 and 1 defining the transparency. 0 = transparent. 1.0 = opaque

    Returns:
        string: string of random rgba color
    """
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)
    a = a
    return f"rgba({r}, {g}, {b}, {a})"

#GENERATE A RANDOM HEX COLOR
def random_hex_color():
    """Generate a random color in hex format and returns it.

    Returns:
        str: string representing a random color in hex format.
    """
    hex_chars = "0123456789ABCDEF"
    color = "#" + "".join(random.choice(hex_chars) for _ in range(6))
    return color

#GENRATE A RANDOM RBGA COLOR
def random_rgb():
    """Generate a random RGB color in the format rgb(r, g, b)
    
    Args:
        a - an integer between 0 and 1 defining the transparency. 
        0 = transparent. 1.0 = opaque

    Returns:
        string: string of random rgb color
    """
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)
    return f"rgb({r}, {g}, {b})"

#INTERPOLTATE Z_COORDS FOR A FLAT SLOPING SURFACE
def interpolate_z_flat(x_coords, y_coords, corner_z_coords):
    """
    Interpolate the z-coordinates for all points in a grid of x,y coordinates based on the z-coordinates of the
    four corners of the rectangular grid. Returns flat lists of the x, y, and z coordinates.

    :param x_coords: a list or 1D numpy array of x-coordinates of the grid
    :param y_coords: a list or 1D numpy array of y-coordinates of the grid
    :param corner_z_coords: a list of four corner z-coordinates in the order top-left, top-right, bottom-right, bottom-left
    :return: a tuple containing the x-coordinates, y-coordinates, and z-coordinates for all points in the grid
    """
    
    top_left_z, top_right_z, bottom_right_z, bottom_left_z = corner_z_coords

    x_coords_flat = []
    y_coords_flat = []
    z_coords_flat = []
    for y in y_coords:
        for x in x_coords:
            # interpolate z-coordinate for point (x,y)
            x_frac = (x - x_coords[0]) / (x_coords[-1] - x_coords[0])
            y_frac = (y - y_coords[0]) / (y_coords[-1] - y_coords[0])
            top_z = top_left_z + (top_right_z - top_left_z) * x_frac
            bottom_z = bottom_left_z + (bottom_right_z - bottom_left_z) * x_frac
            point_z = top_z + (bottom_z - top_z) * y_frac
            
            x_coords_flat.append(x)
            y_coords_flat.append(y)
            z_coords_flat.append(point_z)

    return x_coords_flat, y_coords_flat, z_coords_flat

#INTERPOLATE Z_COORDS FOR A SURFACE AND RETURN JUST a Z-GRID
def interpolate_z_nonflat(x_coords, y_coords, corner_z_coords):
    """
    Interpolate the z-coordinates for all points in a grid of x,y coordinates based on the z-coordinates of the
    four corners of the rectangular grid.

    :param x_coords: a list of x-coordinates of the grid
    :param y_coords: a list of y-coordinates of the grid
    :param corner_z_coords: a list of four corner z-coordinates in the order top-left, top-right, bottom-right, bottom-left
    :return: a list of z-coordinates for all points in the grid
    """
    top_left_z, top_right_z, bottom_right_z, bottom_left_z = corner_z_coords

    z_coords = []
    for y in y_coords:
        row_z = []
        for x in x_coords:
            # interpolate z-coordinate for point (x,y)
            top_z = top_left_z + (top_right_z - top_left_z) * (x - x_coords[0]) / (x_coords[-1] - x_coords[0])
            bottom_z = bottom_left_z + (bottom_right_z - bottom_left_z) * (x - x_coords[0]) / (x_coords[-1] - x_coords[0])
            point_z = top_z + (bottom_z - top_z) * (y - y_coords[0]) / (y_coords[-1] - y_coords[0])
            row_z.append(point_z)
        z_coords.append(row_z)
    return z_coords

#CALCULATE Z_COORDS FOR GIVEN X,Y POINTS BASED ON GRID CORNER Z-VALUES
def calculate_z_on_surface(numrow=50, numcol=50, corner_z_coords=None, xy_points=None)->list:
    """Calculate the z-coordinates for each point in a list of x,y points,
    based on the z-coordinates of each corner of a rectangular grid. The size
    of the grid is defined by 'numrow' and 'numcol' numbers of rows and columns. 

    Args:
        numrow (int): number of rows
        numcol (int): number of columns
        corner_z_coords (list): list of z-coordinates from top left, top right, 
        bottom right, bottom left
        xy_points (list): list of lists of individual x,y points - [[x0,y0], [x1,y1], [etc]]

    Returns:
        list: list of lists of all [x,y,z] coorindates
    """
    
    top_left_z, top_right_z, bottom_right_z, bottom_left_z = corner_z_coords
    xyz_points=[]
    p1=(0, numrow-1, bottom_left_z)
    p2=(numcol-1, numrow-1, bottom_right_z)
    p3=(numcol-1, 0, top_right_z)
    p4=(0, 0, top_left_z)
    # Calculate vectors for two triangles that make up the surface
    v1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    v2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])
    v3 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])
    v4 = (p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2])
    # Calculate normal vectors for the two triangles
    n1 = (v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0])
    n2 = (v3[1]*v4[2] - v3[2]*v4[1], v3[2]*v4[0] - v3[0]*v4[2], v3[0]*v4[1] - v3[1]*v4[0])
    # Calculate the intersection point of the two planes defined by the two triangles
    d1 = -(n1[0]*p1[0] + n1[1]*p1[1] + n1[2]*p1[2])
    d2 = -(n2[0]*p3[0] + n2[1]*p3[1] + n2[2]*p3[2])
    # Loop to interate through all points in the list
    for i in range(len(xy_points)):
        x = xy_points[i][0]
        y = xy_points[i][1]
        z = (-n1[0]*x - n1[1]*y - d1)/n1[2]
        z2 = (-n2[0]*x - n2[1]*y - d2)/n2[2]
        # Average the two z-coordinates calculated from the two planes
        z_avg = (z + z2) / 2
        # Create x,y,z point and append to the xyz_points list
        point_aslist=[x, y, z_avg]
        xyz_points.append(point_aslist)
    return xyz_points

def deleteFigTraces(figure=None, 
                    attr_val_to_search=None,
                    attr_to_search=None,
                    ) -> go.Figure:
    """Delete traces based on attribute value from a plotly figure. Return the new figure.

    Args:
        figure (plotly Figure, optional): The plotly Figure with traes to delte. Defaults to None.
        attr_to_search (str, optional): Attribute to search for. Defaults to None. 
        attr_val_to_search (str, optional): Attribute value to match. Value to match can be a portion of the full attribute value. Function will attempt to match the beginning of the attribute value based on the value provided. All traces mathing that attribute value will be deleted. Defaults to None.

    Returns:
        plotly Figure: Plotly Figure with traces deleted
    """

    listfig = list(figure.data)
    if listfig==None:
        return    
    traces_to_delete = []
    attr_val_to_search = attr_val_to_search
    attr_to_search = attr_to_search
    listfig_len=len(listfig)
    if listfig_len==0:
        return
    for trace in range(listfig_len):
        if listfig[trace][attr_to_search] is not None:
            if listfig[
                trace][attr_to_search][:len(attr_val_to_search)
                                       ]==attr_val_to_search:
                traces_to_delete.append(trace)        
    traces_to_delete.sort(reverse=True)
    if traces_to_delete==None:
        return
    for trace in traces_to_delete:
        del listfig[trace]
    figure.data = tuple(listfig)
    return figure

def numTracesinFig(
    figure=None,
    attr_to_search=None,
    attr_val=None
):
    """Searches plotly figure and return the number of traces that have the given attribute with the given value.

    Args:
        figure (plotly fig, optional): plotly fig. Defaults to None.
        attr_to_search (str, optional): Attribute to search for. Defaults to None.
        attr_val (str, optional): Attribute value to search for. Defaults to None.

    Returns:
        int: number of traces with a matching attribute value
    """
    
    matching_trace_indx=[]
    numtrace=len(figure.data)
    for trace in range(numtrace):
        if figure.data[trace][attr_to_search] is not None:
            if figure.data[trace][attr_to_search][:len(attr_val)]==attr_val:
                matching_trace_indx.append(trace)
    return (len(matching_trace_indx))


def calculate_mesh_contours(vertices, faces, contour_level):
    # Create a 3D array representing the mesh
    x, y, z = np.transpose(vertices)
    volume = np.zeros(np.max(faces) + 1)
    volume[faces] = 1
    volume = np.reshape(volume, (len(x), len(y), len(z)))
    
    # Calculate the contours using the marching cubes algorithm
    contours = measure.find_contours(volume, contour_level)
    
    # Convert the 3D coordinates of the contours to 2D coordinates
    for i in range(len(contours)):
        contour = contours[i]
        contour[:, [0, 1]] = contour[:, [1, 0]]
        contour[:, 0] = x[contour[:, 0].astype(int)]
        contour[:, 1] = y[contour[:, 1].astype(int)]
        contours[i] = contour
    
    return contours

def calculate_contours(x, y, z, contour_level):
    # Create a 3D array representing the point cloud
    volume = np.zeros((len(x), len(y), len(z)))
    for i in range(len(x)):
        volume[i, :, :] = np.abs(y[:, np.newaxis] - y[i]) + np.abs(z[:, np.newaxis] - z[i])

    # Calculate the contours using the marching cubes algorithm
    verts, faces, _, _ = measure.marching_cubes(volume, contour_level=contour_level)

    # Convert the 3D coordinates of the contours to 2D coordinates
    contours = []
    for contour in measure.perimeter_edges(verts[faces]):
        c = np.zeros((len(contour), 2))
        c[:, 0] = y[contour[:, 0].astype(int)]
        c[:, 1] = z[contour[:, 1].astype(int)]
        contours.append(c)

    return contours

def flatten(l):
    return [item for sublist in l for item in sublist]