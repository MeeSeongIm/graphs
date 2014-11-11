# download python imaging library (PIL) for Python 3*: http://www.pythonware.com/products/pil/
# to be used with percolation Monte Carlo py files. 

import numpy, os, random
import Image, ImageDraw
import tempfile

# some graphics for networks

# -----------------------------------------------------------------------

# image display

def Display(image_file='tmpf.jpg'):
    os_name = os.name
    if os_name == 'nt': # Windows
        try:
            os.system('mspaint %s &' % image_file)
        except:
            raise OSError("Cannot display %s with Windows mspaint" % image_file)
    else:
        os_uname = os.uname()
        if os_uname[0] == 'Darwin':  # Mac OS X, assume no X server running
            try:
                os.system('open /Applications/Preview.app %s &' % image_file)
            except:
                raise OSError("Cannot display %s with Preview application" % image_file)

        elif os_name == 'posix':  # Linux, Unix, etc.
            try:
                os.system('display %s &' % image_file)
            except:
                raise OSError("Cannot display %s with ImageMagick display. ImageMagick display requires a running X server." % image_file)
        else:
            raise OSError("no known display function for OS %s" % os_name)


# -----------------------------------------------------------------------

#: circular graph layouts

def GenerateCircleGraphImage(graph, imfile, windowSize=800, 
			     dotsize=4, windowMargin = 0.02):
    """creates an image file in the specified imfile and
    returns the PIL image object, laying out the specified graph in a
    circle."""
    im = Image.new('RGB', (windowSize,windowSize))
    draw = ImageDraw.Draw(im)
    color = (255,255,255)
    center = windowSize/2.
    radius = (1.-2*windowMargin)*windowSize/2.

    all_nodes = graph.GetNodes()
    all_nodes.sort()
    L = len(all_nodes)
    # Create a dictionary that maps nodes to their positions around the circle
    nodePosition = {}
    for index, node in enumerate(all_nodes):
        theta = 2.*numpy.pi*float(index)/L
        x = radius * numpy.cos(theta) + center
        y = radius * numpy.sin(theta) + center
        nodePosition[node] = (x, y)
        draw.ellipse( ((x-dotsize/2, y-dotsize/2), (x+dotsize/2, y+dotsize/2)), fill=color )

    # Draw the lines between the nodes
    for node in all_nodes:
        neighbors = graph.GetNeighbors(node)
        for neighbor in neighbors:
	    # We want to draw bonds only once, even though two bonds connect
	    #   node and neighbor.  We can test "if neighbor > node" to
	    #   implement this.  For any pair of objects, Python will
	    #   consistently define this operation (e.g., if x1<x2 is True,
	    #   then x1>x2 is False).  For arbitrary node IDs, we may not
	    #   know what it means for ID1 to be greater than ID2, but it
	    #   does not matter for the purpose here.  For Python classes,
	    #   one define the special __cmp__ (compare) method that
	    #   indicates how to compare (>,<,==) two instances of the
	    #   class.
                if neighbor > node:
                    draw.line((nodePosition[node], nodePosition[neighbor]), fill=color)

    im.save(imfile)
    return im

def GenerateCircleGraphImageWithWeights(graph, imfile, 
					edge_weights=None, node_weights=None,
					edge_scale=None, node_scale=None,
					windowSize=800, 
					dotscale=4., linescale=2., 
					magnification=4, windowMargin = 0.02):
    """creates an image file in the specified imfile and
    returns the PIL image object, laying out the specified graph in a
    circle. Nodes are drawn with radii proportional to the specified
    node_weights divided by node_scale, and edges with thicknesses
    proportional to the specified edge_weights divided by edge_scale.
    If node_scale and edge_scale are not supplied, they are assumed to be
    equation to the number of nodes in the system.
    node_weights is a dictionary mapping node ids to weights.
    edge_weights is a dictionary mapping edge ids (tupled pairs of node ids)
    to weights."""
    # Make big image, shrink with antialiasing
    m = magnification
    bigSize = (m*windowSize, m*windowSize)
    imBig = Image.new('RGB', bigSize, (255,255,255)) # White background
    draw = ImageDraw.Draw(imBig)
    edge_color = (0,0,0)
    node_color = (255,0,0)
    center = windowSize/2.
    radius = (1.-2*windowMargin)*windowSize/2.
    # Nodes for SWG are integers 0<=node<L
    all_nodes = graph.GetNodes()
    all_nodes.sort()
    L = len(all_nodes)
    if node_scale is None:
        node_scale = float(L)
    if edge_scale is None:
        edge_scale = float(L)
    # Create a dictionary that maps nodes to their positions around the circle
    nodePosition = {}
    for index, node in enumerate(all_nodes):
        theta = 2.*numpy.pi*float(index)/L
        x = radius * numpy.cos(theta) + center
        y = radius * numpy.sin(theta) + center
        nodePosition[node] = (x, y)
    for node in all_nodes:
        x, y = nodePosition[node]
        neighbors = graph.GetNeighbors(node)
        for neighbor in neighbors:
	    # Only draw bonds once (see above)
            if neighbor > node:
                xNbr, yNbr = nodePosition[neighbor]
                # Thickness of line proportional to edge_scale
                if (edge_weights):
                    linewidth = linescale * numpy.sqrt(edge_weights[(node, neighbor)]/edge_scale)
                else:
                    linewidth = linescale
                # Draw rectangles to get thick lines
                perpLength = numpy.sqrt((y-yNbr)**2+(x-xNbr)**2)
                perpx = ((y-yNbr)/perpLength)*linewidth/2
                perpy = (-(x-xNbr)/perpLength)*linewidth/2
                polyFromLine = ((m*(x+perpx), m*(y+perpy)), (m*(x-perpx), m*(y-perpy)), (m*(xNbr-perpx), m*(yNbr-perpy)), (m*(xNbr+perpx), m*(yNbr+perpy)))
                draw.polygon(polyFromLine, fill=edge_color )
    for node in all_nodes:
        x, y = nodePosition[node]
        # Make size of dot proportional to square root of weight/L
        # (e.g., small world network betweenness ~ L)
        if (node_weights):
            dotsize = dotscale * numpy.sqrt(node_weights[node]/node_scale)
        else:
            dotsize = dotscale
        draw.ellipse( ((m*(x-dotsize/2), m*(y-dotsize/2)), (m*(x+dotsize/2), m*(y+dotsize/2))), fill=node_color )
    im = imBig.resize((windowSize,windowSize), Image.ANTIALIAS)
    im.save(imfile)
    return im


def DisplayCircleGraphSimple(graph, windowSize=800, dotsize=4, windowMargin = 0.02):
    import tempfile
    filename = tempfile.mktemp()  # make unique temporary filename in /tmp
    filename += ".png"
    GenerateCircleGraphImage(graph, filename, windowSize, dotsize,
			     windowMargin)
    Display(filename)

def DisplayCircleGraph(graph, 
		       edge_weights=None, node_weights=None,
		       edge_scale=None, node_scale=None,
		       windowSize=800, 
		       dotscale=4., linescale=2.,
		       magnification=4, windowMargin=0.02):
    import tempfile
    filename = tempfile.mktemp()  # make unique temporary filename in /tmp
    filename += ".png"
    GenerateCircleGraphImageWithWeights(graph, imfile=filename, 
    					edge_weights=edge_weights,
					node_weights=node_weights, 
					edge_scale=edge_scale, 
					node_scale=node_scale,
					windowSize=windowSize, 
					dotscale=dotscale,
					linescale=linescale,
					magnification=magnification,
					windowMargin=windowMargin)
    Display(filename)
    
# -----------------------------------------------------------------------

# 2D percolation graphics

def DrawSquareNetworkBonds(graph, nodelists=None,
			   dotsize=None, linewidth=None, 
			   imsize=800, windowMargin=0.02, imfile=None):
    """DrawSquareNetworkBonds(g) will draw an image file of the 2D 
    square--lattice bond percolation network g, with bonds and sites shown,
    and then will display the result.
    DrawSquareNetworkBonds(g,nodelists) for a percolation graph g and a list of
    node clusters nodelists = [[node,node],[node],...] will draw the
    first cluster black, and the rest each in a random color.
    By default, the image file will be stored in a uniquely named png file
    in /tmp, although the image file name can be supplied optionally with
    the imfile argument.
    A node is a tuple (i,j): DrawSquareNetworkBonds will display it with i
    labeling the horizontal axis and j the vertical, with (0,0) in the
    upper right hand corner. [This is the transpose of the matrix
    convention, so (column, row); it's flipping the vertical axis, so
    it's like (x, -y) ]"""
    # Set up cluster of all nodes in network if no clusters given
    if nodelists is None:
        nodelists = [graph.GetNodes()]
    # Set up image file
    if imfile is None:
        imfile = tempfile.mktemp()  # make unique filename in /tmp
        imfile += "_square_network_bonds.png"
    white = (255,255,255)  # background color
    im = Image.new('RGB', (imsize,imsize), color=white)
    draw = ImageDraw.Draw(im)
    # Nodes = (ix, iy) running from (0,0) to (L-1,L-1)
    # Won't always work for site percolation: 
    # Assumes entire row and column of nodes not missing
    L = max(max([node[0] for node in graph.GetNodes()]), \
	    max([node[1] for node in graph.GetNodes()])) + 1.0
    # Default dot size and line width depends on L
    if dotsize is None:
    	dotsize = max((1-2*windowMargin)*imsize/(4*L),1)
    if linewidth is None:
    	linewidth = max((1-2*windowMargin)*imsize/(10*L),1)
    # Start colors with black
    color = (0,0,0)
    # Draw clusters
    for cluster in nodelists:
        def ScreenPos(i):
            return (windowMargin + ((i+0.5)/L)*(1-2*windowMargin))*imsize
        # Find screen location (sx,sy) for node
        for node in cluster:
            ix, iy = node;	# node = (ix,iy) running from (0,0) to (L-1,L-1)
            sx = ScreenPos(ix)
            sy = ScreenPos(iy)
            draw.ellipse( ((sx-dotsize/2,sy-dotsize/2), (sx+dotsize/2,sy+dotsize/2)), fill=color )
        # Define function to draw thick line
            def DrawThickLine(sx1, sy1, sx2, sy2):
                perpLength = numpy.sqrt((sy2-sy1)**2+(sx2-sx1)**2)
                perpx = ((sy2-sy1)/perpLength)*linewidth/2
                perpy = (-(sx2-sx1)/perpLength)*linewidth/2
                polyFromLine = ((sx1+perpx, sy1+perpy),
				(sx1-perpx, sy1-perpy),
				(sx2-perpx, sy2-perpy),
				(sx2+perpx, sy2+perpy))
                draw.polygon(polyFromLine, fill=color)
            neighbors = graph.GetNeighbors(node) # Find neighbors
            for neighbor in neighbors:
                if neighbor <= node: # Draw each bond once: only if i>j
                    continue
                # Find screen location (sxNbr,syNbr)) for edge
                ixNbr, iyNbr = neighbor
                sxNbr = (windowMargin+((ixNbr+0.5)/L)*(1-2*windowMargin))*imsize
                syNbr = (windowMargin+((iyNbr+0.5)/L)*(1-2*windowMargin))*imsize
                
                if (ix == 0) & (ixNbr == L-1):      # periodic boundary conditions
                    sxMinusHalf = ScreenPos(ix-0.5) # Only nearest neighbor bonds implemented
                    DrawThickLine(sx,sy, sxMinusHalf,syNbr)
                    sxNbrPlusHalf = ScreenPos(ixNbr+0.5)
                    DrawThickLine(sxNbrPlusHalf,sy, sxNbr,syNbr)
                elif (ix==L-1) & (ixNbr == 0):
                    sxPlusHalf = ScreenPos(ix+0.5)
                    DrawThickLine(sx,sy, sxPlusHalf,syNbr)
                    sxNbrMinusHalf = ScreenPos(ixNbr-0.5)
                    DrawThickLine(sxNbrMinusHalf,sy, sxNbr,syNbr)
                elif (iy==0) & (iyNbr == L-1):
                    syMinusHalf = ScreenPos(iy-0.5)
                    DrawThickLine(sx,sy, sxNbr,syMinusHalf)
                    syNbrPlusHalf = ScreenPos(iyNbr+0.5)
                    DrawThickLine(sx,syNbrPlusHalf, sxNbr,syNbr)
                elif (iy==L-1) & (iyNbr==0):
                    syPlusHalf = ScreenPos(iy+0.5)
                    DrawThickLine(sx,sy, sxNbr, syPlusHalf)
                    syNbrMinusHalf = ScreenPos(iyNbr-0.5)
                    DrawThickLine(sx,syNbrMinusHalf, sxNbr,syNbr)
                else:
                    DrawThickLine(sx,sy, sxNbr,syNbr)
        colorRange = (0, 200) # pick random color
    color = (random.randint(*colorRange), random.randint(*colorRange), random.randint(*colorRange))
    im.save(imfile)
    Display(imfile)
    return im

def DrawSquareNetworkSites(graph, nodelists=None, scale=0, imsize=800, imfile=None):
    """DrawSquareNetworkSites(g) will draw an image file of the 2D 
    square--lattice percolation network g, with sites shown,
    and then will display the result.
    DrawSquareNetworkSites(g, nodelists) for a percolation graph g and
    a list of node clusters nodelists = [[node,node],[node],...] will draw the
    first cluster black, and the rest each in a random color.
    By default, the image file will be stored in a uniquely named png file
    in /tmp, although the image file name can be supplied optionally with
    the imfile argument.
    A node is a tuple (i,j): DrawSquareNetworkSites will display it with i
    labeling the horizontal axis and j the vertical, with (0,0) in the
    upper right hand corner. [This is the transpose of the matrix
    convention, so (column, row); it's flipping the vertical axis, so
    it's like (x, -y) ]"""
    # Set up cluster of all nodes in network if no clusters given
    if nodelists is None:
        nodelists = [graph.GetNodes()]
    # Set up image file
    if imfile is None:
        imfile = tempfile.mktemp()  # make unique filename in /tmp
        imfile += "_square_network_sites.png"
    L = max(max([node[0] for node in graph.GetNodes()]), max([node[1] for node in graph.GetNodes()])) + 1
    if (scale==0):
        scale = max(1,int(imsize/L)) # Size of squares for each node
    # make the background white 
    white = (255,255,255)
    im = Image.new('RGB', (scale*L, scale*L), white)
    if (scale>1):
        draw = ImageDraw.Draw(im)
    # Nodes = (ix, iy) running from (0,0) to (L-1,L-1)
    # Won't always work for site percolation: 
    # Assumes entire row and column of nodes not missing
    color = (0,0,0) # starting color
    # Draw clusters
    for cluster in nodelists:
        if (scale==1):
            for node in cluster:
                im.putpixel(node, color)
        else:
            for node in cluster:
                x = node[0]*scale
                y = node[1]*scale
                draw.rectangle(((x,y),(x+scale,y+scale)), fill=color)
        colorRange = (0, 200)  # pick random color for next cluster 
        color = (random.randint(*colorRange), random.randint(*colorRange), random.randint(*colorRange))
    im.save(imfile)
    Display(imfile)
    return im

def DrawTriangularNetworkSites(graph, nodelists=None, L=0, 
			     scale=0, magnification=4, imsize=800,
			     imfile=None):
    """DrawTriangularNetworkSites(g) will draw an image file of the 2D 
    triangle--lattice percolation network g, with sites shown,
    and then will display the result.
    DrawTriangularNetworkSites(g, nodelists) for a percolation graph g and
    a list of node clusters nodelists = [[node,node],[node],...] will draw the
    first cluster black, and the rest each in a random color.
    By default, the image file will be stored in a uniquely named png file
    in /tmp, although the image file name can be supplied optionally with
    the imfile argument."""
    if nodelists is None:           # set-upcluster of all nodes in network if no clusters are given 
        nodelists = [graph.GetNodes()]
    if imfile is None:              # set-up image file
        imfile = tempfile.mktemp()  # make unique filename in /tmp
        imfile += "_square_network_sites.png"
    if L==0:                        # if L =0, guess L from size of array
        L = max(max([node[0] for node in graph.GetNodes()]), max([node[1] for node in graph.GetNodes()])) + 1
    # if scale is zero, set it to a size so the system is roughly imsize
    if (scale==0):
        scale = max(1,int(imsize/L)) # Size of triangles for each node
    # Background white (missing nodes)
    white = (255,255,255)
    # Basic idea: draw system magnification*big, and shrink at end
    bigHeightx = int(round(magnification*scale*L))
    bigHeighty = int(round(magnification*scale*L*numpy.sqrt(3.)/2.))
    imBig = Image.new('RGB', (bigHeightx, bigHeighty), white)
    draw = ImageDraw.Draw(imBig)
    # Nodes = (ix, iy) running from (0,0) to (L-1,L-1)
    # Won't always work for site percolation: 
    # Assumes entire row and column of nodes not missing
    color = (0,0,0) # starting color
    # Basis for center of hexagon
    xhat = magnification*scale*numpy.array([1.,0.])
    yhat = magnification*scale*numpy.array([1./2.,numpy.sqrt(3.)/2.])
    # Shape of polygon
    up = magnification*scale*\
    		numpy.array([0.,1./numpy.sqrt(3.)])
    upright = magnification*scale*\
    		numpy.array([1./2.,1./(2.*numpy.sqrt(3.))])
    upleft = magnification*scale*\
    		numpy.array([-1./2.,1./(2.*numpy.sqrt(3.))])
    hexagon = numpy.array([upright,up,upleft,-upright,-up,-upleft])
    # Periodic boundary condition
    pbc = magnification*scale*L*numpy.array([1., numpy.sqrt(3.)/2.])
    # Rounding error check
    eps = magnification*scale*1.e-6
    # Draw clusters
    for cluster in nodelists:
        for node in cluster:
            xc = (node[0]*xhat+node[1]*yhat) % pbc
            vertices = hexagon + xc
            draw.polygon(tuple(map(tuple, vertices)), fill=color)
            # Check if polygon crosses boundary
            # Modulo doesn't work for negatives
            if (not numpy.allclose( vertices%pbc, vertices )) or (vertices.min() < 0.):
                xcPBC = xc + numpy.array([pbc[0],0.])           # plot polygon plus pbc vectors
                vertices = hexagon + xcPBC
                draw.polygon(tuple(map(tuple, vertices)), fill=color)
                xcPBC = xc + numpy.array([pbc[0]/2.,pbc[1]])
                vertices = hexagon + xcPBC
                draw.polygon(tuple(map(tuple, vertices)), fill=color)
                xcPBC = xc + numpy.array([-pbc[0]/2.,pbc[1]])
                vertices = hexagon + xcPBC
                draw.polygon(tuple(map(tuple, vertices)), fill=color)
        colorRange = (50, 225)                                  # pick random color for next cluster
        color = (random.randint(*colorRange), random.randint(*colorRange), random.randint(*colorRange))
    heightx = int(round(scale*L))
    heighty = int(round(scale*L*numpy.sqrt(3.)/2.))
    im = imBig.resize((heightx,heighty), Image.ANTIALIAS)
    im.save(imfile)
    Display(imfile)
    return im
 
