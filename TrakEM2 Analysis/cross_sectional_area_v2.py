##Jason Pipkin for use in the FIJI TrakEM2 environment
##2018-12-04
##Marder Lab

import math
import csv
from ij.text import TextWindow

#Works for the first TrakEM2 project open only!
display = Display.getFront()
project = Project.getProjects()[0]
display = Display.getFront()
layerset = display.getLayerSet()
calibration = layerset.getCalibration()
layers = layerset.getLayers()


def roi_from_area(arealist):
  #Return list of ROIs matched to the layer they came from for a given arealist
  roi_list = []
  for layer in arealist.getLayerRange():
    area = arealist.getAreaAt(layer)
    if area is not None and not area.isEmpty():
      roi = ShapeRoi(area)
      roi_list.append([roi,layer])

  return roi_list


#https://stackoverflow.com/questions/1984799/cross-product-of-two-vectors-in-python
def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c

#My own
def dot(coord_1, coord_2):
  return coord_1[0]*coord_2[0] + coord_1[1]*coord_2[1] + coord_1[2]*coord_2[2]



#https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
def PolygonArea(corners):
    if len(corners) < 3:
      return 0
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area



#https://stackoverflow.com/questions/451426/how-do-i-calculate-the-area-of-a-2d-polygon/717367#717367
def area(p):
    return 0.5 * abs(sum(x0*y1 - x1*y0
                         for ((x0, y0), (x1, y1)) in segments(p)))

def segments(p):
    return zip(p, p[1:] + [p[0]])

#https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
# Ray tracing, this function not actually used, but one could conceivably use it to test if node is inside cross-sectional polygon
def ray_tracing_method(x,y,poly):

    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside



##Manually change these values for the given branch
branchID = 973
treeID = 975
branch = project.findById(branchID)

#Un-comment the below line, choose a ball object, and put its id there.  This is for if you want to visualize the cross-sections in TrakEM2.
#testob = project.findById(992)  
t = project.findById(treeID)
dist_thresh = 3 #Good threshold.  Threshold for how far an xyz point can be from the node in order to be included.
point_thresh = 0.5 #Can go smaller, but this is fine.  Threshold for how far an xyz point can be from the uv plane in order to be included.



endNodes = t.getEndNodes()
root = t.getRoot() #Gets root of the tree
nd_list = []

#Get the list of nodes arranged as paths from tips to soma one after the other.
#Coords in xyz also included as well as tags.
for i, end in enumerate(endNodes):
  
  nodes = end.findPath(end,root)
  for nd in nodes:
    
    try:
      tag = nd.getTags().toString()
    except:
      tag = 'NaN'
      
    nd_list.append([i, (nd.getX() + t.getX())*calibration.pixelWidth, (nd.getY() + t.getY())*calibration.pixelHeight, nd.getLayer().getZ()*calibration.pixelWidth, nd, end.getTags().toString(), tag])

cross_sections = []
  


##Now we obtain the list of coordinates that defines the perimeter of the branch
##This might not include every pixel that is part of the perimeter in the case of straight lines
##Might need to adjust this to be contour-specific
#Get all the ROIs
rois = roi_from_area(branch)

contours = []
contour_id = 0

print len(rois)

#Convert each ROI to an array of points, parse that array and make a list of contours out of it
for item in rois:
  contour_coded_coords = []
  r = item[0].getShapeAsArray()
  
  contour_count = r.count(4)
  begins = []
  ends = []
  
  for idx, thing in enumerate(r):
    if thing == 0.0:
      begins.append(idx)
    if thing == 4.0:
      ends.append(idx)
    
  for j in xrange(contour_count):

    contour = []
    begin = begins[j]
    end = ends[j]
    
    for idx in xrange(begin,end):
      contour.append(r[idx])
          
    split_contour = [contour[i:i + 3] for i in xrange(0, len(contour), 3)]
    
    for entry in split_contour:
      entry = entry[1:]
      contour_coded_coords.append([contour_id]+entry+[item[1]])

    contour_id = contour_id + 1
  contours.append(contour_coded_coords)

#Organize those contours properly
print len(contours)

final_contours = []
all_points = []
for contour in contours:
  for point in contour:
    all_points.append(point)

print len(all_points)

for i in xrange(contour_id+1):
  contour = []
  for point in all_points:
    if point[0] == i:
      contour.append(point)
  final_contours.append(contour)



#Add a "distance" to each point in each contour (the distance from one point to its neighbor)
for contour in final_contours:
  for idx, point in enumerate(contour):
    x1 = point[1]
    y1 = point[2]

    if idx == 0:
      x2 = contour[-1][1]
      y2 = contour[-1][2]

    else:
      x2 = contour[idx-1][1]
      y2 = contour[idx-1][2]

    
    distance = math.sqrt((x1-x2)**2+(y1-y2)**2)
    
    point.append(distance)


print len(final_contours)

#Now unzip it all again and put into one big list
coords = []
for contour in final_contours:
  for point in contour:
    coords.append(point)

nd_cross_sects = []
##Now it's time to go through the node list and get ourselves some cross-sectional areas
for i, nd in enumerate(nd_list[:-1]):
  
  child_x = nd[1]
  child_y = nd[2]
  child_z = nd[3]

  parent_x = nd_list[i+1][1]
  parent_y = nd_list[i+1][2]
  parent_z = nd_list[i+1][3]

  parent = nd_list[i+1][4]

  #Get the path length from the parent node to the root
  path = parent.findPath(parent,root)
  path_dist = 0
  for k, path_node in enumerate(path[:-1]):
    x_i = (path_node.getX() + t.getX())*calibration.pixelWidth
    y_i = (path_node.getY() + t.getY()) * calibration.pixelHeight
    z_i = path_node.getLayer().getZ() * calibration.pixelWidth

    x_j = (path[k+1].getX() + t.getX())*calibration.pixelWidth
    y_j = (path[k+1].getY() + t.getY()) * calibration.pixelHeight
    z_j = path[k+1].getLayer().getZ() * calibration.pixelWidth

    path_dist += math.sqrt((x_j - x_i)**2 + (y_j - y_i)**2 + (z_j - z_i)**2)
  print path_dist
  
  #Put the display on the node currently being worked on - helpful for visualizing progress and fixing errors.
  coord = Coordinate(parent.getX() + t.getX(), parent.getY() + t.getY(), parent.getLayer(), parent)
  display.centerAt(coord)

  #Get a,b,c of the normal vector
  #Get plane equation stuff, and get normalized vectors describing the axes of the plane
  a = child_x - parent_x
  b = child_y - parent_y
  c = child_z - parent_z
  d = -1*(a*parent_x + b*parent_y + c*parent_z)
  magnitude = (a**2 + b**2 + c**2)**.5
  u_magnitude = (b**2 + a**2)**.5

  normal = [a/magnitude, b/magnitude, c/magnitude]
  u = [b/u_magnitude, -a/u_magnitude, 0]
  v = cross(normal, u)

  parent_u = dot(u,[parent_x,parent_y,parent_z])
  parent_v = dot(v,[parent_x,parent_y,parent_z])




  ##Now, site-wise, we find out how many of the perimeter coordinates are within 12 microns of the parent node
  membrane_coords = []
   
  for coord in coords:
    x = coord[1]*calibration.pixelWidth
    y = coord[2]*calibration.pixelHeight
    z = coord[3].getZ()*calibration.pixelWidth
    dist = math.sqrt((x-parent_x)**2 + (y-parent_y)**2 + (z-parent_z)**2)
  
    #Calculate the distance from the plane if the point is close enough to the node, then add to list
    #Note here we start with a big distance, this gets anything conceivably part of the branch
    if dist < 15:
      dist_from_plane = (a*x + b*y + c*z - (a*parent_x + b*parent_y + c*parent_z))/math.sqrt(a**2 + b**2 + c**2) 
      membrane_coords.append(coord + [x] + [y] + [z] + [dist_from_plane] + [dist])


  #print membrane_coords[0]
  #print len(membrane_coords)

  ##Now collect the membrane coords that are "on" the plane

  planar_coords = []
  simple_planar_coords = []
  uv_planar_coords = []
  
  #Initialize our max_dist value, we will add 2 back to it again just below.
  max_dist = dist_thresh - 2
  
  #This is why we got tags - make sure your thick branches are tagged 'thick' !
  #Idea is to use a looser distance threshold for these branches.
  if nd[6] == '[thick]':
    max_dist = 9 
  loop = 1

  while loop  == 1:
    #Each loop, we look out a little bit further if it hasn't found a workable polygon yet
    max_dist += 2
    for membrane_coord in membrane_coords:
  
      if abs(membrane_coord[8]) < point_thresh and abs(membrane_coord[9]) < max_dist:
        #Find the coordinate of the point on the plane that's closest to the pixel
        x = membrane_coord[5]
        y = membrane_coord[6]
        z = membrane_coord[7]
        n = (-a*x - b*y - c*z - d) / (a**2 + b**2 + c**2)
        #Coords in x,y,z of dataset
        p_x = x + a*n
        p_y = y + b*n
        p_z = z + c*n
        #Coords in u,v of the cross-sectional plane
        p_u = dot(u,[p_x,p_y,p_z])
        p_v = dot(v,[p_x,p_y,p_z])
    
        planar_coords.append(membrane_coord + [p_x] + [p_y] + [p_z] + [p_u] + [p_v])
        simple_planar_coords.append((p_x,p_y,p_z))
        uv_planar_coords.append([p_u,p_v])

    try:
    
      
      ##https://stackoverflow.com/questions/10846431/ordering-shuffled-points-that-can-be-joined-to-form-a-polygon-in-python
      # compute centroid
      cent=(sum([p[0] for p in uv_planar_coords])/len(uv_planar_coords),sum([p[1] for p in uv_planar_coords])/len(uv_planar_coords))
      # sort by polar angle
      # this re-orders stuff so that the points of the polygon are in order; necessary for area compute algorithms
      uv_planar_coords.sort(key=lambda p: math.atan2(p[1]-cent[1],p[0]-cent[0]))

  

      cross_area = PolygonArea(uv_planar_coords)
      print max_dist, cross_area
      diameter = 2*math.sqrt(cross_area/math.pi)
  
  
      ##Only append the cross-sectional area if the two nodes are in the same branch.  Removes Root as child -> end of new branch as parent situations.
      if nd_list[i][0] == nd_list[i+1][0]:
        cross_sections.append(nd_list[i+1] + [path_dist] + [cross_area] + [diameter])

        for coord in uv_planar_coords:
          nd_cross_sects.append([i, parent_u, parent_v, coord[0], coord[1], path_dist, cross_area, diameter])
          if i == 9:
            print coord[0],coord[1]
      print "Done with " + str(i+1) + " / " + str(len(nd_list))
      loop = 0
      

    except:
      continue
#Uncomment the below if you want to add balls
'''  for coord in planar_coords:
  
    
    x = coord[1]
    y = coord[2]
    layer = coord[3]

    print coord
    
    #testob.addBall(x,y,4.0,layer.getId())'''


# Convert the data into a string and output to text window from which a .csv can be manually saved
#This first one is for the cross-sectional area and diameter data
cross_data = []
for datapoint in cross_sections:
  cross_data.append(",".join(str(d) for d in datapoint))
    
csv = "\n".join(cross_data)

t = TextWindow("Cross Sections", csv, 400,400) #csv can be saved from this window.

#This second one is for the coordinates of each cross-section in the UV plane if you want to double-check that they're sensible shapes
cross_sects = []
for datapoint in nd_cross_sects:
  cross_sects.append(",".join(str(d) for d in datapoint))
    
csv = "\n".join(cross_sects)

t = TextWindow("Cross Sections", csv, 400,400)

'''for section in cross_sections:
  print section'''

'''for coord in planar_coords:
  
    
    x = coord[1]
    y = coord[2]
    layer = coord[3]
    
    testob.addBall(x,y,2.0,layer.getId())

  print "Done!"'''



'''planar_layers = []

for planar_coord in planar_coords:
  planar_layers.append(planar_coord[3])

unique_planar_layers = set(planar_layers)
print unique_planar_layers

closest_planar_coords = []

for i in xrange(len(unique_planar_layers)):
  closet_planar_coords.append([10,10])

for i,layer in enumerate(unique_planar_layers):
  for planar_coord in planar_coords:
    if planar_coord[3] == layer:
      if closest_planar_coords[i][0] == 10:'''

'''  ##Now we obtain the list of coordinates that defines the perimeter of the branch
  ##This might not include every pixel that is part of the perimeter in the case of straight lines
  ##Might need to adjust this to be contour-specific
  #Get all the ROIs
  rois = roi_from_area(branch)

  contours = []
  contour_id = 0

  print len(rois)

  #Convert each ROI to an array of points, parse that array and make a list of contours out of it
  for item in rois:
    contour_coded_coords = []
    r = item[0].getShapeAsArray()
  
    contour_count = r.count(4)
    begins = []
    ends = []
  
    for idx, thing in enumerate(r):
      if thing == 0.0:
        begins.append(idx)
      if thing == 4.0:
        ends.append(idx)
    
    for j in xrange(contour_count):

      contour = []
      begin = begins[j]
      end = ends[j]
    
      for idx in xrange(begin,end):
        contour.append(r[idx])
          
      split_contour = [contour[i:i + 3] for i in xrange(0, len(contour), 3)]
    
      for entry in split_contour:
        entry = entry[1:]
        contour_coded_coords.append([contour_id]+entry+[item[1]])

      contour_id = contour_id + 1
    contours.append(contour_coded_coords)

  #Organize those contours properly
  print len(contours)

  final_contours = []
  all_points = []
  for contour in contours:
    for point in contour:
      all_points.append(point)

  print len(all_points)

  for i in xrange(contour_id+1):
    contour = []
    for point in all_points:
      if point[0] == i:
        contour.append(point)
    final_contours.append(contour)



  #Add a "distance" to each point in each contour (the distance from one point to its neighbor)
  for contour in final_contours:
    for idx, point in enumerate(contour):
      x1 = point[1]
      y1 = point[2]

      if idx == 0:
        x2 = contour[-1][1]
        y2 = contour[-1][2]

      else:
        x2 = contour[idx-1][1]
        y2 = contour[idx-1][2]

    
      distance = math.sqrt((x1-x2)**2+(y1-y2)**2)
    
      point.append(distance)


  print len(final_contours)

  #Now unzip it all again and put into one big list
  coords = []
  for contour in final_contours:
    for point in contour:
      coords.append(point)'''

'''##STOLE THE FUCK OUT OF THE CODE BELOW 
##from https://stackoverflow.com/questions/12642256/python-find-area-of-polygon-from-xyz-coordinates

#determinant of matrix a
def det(a):
    return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[0][0]*a[1][2]*a[2][1]



#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = det([[1,a[1],a[2]],
             [1,b[1],b[2]],
             [1,c[1],c[2]]])
    y = det([[a[0],1,a[2]],
             [b[0],1,b[2]],
             [c[0],1,c[2]]])
    z = det([[a[0],a[1],1],
             [b[0],b[1],1],
             [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#dot product of vectors a and b
def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

#cross product of vectors a and b
def cross(a, b):
    x = a[1] * b[2] - a[2] * b[1]
    y = a[2] * b[0] - a[0] * b[2]
    z = a[0] * b[1] - a[1] * b[0]
    return (x, y, z)

#area of polygon poly
def area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0

    total = [0, 0, 0]
    for i in range(len(poly)):
        vi1 = poly[i]
        if i is len(poly)-1:
            vi2 = poly[0]
        else:
            vi2 = poly[i+1]
        prod = cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)'''

'''#Second, select the end-node of a branch, or wherever you want to start going backwards
lastVisited = t.getLastVisited() #Gets the actively selected node


#Get the parent of the tree, assuming you're not already at the end.  ***ADD THING TO SKIP REST OF CODE IF YOU ARE AT ROOT LATER***
if lastVisited is not root:
  parent = lastVisited.getParent()


coord = Coordinate(parent.getX() + t.getX(), parent.getY() + t.getY(), parent.getLayer(), parent)
display.centerAt(coord)

#Get child coordinates in real units
child_x = (lastVisited.getX() + t.getX()) * calibration.pixelWidth
child_y = (lastVisited.getY() + t.getY()) * calibration.pixelHeight
child_z = lastVisited.getLayer().getZ() * calibration.pixelWidth

#Get parent coordinates in real units
parent_x = (parent.getX() + t.getX()) * calibration.pixelWidth
parent_y = (parent.getY() + t.getY()) * calibration.pixelHeight
parent_z = parent.getLayer().getZ() * calibration.pixelWidth

#print child_x, child_y, child_z
#print parent_x, parent_y, parent_z'''