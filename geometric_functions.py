import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.interpolate import interp2d
from stl import mesh
import stl
import time
import random

def naca_4digit(number,res):


    # Get Airfoil Parameters

    m = float(number[0])/100.0
    p = float(number[1])/10.0
    t = float(number[2:])/100.0

    # Cosine Point Spacing

    b_front = np.linspace(0,np.arccos(1-2*p),res)
    b_back = np.linspace(np.arccos(1-2*p),np.pi,res)[1:]

    front = (1-np.cos(b_front))/2.0
    back = (1-np.cos(b_back))/2.0

    x = np.append(front,back)

    # Calculate Thickness Profile

    a0 = 0.2969
    a1 = 0.126
    a2 = 0.3516
    a3 = 0.2843
    a4 = 0.1036

    yt = (t/0.2)*(a0*np.sqrt(x) - a1*x - a2*x**2 + a3*x**3 - a4*x**4)

    # Generate Mean Camber Line

    yc_front = (m/p**2)*(2*p*front - front**2)
    yc_back = (m/(1-p)**2)*((1-2*p) + 2*p*back - back**2)

    yc = np.append(yc_front,yc_back)

    # Generate Slope Profile

    dy_front = ((2*m)/p**2)*(p-front)
    dy_back = ((2*m)/(1-p)**2)*(p-back)

    dy_dx = np.append(dy_front,dy_back)
    theta = np.arctan(dy_dx)

    # Calculate Airfoil Coordinates

    x_u = x - yt*np.sin(theta)
    y_u = yc + yt*np.cos(theta)

    x_l = x + yt*np.sin(theta)
    y_l = yc - yt*np.cos(theta)

    x_u[-1] = x_l[-1]
    y_u[-1] = y_l[-1]

    return np.array([x_u,y_u,x_l,y_l])


def transform(points,chord,AoA,sweep,dihedral,distance,apply_trans,trans):

    # Scale Airfoil by Desired Chord Length

    points = chord*points

    # Rotate Airfoil by Angle of Attack


        # Enable these lines if wanting to rotate about point other than origin
        #xloc = 0.25
        #points[0,:] = points[0,:] - xloc
        #points[2,:] = points[2,:] - xloc

    theta = np.deg2rad(AoA)

    R = np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])

    upper = np.dot(R,points[0:2,:])
    lower = np.dot(R,points[2:,:])

        #upper[0] = upper[0] + 0.25
        #lower[0] = lower[0] + 0.25

    # Translate Section to Accommodate Previous section translations

    if apply_trans:
        upper[0, :] = upper[0, :] + trans[0]
        lower[0, :] = lower[0, :] + trans[0]

        upper[1, :] = upper[1, :] + trans[1]
        lower[1, :] = lower[1, :] + trans[1]


    # Translate Airofil According to Sweep Angle (modifies x coordinates only)

    x_trans = distance*np.sin(np.deg2rad(sweep))

    upper[0,:] = upper[0,:] + x_trans
    lower[0,:] = lower[0,:] + x_trans

    # Translate Airfoil According to Dihedral Angle (modifies y coordinates only)

    y_trans = distance*np.sin(np.deg2rad(dihedral))

    upper[1,:] = upper[1,:] + y_trans
    lower[1, :] = lower[1, :] + y_trans

    # Apply Z Coordinate for points

    Z = [distance]*len(upper[0,:]) + trans[2] if apply_trans else [distance]*len(upper[0,:])

    # Format and Return

    points = np.array([upper[0],upper[1],lower[0],lower[1],Z])


    return points, np.array([x_trans,y_trans,distance])


def find_xy(p1, p2, z):
    x1, y1, z1 = p1
    x2, y2, z2 = p2

    if z2 < z1:
        return find_xy(p2, p1, z)

    x = np.interp(z, (z1, z2), (x1, x2))
    y = np.interp(z, (z1, z2), (y1, y2))

    return x, y


def point_generate(sections,divisions=10,top=True):

    if top:

        ind1 = 0
        ind2 = 1

    else:

        ind1 = 2
        ind2 = 3

    interp_points = list()

    # Generate Interpolated Points

    i = 0
    cnt = 0

    for j in range(0,len(sections[i][0,:])):

        z = np.linspace(sections[i][4,j],sections[i+1][4,j],divisions)

        x,y = find_xy([sections[i][ind1,j],sections[i][ind2,j],sections[i][4,j]],
                      [sections[i+1][ind1,j],sections[i+1][ind2,j],sections[i+1][4,j]],z)

        tmp = np.array([x,y,z]).T

        if cnt == 0:

            interp_points = tmp

        else:

            interp_points = np.concatenate((interp_points,tmp))

        cnt +=1

    return interp_points



def tri_generate(vertices,divisions):


    n = divisions - 1

    rows = divisions*np.arange(0,int(len(vertices)/divisions))[:-1]

    for i in rows:

        for j in range(i,i+n):

            tmp = [[j,j+divisions,j+1],[j+divisions,j+divisions+1,j+1]]

            if j == 0:

                faces = tmp

            else:

                faces = np.concatenate((faces,tmp))

    return faces

def generate_stl(vertices,faces,name):


    wing = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype),name=name)
    for i, f in enumerate(faces):
        for j in range(3):
            wing.vectors[i][j] = vertices[f[j], :]

    #wing.save(sys_path+filename,mode=stl.Mode.ASCII)

    return wing


def solid_names(names,file_path):

    sys_path = file_path

    for i in range(0,len(names)):

        name = names[i]

        filename = name + '.stl'

        with open(sys_path+filename,'r') as file:

            data = file.readlines()

        data[0] = 'solid '+name+'\n'
        data[-1] = 'endsolid '+name

        with open(sys_path+filename,'w') as file:

            file.writelines(data)



def generate_cap(airfoil):

    topx = airfoil[0][:-1]
    topy = airfoil[1][:-1]
    bottomx = airfoil[2][1:]
    bottomy = airfoil[3][1:]
    z = airfoil[4]

    concatx = [None]*(len(topx)+len(bottomx))
    concaty = [None] * (len(topx) + len(bottomx))

    concatx[1::2] = topx
    concatx[::2] = bottomx

    concaty[1::2] = topy
    concaty[::2] = bottomy

    concatz = [z[0]]*len(concatx)

    num_tri = len(concatx) - 2

    for i in range(0,num_tri,2):

        tmp = [[i+2,i+1,i],[i+2,i+3,i+1]]

        if i == 0:

            faces = tmp

        else:

            faces = np.concatenate((faces,tmp))

    return np.array([concatx,concaty,concatz]),faces

def refinement_regions(sections,spacing):

    section0 = sections[0]
    section1 = sections[1]
    section2 = sections[2]

    x0 = [min(np.concatenate((section0[0], section0[2]))) - spacing,
          max(np.concatenate((section0[0], section0[2]))) + spacing]
    y0 = [min(np.concatenate((section0[1], section0[3]))) - spacing,
          max(np.concatenate((section0[1], section0[3]))) + spacing]

    x1 = [min(np.concatenate((section1[0], section1[2]))) - spacing,
          max(np.concatenate((section1[0], section1[2]))) + spacing]
    y1 = [min(np.concatenate((section1[1], section1[3]))) - spacing,
          max(np.concatenate((section1[1], section1[3]))) + spacing]

    x2 = [min(np.concatenate((section2[0], section2[2]))) - spacing,
          max(np.concatenate((section2[0], section2[2]))) + spacing]
    y2 = [min(np.concatenate((section2[1], section2[3]))) - spacing,
          max(np.concatenate((section2[1], section2[3]))) + spacing]

    z0 = 0
    z1 = section1[4][0]
    z2 = section2[4][0] + spacing

    #Vertices on Plane 1

    v1 = [x0[0],y0[0],z0]
    v2 = [x0[1],y0[0],z0]
    v3 = [x0[1],y0[1],z0]
    v4 = [x0[0],y0[1],z0]

    vert1 = np.array([v1,v2,v3,v4])

    #Vertices on Plane 2

    v1 = [x1[0], y1[0], z1]
    v2 = [x1[1], y1[0], z1]
    v3 = [x1[1], y1[1], z1]
    v4 = [x1[0], y1[1], z1]

    vert2 = np.array([v1,v2,v3,v4])

    # Vertices on Plane 3

    v1 = [x2[0], y2[0], z2]
    v2 = [x2[1], y2[0], z2]
    v3 = [x2[1], y2[1], z2]
    v4 = [x2[0], y2[1], z2]

    vert3 = np.array([v1, v2, v3, v4])

    # Define Vertices for 2 Wing Sections

    v_span1 = np.concatenate((vert1,vert2))
    v_span2 = np.concatenate((vert2,vert3))

    ind1 = np.array([[0,2,1],
                   [0,3,2],
                   [0,1,5],
                   [0,5,4],
                   [1,2,6],
                    [1,6,5],
                   [2,7,6],
                   [2,3,7],
                   [0,4,7],
                   [0,7,3]])

    ind2 = np.array([[0,1,5],
                   [0,5,4],
                   [1,2,6],
                    [1,6,5],
                   [2,7,6],
                   [2,3,7],
                   [0,4,7],
                   [0,7,3],
                   [5,6,4],
                   [4,6,7]])

    region1 = mesh.Mesh(np.zeros(ind1.shape[0], dtype=mesh.Mesh.dtype),name='refine1')
    for i, f in enumerate(ind1):
        for j in range(3):
            region1.vectors[i][j] = v_span1[f[j], :]

    region2 = mesh.Mesh(np.zeros(ind2.shape[0], dtype=mesh.Mesh.dtype), name='refine2')
    for i, f in enumerate(ind2):
        for j in range(3):
            region2.vectors[i][j] = v_span2[f[j], :]

    total_region = mesh.Mesh(np.concatenate([region1.data,region2.data]))

    return total_region



def wing_design_pipeline(dna,regionLevels,file_path,airfoil_res=25,span_res=10):

    NACA1 = dna[0]
    NACA2 = dna[1]
    NACA3 = dna[2]
    length = dna[3]
    chord = dna[4]
    percent_l = dna[5]
    sweep1 = dna[6]
    sweep2 = dna[7]
    dihed1 = dna[8]
    dihed2 = dna[9]
    taper1 = dna[10]
    taper2 = dna[11]
    AoA1 = dna[12]
    AoA2 = dna[13]
    AoA3 = dna[14]

    # Generate Airfoil Sections

    foil1 = naca_4digit(NACA1,airfoil_res)
    foil2 = naca_4digit(NACA2,airfoil_res)
    foil3 = naca_4digit(NACA3,airfoil_res)

    # Translate and Rotate Airfoil Sections

    section0,empty = transform(foil1,chord,AoA1,sweep1,dihed1,0,False,[])
    section1, trans1 = transform(foil2, taper1*chord, AoA2, sweep2,
                                 dihed1, percent_l*length, False, [])
    section2,trans2 = transform(foil3,chord*taper1*taper2,AoA3,sweep2,
                                dihed2,(1-percent_l)*length,True,trans1)

    # Prepare Data for Mesh Generation

    set1 = [section0,section1]
    set2 = [section1,section2]

    # Generate Mesh

    vertices1 = np.array(point_generate(set1, span_res))
    faces1 = np.array(tri_generate(vertices1, span_res))

    vertices2 = np.array(point_generate(set2, span_res))
    faces2 = np.array(tri_generate(vertices2, span_res))

    vertices3 = np.array(point_generate(set1, span_res, False))
    faces3 = np.array(tri_generate(vertices3, span_res))
    vertices4 = np.array(point_generate(set2, span_res, False))
    faces4 = np.array(tri_generate(vertices4, span_res))

    # Generate Cap Mesh

    cap_vertices1,cap_faces1 = generate_cap(section0)
    cap_vertices2,cap_faces2 = generate_cap(section2)

    # Generate STL Files from Triangulated Meshes

    names = ['upper1','upper2','lower1','lower2','cap1','cap2']
    mesh_filename = ['wing']
    regions = ['region1','region2','region3']

    meshes = []

    mesh1 = generate_stl(vertices1, faces1, names[0])
    mesh2 = generate_stl(vertices2, faces2, names[1])
    mesh3 = generate_stl(vertices3, faces3, names[2])
    mesh4 = generate_stl(vertices4, faces4, names[3])
    mesh5 = generate_stl(cap_vertices1.T,cap_faces1, names[4])
    mesh6 = generate_stl(cap_vertices2.T,cap_faces2, names[5])


    total_mesh = mesh.Mesh(np.concatenate([mesh1.data,mesh2.data,mesh3.data,mesh4.data,mesh5.data,mesh6.data]))

    sys_path = file_path
    filename = mesh_filename[0] + '.stl'

    total_mesh.save(sys_path+filename,mode=stl.Mode.ASCII)

    for i in range(0,len(regionLevels)):

        total_region1 = refinement_regions([section0,section1,section2],regionLevels[i])
        total_region1.save(sys_path+'region'+str(i+1)+'.stl',mode=stl.Mode.ASCII)

    # Rename Solids in First and Last Lines of STL files

    solid_names(mesh_filename,file_path)
    solid_names(regions,file_path)


#NACA1 = dna[0]
#NACA2 = dna[1]
#NACA3 = dna[2]
#length = dna[3]
#chord = dna[4]
#percent_l = dna[5]
#sweep1 = dna[6]
#sweep2 = dna[7]
#dihed1 = dna[8]
#dihed2 = dna[9]
#taper1 = dna[10]
#taper2 = dna[11]
#AoA1 = dna[12]
#AoA2 = dna[13]
#AoA3 = dna[14]


