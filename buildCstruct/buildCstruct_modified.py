'''
buildCstruct - Modified by: mdanh
Build armchair carbon nanotube structures with distance-based bond detection
'''

from sys import argv, exit
from optparse import OptionParser as OP
try:
    from numpy import zeros, pi, sin, cos, modf, sqrt
except:
    print("Numpy not installed. Exiting...")
    exit(10)
from os import path as path_os

def getdist(at1, at2):
    '''Calculate distance between two particles'''
    dist_at = sqrt((at2[0]-at1[0])**2+(at2[1]-at1[1])**2+(at2[2]-at1[2])**2)
    return dist_at

def filecheck(file):
    '''Check if file exists'''
    return path_os.isfile(file)

def parsecmd():
    description = "Build armchair carbon nanotube structures"
    usage = "usage: %prog [options] output_file"
    parser = OP(version='%prog 1.2-distance', description=description, usage=usage)
    
    parser.add_option('-g', '--geometry', dest='geometry', nargs=2, type='float',
                      help='define the geometry: index_n cnt_length')
    parser.add_option('-f', '--funct', dest='functionalization', default='none',
                      help='define functionalization: none, oh, cooh, coo')
    parser.add_option('--xyz', dest='xyz', action='store_true',
                      help='write xyz file')
    parser.add_option('--gro', dest='gro', action='store_true',
                      help='write gromacs gro file')
    parser.add_option('--mol2', dest='mol2', action='store_true',
                      help='write mol2 file')
    (options, args) = parser.parse_args(argv[1:])
    
    if len(args) == 0:
        parser.exit(parser.print_help())
    
    if len(args) > 1:
        parser.error('Too many arguments provided')
    
    return options, args

def armcnt(n, l, ccbond, funct):
    '''Build armchair carbon nanotube'''
    atc = []
    circ1 = []
    circ2 = []
    if funct == "oh":
        c1charge = -0.28
        c1acharge = 0.24
        c2charge = 0.01
    elif funct == "coo":
        c1charge = -0.34
        c1acharge = -0.09
        c2charge = 0.03
    elif funct == "cooh":
        c1charge = -0.12
        c1acharge = -0.1
        c2charge = 0.03
    else: 
        c1charge = -0.16
        c1acharge = -0.16
        c2charge = 0.03
    cmcharge = 0.0
    dx = ccbond*cos(120/2*(pi/180.0))
    dy = ccbond*sin(120/2*(pi/180.0))
    radius = (n*(2*dx+ccbond)+n*ccbond)/(2*pi)
    ycoord = +dy
    natoms = 2*n
    
    for i in range(n):
        circ1.append(2*dx+ccbond)
        circ1.append(ccbond)
        circ2.append(ccbond)
        circ2.append(2*dx+ccbond)
    
    circ1.insert(0, 0.0)
    circ1.pop()
    circ2.insert(0, dx)
    circ2.pop()
    
    while ycoord > -l:
        ycoord -= dy
        arc = 0.0
        if ycoord == 0:
            ccharge_circ1a = c1acharge
            ccharge_circ1 = c1charge
            ccharge_circ2a = c2charge
            ccharge_circ2 = c2charge
        elif ycoord < -l+dy:
            ccharge_circ1a = c2charge
            ccharge_circ1 = c2charge
            ccharge_circ2a = c1acharge
            ccharge_circ2 = c1charge
        else:
            ccharge_circ1 = cmcharge
            ccharge_circ2 = cmcharge   
            ccharge_circ1a = cmcharge
            ccharge_circ2a = cmcharge 
        
        for i in range(natoms):
            if modf(float(i)/2.0)[0] == 0:
                newccharge = ccharge_circ1a
            else:
                newccharge = ccharge_circ1
            tmpcoords = ['C']
            arc += circ1[i]
            theta = arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
        ycoord -= dy
        arc = 0.0
        for i in range(natoms):
            if modf(float(i)/2.0)[0] == 0:
                newccharge = ccharge_circ2a
            else:
                newccharge = ccharge_circ2
            tmpcoords = ['C']
            arc += circ2[i]
            theta = arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
    
    pbc_l = abs(ycoord)+dy
    print('\n*******************************')
    print('armchair CNT: n= ', n, ' l (ang)= ', abs(ycoord))
    print('periodicity (if apply) (ang)= ', pbc_l)
    print('diameter (ang): ', 2*radius)

    return atc, natoms, pbc_l, len(atc)

def add_COO(coords, natx, is_protonated):
    '''Add COO- groups and hydrogens to armchair CNT'''
    Hcov_r = 0.32 
    Ocov_r = 0.66
    Ccov_r = 0.77

    chbond = Hcov_r+Ccov_r
    cobond = Ccov_r+Ocov_r
    ohbond = Ocov_r+Hcov_r
    ccbond = 2*Ccov_r

    if is_protonated:
        hocharge = 0.44
        h1charge = 0.17
        c0charge = 0.7
        o1charge = -0.55
        o2charge = -0.6
    else:
        h1charge = 0.22
        c0charge = 0.83
        o1charge = -0.84
        o2charge = -0.84

    Hxz = chbond*cos(120/2*pi/180)
    Hy1 = chbond*sin(120/2*pi/180)
    Cxz = ccbond*cos(120/2*pi/180)
    Cy1 = ccbond*sin(120/2*pi/180)
        
    for i in range(natx):
        Hx1 = Hxz*cos(i*2*pi/natx)
        Hz1 = Hxz*sin(i*2*pi/natx)
        Cx1 = Cxz*cos(i*2*pi/natx)
        Cz1 = Cxz*sin(i*2*pi/natx)
        Ox1 = Cx1+cobond*cos(i*2*pi/natx - pi/4)
        Oz1 = Cz1+cobond*sin(i*2*pi/natx - pi/4)
        Oy1 = Cy1
        Ox2 = Cx1+cobond*cos(i*2*pi/natx + pi/4)
        Oz2 = Cz1+cobond*sin(i*2*pi/natx + pi/4)
        Oy2 = Cy1
        HOx = ohbond*cos(i*2*pi/natx+pi/4)
        HOz = ohbond*sin(i*2*pi/natx+pi/4)
        if modf(float(i)/2.0)[0] == 0:
            tmpcoords = ['C']
            tmpcoords.append(coords[i][1]+Cx1) 
            tmpcoords.append(coords[i][2]+Cy1)
            tmpcoords.append(coords[i][3]+Cz1)
            tmpcoords.append('C.2')
            tmpcoords.append(c0charge)
            coords.append(tmpcoords)
            tmpcoords = ['O']
            tmpcoords.append(coords[i][1]+Ox1) 
            tmpcoords.append(coords[i][2]+Oy1)
            tmpcoords.append(coords[i][3]+Oz1)
            tmpcoords.append('O.co2')
            tmpcoords.append(o1charge)
            coords.append(tmpcoords)
            tmpcoords = ['O']
            tmpcoords.append(coords[i][1]+Ox2) 
            tmpcoords.append(coords[i][2]+Oy2)
            tmpcoords.append(coords[i][3]+Oz2)
            tmpcoords.append('O.co2')
            tmpcoords.append(o2charge)
            coords.append(tmpcoords)
            if is_protonated:
                tmpcoords = ['H']
                tmpcoords.append(coords[i][1]+Ox2+HOx) 
                tmpcoords.append(coords[i][2]+Oy2)
                tmpcoords.append(coords[i][3]+Oz2+HOz)
                tmpcoords.append('H')
                tmpcoords.append(hocharge)
                coords.append(tmpcoords)
        else:
            tmpcoords = ['H']
            tmpcoords.append(coords[i][1])
            tmpcoords.append(coords[i][2]+chbond)
            tmpcoords.append(coords[i][3])
            tmpcoords.append('H')
            tmpcoords.append(h1charge)
            coords.append(tmpcoords)

    if is_protonated:
        added_atoms = int(2.5*natx)
    else:
        added_atoms = int(2*natx)

    for i in range(len(coords)-added_atoms-natx, len(coords)-added_atoms):
        Hx1 = Hxz*cos(i*2*pi/natx)
        Hz1 = Hxz*sin(i*2*pi/natx)
        Cx1 = Cxz*cos(i*2*pi/natx)
        Cz1 = Cxz*sin(i*2*pi/natx)
        Ox1 = Cx1+cobond*cos(i*2*pi/natx - pi/4)
        Oz1 = Cz1+cobond*sin(i*2*pi/natx - pi/4)
        Oy1 = Cy1
        Ox2 = Cx1+cobond*cos(i*2*pi/natx + pi/4)
        Oz2 = Cz1+cobond*sin(i*2*pi/natx + pi/4)
        Oy2 = Cy1
        HOx = ohbond*cos(i*2*pi/natx+pi/4)
        HOz = ohbond*sin(i*2*pi/natx+pi/4)
        if modf(float(i)/2.0)[0] == 0:
            tmpcoords = ['C']
            tmpcoords.append(coords[i][1]+Cx1)
            tmpcoords.append(coords[i][2]-Cy1)
            tmpcoords.append(coords[i][3]+Cz1)
            tmpcoords.append('C.2')
            tmpcoords.append(c0charge)
            coords.append(tmpcoords)
            tmpcoords = ['O']
            tmpcoords.append(coords[i][1]+Ox1) 
            tmpcoords.append(coords[i][2]-Oy1)
            tmpcoords.append(coords[i][3]+Oz1)
            tmpcoords.append('O.co2')
            tmpcoords.append(o1charge)
            coords.append(tmpcoords)
            tmpcoords = ['O']
            tmpcoords.append(coords[i][1]+Ox2) 
            tmpcoords.append(coords[i][2]-Oy2)
            tmpcoords.append(coords[i][3]+Oz2)
            tmpcoords.append('O.co2')
            tmpcoords.append(o2charge)
            coords.append(tmpcoords)
            if is_protonated:
                tmpcoords = ['H']
                tmpcoords.append(coords[i][1]+Ox2+HOx) 
                tmpcoords.append(coords[i][2]-Oy2)
                tmpcoords.append(coords[i][3]+Oz2+HOz)
                tmpcoords.append('H')
                tmpcoords.append(hocharge)
                coords.append(tmpcoords)
        else:
            tmpcoords = ['H']
            tmpcoords.append(coords[i][1])
            tmpcoords.append(coords[i][2]-chbond)
            tmpcoords.append(coords[i][3])
            tmpcoords.append('H')
            tmpcoords.append(h1charge)
            coords.append(tmpcoords)

def add_H(coords, natx, funct_OH):
    '''Add hydrogens to armchair CNT'''
    Hcov_r = 0.32 
    Ocov_r = 0.66
    Ccov_r = 0.77
    
    chbond = Hcov_r+Ccov_r
    cobond = Ccov_r+Ocov_r
    ohbond = Ocov_r+Hcov_r

    Hxz = chbond*cos(120/2*pi/180)
    Hy1 = chbond*sin(120/2*pi/180)
    Oxz = cobond*cos(120/2*pi/180)
    Oy1 = cobond*sin(120/2*pi/180)

    if funct_OH:
        h1charge = 0.18
        o1charge = -0.53
        h2charge = 0.37
    else:
        h1charge = 0.13

    for i in range(natx):
        Hx1 = Hxz*cos(i*2*pi/natx)
        Hz1 = Hxz*sin(i*2*pi/natx)
        Ox1 = Oxz*cos(i*2*pi/natx)
        Oz1 = Oxz*sin(i*2*pi/natx)
        Hx2 = Ox1+ohbond*cos(i*2*pi/natx)
        Hz2 = Oz1+ohbond*sin(i*2*pi/natx)
        Hy2 = Oy1

        if modf(float(i)/2.0)[0] == 0 and funct_OH:
            tmpcoords = ['O']
            tmpcoords.append(coords[i][1]+Ox1) 
            tmpcoords.append(coords[i][2]+Oy1)
            tmpcoords.append(coords[i][3]+Oz1)
            tmpcoords.append('O.3')
            tmpcoords.append(o1charge)
            coords.append(tmpcoords)
            tmpcoords = ['H']
            tmpcoords.append(coords[i][1]+Hx2) 
            tmpcoords.append(coords[i][2]+Hy2)
            tmpcoords.append(coords[i][3]+Hz2)
            tmpcoords.append('H')
            tmpcoords.append(h2charge)
            coords.append(tmpcoords)
        else:
            tmpcoords = ['H']
            tmpcoords.append(coords[i][1]+Hx1) 
            tmpcoords.append(coords[i][2]+Hy1)
            tmpcoords.append(coords[i][3]+Hz1)
            tmpcoords.append('H')
            tmpcoords.append(h1charge)
            coords.append(tmpcoords)

    if funct_OH:
        added_atoms = int(1.5*natx)
    else:
        added_atoms = natx

    for i in range(len(coords)-added_atoms-natx, len(coords)-added_atoms):
        Hx1 = Hxz*cos(i*2*pi/natx)
        Hz1 = Hxz*sin(i*2*pi/natx)
        Ox1 = Oxz*cos(i*2*pi/natx)
        Oz1 = Oxz*sin(i*2*pi/natx)
        Hx2 = Ox1+ohbond*cos(i*2*pi/natx)
        Hz2 = Oz1+ohbond*sin(i*2*pi/natx)
        if modf(float(i)/2.0)[0] == 0 and funct_OH:
            tmpcoords = ['O']
            tmpcoords.append(coords[i][1]+Ox1)
            tmpcoords.append(coords[i][2]-Oy1)
            tmpcoords.append(coords[i][3]+Oz1)
            tmpcoords.append('O.3')
            tmpcoords.append(o1charge)
            coords.append(tmpcoords)
            tmpcoords = ['H']
            tmpcoords.append(coords[i][1]+Hx2)
            tmpcoords.append(coords[i][2]-Oy1)
            tmpcoords.append(coords[i][3]+Hz2)
            tmpcoords.append('H')
            tmpcoords.append(h2charge)
            coords.append(tmpcoords)
        else:
            tmpcoords = ['H']
            tmpcoords.append(coords[i][1]+Hx1)
            tmpcoords.append(coords[i][2]-Hy1)
            tmpcoords.append(coords[i][3]+Hz1)
            tmpcoords.append('H')
            tmpcoords.append(h1charge)
            coords.append(tmpcoords)

def connect_distance(coords, natx, nohcoords):
    '''Build connectivity using distance-based method with atom type checking'''
    Ccov_r = 0.77
    Hcov_r = 0.32
    Ocov_r = 0.66

    btollcc = (2*Ccov_r)*15/100
    btollch = (Hcov_r+Ccov_r)*5/100
    btolloh = (Hcov_r+Ocov_r)*5/100
    btollco = (Ccov_r+Ocov_r)*5/100
    bondcc = [2*Ccov_r-btollcc, 2*Ccov_r+btollcc]
    bondch = [(Hcov_r+Ccov_r)-btollch, (Hcov_r+Ccov_r)+btollch]
    bondco = [(Ocov_r+Ccov_r)-btollco, (Ocov_r+Ccov_r)+btollco]
    bondoh = [(Hcov_r+Ocov_r)-btolloh, (Hcov_r+Ocov_r)+btolloh]
    connect = zeros((len(coords), 3), int)
    bondlist = []
    bondnumber = 0
    
    for i in range(len(coords)):
        for j in range(i+1, i+2*natx):
            if j < nohcoords:
                at1 = [coords[i][1], coords[i][2], coords[i][3]]
                at2 = [coords[j][1], coords[j][2], coords[j][3]]
                bond = getdist(at1, at2)
                if bond >= bondcc[0] and bond < bondcc[1]:
                    if i < j:
                        bondnumber += 1
                        tmpbond = [bondnumber, i+1, j+1, 'ar']
                        bondlist.append(tmpbond)
                    for k in range(3):
                        if connect[i][k] == 0:
                            connect[i][k] = j+1
                            break
                    for k in range(3):
                        if connect[j][k] == 0:
                            connect[j][k] = i+1
                            break
        if len(coords) != nohcoords:
            for j in range(nohcoords, len(coords)):
                at1 = [coords[i][1], coords[i][2], coords[i][3]]
                at2 = [coords[j][1], coords[j][2], coords[j][3]]
                bond = getdist(at1, at2)
                
                atom1_type = coords[i][0]
                atom2_type = coords[j][0]
                
                valid_bond = False
                
                if ((atom1_type == 'C' and atom2_type == 'H') or 
                    (atom1_type == 'H' and atom2_type == 'C')) and \
                   (bond >= bondch[0] and bond < bondch[1]):
                    valid_bond = True
                
                elif ((atom1_type == 'C' and atom2_type == 'O') or 
                      (atom1_type == 'O' and atom2_type == 'C')) and \
                     (bond >= bondco[0] and bond < bondco[1]):
                    valid_bond = True
                
                elif ((atom1_type == 'O' and atom2_type == 'H') or 
                      (atom1_type == 'H' and atom2_type == 'O')) and \
                     (bond >= bondoh[0] and bond < bondoh[1]):
                    valid_bond = True
                
                elif (atom1_type == 'C' and atom2_type == 'C') and \
                     (bond >= bondcc[0] and bond < bondcc[1]):
                    valid_bond = True
                
                if valid_bond:
                    if i < j:
                        bondnumber += 1
                        tmpbond = [bondnumber, i+1, j+1, '1']
                        bondlist.append(tmpbond)
                    for k in range(3):
                        if connect[i][k] == 0:
                            connect[i][k] = j+1
                            break
                    for k in range(3):
                        if connect[j][k] == 0:
                            connect[j][k] = i+1
                            break
    
    return connect, bondlist

def write_xyz(file, data):
    '''Write xyz file'''
    file.write(" "+str(len(data))+"\nGenerated by buildCstruct - Armchair CNT\n")
    for line in data:
       outline = "%-3s%12.6f%12.6f%12.6f" % (line[0], float(line[1]), float(line[2]), float(line[3]))
       file.write(outline+"\n")

def write_gro(file, data, pbc1=""):
    '''Write gromacs gro file'''
    file.write("Generated by buildCstruct - Armchair CNT\n "+str(len(data))+"\n")
    for index, line in enumerate(data):
       outline = "%5i%-5s%5s%5i%8.3f%8.3f%8.3f" % (1, "CNT1", line[0], index, float(line[1])/10.0, float(line[2])/10.0, float(line[3])/10.0)
       file.write(outline+"\n")
    if pbc1 == "":
        outline = "  10   10   10\n"
    else:
        outline = "  10   "+str(float(pbc1)/10.0)+"   10\n"
    file.write(outline+"\n")

def write_mol2(file, data, bondlist):
    '''Write mol2 file'''
    file.write("@<TRIPOS>MOLECULE\nCNT_Distance\n "+str(len(data))+" "+str(len(bondlist))+" 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
    for index, line in enumerate(data):
       outline = "%7i %5s %8.3f %8.3f %8.3f %7s %7i %7s %8.3f" % (index+1, line[0], float(line[1]), float(line[2]), float(line[3]), line[4], 1, "CNT1", float(line[5]))
       file.write(outline+"\n")

    file.write("@<TRIPOS>BOND\n")
    for line in bondlist:
       outline = "%7i %7i %7i %7s" % (line[0], line[1], line[2], line[3])
       file.write(outline+"\n")

def main():
    '''Main function'''
    ccbond = 1.3874
    
    (options, args) = parsecmd()
    ofile = args[0]
    funct = options.functionalization.lower()

    coords, natx, pbc_l, nohcoords = armcnt(int(options.geometry[0]), float(options.geometry[1]), ccbond, funct)
    
    if funct == "coo":
        add_COO(coords, natx, False)
    elif funct == "cooh":
        add_COO(coords, natx, True)
    elif funct == "oh":
        add_H(coords, natx, True)        
    else:
        add_H(coords, natx, False)
    
    print('Atoms: ', len(coords))
    
    if not (options.xyz or options.gro):
        conn, bondlist = connect_distance(coords, natx, nohcoords)
    
    print('*******************************')
    OUT = open(ofile, 'w')
    
    if options.xyz:
        write_xyz(OUT, coords)
    elif options.gro:
        write_gro(OUT, coords, pbc_l)
    else:
        if 'bondlist' not in locals():
            conn, bondlist = connect_distance(coords, natx, nohcoords)
        write_mol2(OUT, coords, bondlist)
   
    OUT.close()
    exit(0)
    
if __name__ == "__main__":
    main() 