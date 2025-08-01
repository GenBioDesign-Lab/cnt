#!/usr/bin/env python3
'''
Build armchair carbon nanotube structures with -OH functionalization
Output: mol2 format only
'''

from sys import argv, exit
from optparse import OptionParser as OP
try:
    from numpy import zeros, pi, sin, cos, modf, sqrt
except:
    print("Numpy not installed. Exiting...")
    exit(10)

def getdist(at1, at2):
    '''Calculate distance between two particles'''
    dist_at = sqrt((at2[0]-at1[0])**2+(at2[1]-at1[1])**2+(at2[2]-at1[2])**2)
    return dist_at

def parsecmd():
    description = "Build armchair carbon nanotube structures with OH functionalization"
    usage = "usage: %prog [options] output_file.mol2"
    parser = OP(version='%prog 1.0-OH', description=description, usage=usage)
    
    parser.add_option('-n', '--index', dest='n_index', type='int', default=8,
                      help='nanotube index n (default: 8)')
    parser.add_option('-l', '--length', dest='length', type='float', default=20.0,
                      help='nanotube length in Angstroms (default: 20.0)')
    
    (options, args) = parser.parse_args(argv[1:])
    
    if len(args) == 0:
        parser.exit(parser.print_help())
    
    if len(args) > 1:
        parser.error('Too many arguments provided')
    
    # Ensure output file has .mol2 extension
    output_file = args[0]
    if not output_file.endswith('.mol2'):
        output_file += '.mol2'
    
    return options, output_file

def create_oh_cnt(n, l, ccbond):
    '''Build armchair carbon nanotube with OH functionalization charges'''
    atc = []
    circ1 = []
    circ2 = []
    
    # Charges for OH functionalization
    c1charge = -0.28    # End carbon with OH
    c1acharge = 0.24    # Adjacent to end carbon
    c2charge = 0.01     # Regular carbon
    cmcharge = 0.0      # Middle carbons
    
    dx = ccbond * cos(120/2 * (pi/180.0))
    dy = ccbond * sin(120/2 * (pi/180.0))
    radius = (n * (2*dx + ccbond) + n * ccbond) / (2*pi)
    ycoord = +dy
    natoms = 2 * n
    
    # Build the circular patterns
    for i in range(n):
        circ1.append(2*dx + ccbond)
        circ1.append(ccbond)
        circ2.append(ccbond)
        circ2.append(2*dx + ccbond)
    
    circ1.insert(0, 0.0)
    circ1.pop()
    circ2.insert(0, dx)
    circ2.pop()
    
    # Build the nanotube structure
    while ycoord > -l:
        ycoord -= dy
        arc = 0.0
        
        # Assign charges based on position
        if ycoord == 0:  # Top end
            ccharge_circ1a = c1acharge
            ccharge_circ1 = c1charge
            ccharge_circ2a = c2charge
            ccharge_circ2 = c2charge
        elif ycoord < -l + dy:  # Bottom end
            ccharge_circ1a = c2charge
            ccharge_circ1 = c2charge
            ccharge_circ2a = c1acharge
            ccharge_circ2 = c1charge
        else:  # Middle section
            ccharge_circ1 = cmcharge
            ccharge_circ2 = cmcharge   
            ccharge_circ1a = cmcharge
            ccharge_circ2a = cmcharge 
        
        # First ring
        for i in range(natoms):
            if modf(float(i)/2.0)[0] == 0:
                newccharge = ccharge_circ1a
            else:
                newccharge = ccharge_circ1
            tmpcoords = ['C']
            arc += circ1[i]
            theta = arc/radius
            tmpcoords.append(radius * cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius * sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
        
        ycoord -= dy
        arc = 0.0
        
        # Second ring
        for i in range(natoms):
            if modf(float(i)/2.0)[0] == 0:
                newccharge = ccharge_circ2a
            else:
                newccharge = ccharge_circ2
            tmpcoords = ['C']
            arc += circ2[i]
            theta = arc/radius
            tmpcoords.append(radius * cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius * sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
    
    pbc_l = abs(ycoord) + dy
    print('\n*******************************')
    print('OH-functionalized armchair CNT:')
    print('  n =', n)
    print('  length (Å) =', abs(ycoord))
    print('  diameter (Å) =', 2*radius)
    print('  periodicity (Å) =', pbc_l)
    print('*******************************')

    return atc, natoms, pbc_l, len(atc)

def add_oh_groups(coords, natx):
    '''Add OH groups to ALL carbons at both ends of armchair CNT'''
    Hcov_r = 0.32 
    Ocov_r = 0.66
    Ccov_r = 0.77
    
    cobond = Ccov_r + Ocov_r
    ohbond = Ocov_r + Hcov_r

    Oxz = cobond * cos(120/2 * pi/180)
    Oy1 = cobond * sin(120/2 * pi/180)

    # Charges for OH functionalization
    o1charge = -0.53  # O in OH group
    h2charge = 0.37   # H in OH group

    # Add OH groups to ALL carbons at top end
    for i in range(natx):
        Ox1 = Oxz * cos(i * 2 * pi / natx)
        Oz1 = Oxz * sin(i * 2 * pi / natx)
        Hx2 = Ox1 + ohbond * cos(i * 2 * pi / natx)
        Hz2 = Oz1 + ohbond * sin(i * 2 * pi / natx)
        Hy2 = Oy1

        # Add oxygen
        tmpcoords = ['O']
        tmpcoords.append(coords[i][1] + Ox1) 
        tmpcoords.append(coords[i][2] + Oy1)
        tmpcoords.append(coords[i][3] + Oz1)
        tmpcoords.append('O.3')
        tmpcoords.append(o1charge)
        coords.append(tmpcoords)
        
        # Add hydrogen on oxygen
        tmpcoords = ['H']
        tmpcoords.append(coords[i][1] + Hx2) 
        tmpcoords.append(coords[i][2] + Hy2)
        tmpcoords.append(coords[i][3] + Hz2)
        tmpcoords.append('H')
        tmpcoords.append(h2charge)
        coords.append(tmpcoords)

    added_atoms = 2 * natx  # 2 atoms (O+H) per carbon

    # Add OH groups to ALL carbons at bottom end
    for i in range(len(coords) - added_atoms - natx, len(coords) - added_atoms):
        Ox1 = Oxz * cos(i * 2 * pi / natx)
        Oz1 = Oxz * sin(i * 2 * pi / natx)
        Hx2 = Ox1 + ohbond * cos(i * 2 * pi / natx)
        Hz2 = Oz1 + ohbond * sin(i * 2 * pi / natx)
        
        # Add oxygen
        tmpcoords = ['O']
        tmpcoords.append(coords[i][1] + Ox1)
        tmpcoords.append(coords[i][2] - Oy1)
        tmpcoords.append(coords[i][3] + Oz1)
        tmpcoords.append('O.3')
        tmpcoords.append(o1charge)
        coords.append(tmpcoords)
        
        # Add hydrogen on oxygen
        tmpcoords = ['H']
        tmpcoords.append(coords[i][1] + Hx2)
        tmpcoords.append(coords[i][2] - Oy1)
        tmpcoords.append(coords[i][3] + Hz2)
        tmpcoords.append('H')
        tmpcoords.append(h2charge)
        coords.append(tmpcoords)

def create_bonds(coords, natx, nohcoords):
    '''Build connectivity using distance-based method'''
    Ccov_r = 0.77
    Hcov_r = 0.32
    Ocov_r = 0.66

    # Bond tolerances (5-15% tolerance)
    btollcc = (2 * Ccov_r) * 15/100
    btollch = (Hcov_r + Ccov_r) * 5/100
    btolloh = (Hcov_r + Ocov_r) * 5/100
    btollco = (Ccov_r + Ocov_r) * 5/100
    
    # Bond distance ranges
    bondcc = [2*Ccov_r - btollcc, 2*Ccov_r + btollcc]
    bondch = [(Hcov_r + Ccov_r) - btollch, (Hcov_r + Ccov_r) + btollch]
    bondco = [(Ocov_r + Ccov_r) - btollco, (Ocov_r + Ccov_r) + btollco]
    bondoh = [(Hcov_r + Ocov_r) - btolloh, (Hcov_r + Ocov_r) + btolloh]
    
    connect = zeros((len(coords), 3), int)
    bondlist = []
    bondnumber = 0
    
    # Create bonds between carbon atoms (CNT structure)
    for i in range(len(coords)):
        for j in range(i+1, i + 2*natx):
            if j < nohcoords:
                at1 = [coords[i][1], coords[i][2], coords[i][3]]
                at2 = [coords[j][1], coords[j][2], coords[j][3]]
                bond = getdist(at1, at2)
                if bond >= bondcc[0] and bond < bondcc[1]:
                    if i < j:
                        bondnumber += 1
                        tmpbond = [bondnumber, i+1, j+1, 'ar']
                        bondlist.append(tmpbond)
                    # Update connectivity
                    for k in range(3):
                        if connect[i][k] == 0:
                            connect[i][k] = j+1
                            break
                    for k in range(3):
                        if connect[j][k] == 0:
                            connect[j][k] = i+1
                            break
        
        # Create bonds with functional groups
        if len(coords) != nohcoords:
            for j in range(nohcoords, len(coords)):
                at1 = [coords[i][1], coords[i][2], coords[i][3]]
                at2 = [coords[j][1], coords[j][2], coords[j][3]]
                bond = getdist(at1, at2)
                
                atom1_type = coords[i][0]
                atom2_type = coords[j][0]
                
                valid_bond = False
                
                # Check for valid bond types
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
                    # Update connectivity
                    for k in range(3):
                        if connect[i][k] == 0:
                            connect[i][k] = j+1
                            break
                    for k in range(3):
                        if connect[j][k] == 0:
                            connect[j][k] = i+1
                            break
    
    return connect, bondlist

def write_mol2(file, data, bondlist):
    '''Write mol2 file for OH-functionalized CNT'''
    file.write("@<TRIPOS>MOLECULE\n")
    file.write("OH_CNT\n")
    file.write(" %d %d 0 0 0\n" % (len(data), len(bondlist)))
    file.write("SMALL\n")
    file.write("USER_CHARGES\n\n")
    file.write("@<TRIPOS>ATOM\n")
    
    for index, line in enumerate(data):
        outline = "%7d %5s %8.3f %8.3f %8.3f %7s %7d %7s %8.3f\n" % (
            index+1, line[0], float(line[1]), float(line[2]), float(line[3]), 
            line[4], 1, "CNT1", float(line[5]))
        file.write(outline)

    file.write("@<TRIPOS>BOND\n")
    for line in bondlist:
        outline = "%7d %7d %7d %7s\n" % (line[0], line[1], line[2], line[3])
        file.write(outline)

def main():
    ccbond = 1.3874  # C-C bond length in CNT
    
    (options, output_file) = parsecmd()
    
    # Create CNT structure with OH functionalization charges
    coords, natx, pbc_l, nohcoords = create_oh_cnt(options.n_index, options.length, ccbond)
    
    # Add OH groups to all carbon atoms
    add_oh_groups(coords, natx)
    
    print('Total atoms:', len(coords))
    
    # Create bond connectivity
    conn, bondlist = create_bonds(coords, natx, nohcoords)
    print('Total bonds:', len(bondlist))
    
    # Write mol2 file
    with open(output_file, 'w') as outfile:
        write_mol2(outfile, coords, bondlist)
    
    print('\nOH-functionalized CNT saved to:', output_file)

if __name__ == "__main__":
    main()