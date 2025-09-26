#!/usr/bin/env python3
'''
Build armchair carbon nanotube structures with amide (-CONH₂) functionalization
The -l parameter specifies the CNT backbone length BEFORE adding amide groups
Output: mol2 format
'''

from sys import argv, exit
import argparse
from numpy import zeros, pi, sin, cos, modf, sqrt, array, cross, linalg

# ================================================================================
# CHEMICAL PARAMETERS
# ================================================================================

# CNT structural parameters
CNT_PARAMS = {
    'cc_bond_length': 1.3874,  # C-C bond length in CNT backbone (Å)
}

# Amide functionalization parameters
AMIDE_PARAMS = {
    'bond_lengths': {
        'c_cnt_c_amide': 1.47,    # C(CNT)-C(amide) single bond (Å)
        'c_o_double': 1.243,      # C=O double bond in amide (Å)  
        'c_n_amide': 1.325,       # C-N amide bond (Å)
        'n_h': 1.00,              # N-H bond in amide (Å)
        'c_h_passivate': 1.087,   # C-H bond for passivating hydrogens (Å)
    },
    'charges': {
        'c_amide': 0.55,          # Carbonyl carbon charge
        'o_carbonyl': -0.55,      # Carbonyl oxygen charge
        'n_amide': -0.64,         # Amide nitrogen charge
        'h_amide': 0.32,          # Amide hydrogen charge
        # Charges for passivating C-H bonds
        'c_passivate': -0.115,    # Carbon bonded to passivating hydrogen
        'h_passivate': 0.115,     # Passivating hydrogen charge
    },
    'geometry': {
        'default_tilt_angle': 60.0,  # Default tilt angle for amide bond axis (degrees)
        'amide_plane_angles': {
            'c_o_n': 120.0,          # O=C-N angle in amide plane (degrees)
            'n_h_h': 120.0,          # H-N-H angle (degrees)
        }
    }
}

# Bond detection tolerances for connectivity analysis
BOND_TOLERANCES = {
    'cc_aromatic': 0.15,      # 15% tolerance for aromatic C-C bonds
    'ch_single': 0.05,        # 5% tolerance for C-H bonds
    'co_single': 0.05,        # 5% tolerance for C-O single bonds
    'oh_single': 0.05,        # 5% tolerance for O-H bonds
    'cn_amide': 0.05,         # 5% tolerance for C-N amide bonds
    'nh_amide': 0.05,         # 5% tolerance for N-H bonds
    'co_double': 0.04,        # 4% tolerance for C=O double bonds
    'cc_amide': 0.07,         # 7% tolerance for C(CNT)-C(amide) bonds
}

# Covalent radii for bond detection (Å)
COVALENT_RADII = {
    'C': 0.77,
    'H': 0.32,
    'O': 0.66,
    'N': 0.70,
}

def getdist(at1, at2):
    '''Calculate distance between two particles'''
    dist_at = sqrt((at2[0]-at1[0])**2+(at2[1]-at1[1])**2+(at2[2]-at1[2])**2)
    return dist_at

def validate_amide_geometry(coords, amide_h_indices):
    '''Validate that amide hydrogens have proper geometry'''
    validation_results = {
        'total_amide_h': len(amide_h_indices),
        'min_h_c_distance': float('inf'),
        'problematic_h_indices': []
    }
    
    for h_idx in amide_h_indices:
        h_pos = array([coords[h_idx][1], coords[h_idx][2], coords[h_idx][3]])
        
        # Check distances to all carbon atoms
        for i, atom in enumerate(coords):
            if atom[0] == 'C':
                c_pos = array([atom[1], atom[2], atom[3]])
                distance = linalg.norm(h_pos - c_pos)
                validation_results['min_h_c_distance'] = min(validation_results['min_h_c_distance'], distance)
                
                if distance < 1.3:  # Too close to carbon
                    validation_results['problematic_h_indices'].append((h_idx, i, distance))
    
    return validation_results

def validate_amide_bonding(coords, bondlist, amide_h_indices):
    '''Validate that amide groups have correct bonding patterns'''
    results = {
        'cn_bonds': 0,
        'co_double_bonds': 0, 
        'nh_bonds': 0,
        'spurious_ch_bonds': 0
    }
    
    # Analyze each bond
    for bond in bondlist:
        bond_id, atom1_idx, atom2_idx, bond_type = bond[0], bond[1]-1, bond[2]-1, bond[3]
        
        atom1_type = coords[atom1_idx][0]
        atom2_type = coords[atom2_idx][0]
        
        # Count C-N bonds
        if ((atom1_type == 'C' and atom2_type == 'N') or 
            (atom1_type == 'N' and atom2_type == 'C')) and bond_type == '1':
            results['cn_bonds'] += 1
            
        # Count C=O double bonds  
        elif ((atom1_type == 'C' and atom2_type == 'O') or 
              (atom1_type == 'O' and atom2_type == 'C')) and bond_type == '2':
            results['co_double_bonds'] += 1
            
        # Count N-H bonds
        elif ((atom1_type == 'N' and atom2_type == 'H') or 
              (atom1_type == 'H' and atom2_type == 'N')) and bond_type == '1':
            results['nh_bonds'] += 1
            
        # Check for spurious C-H bonds with amide hydrogens
        elif ((atom1_type == 'C' and atom2_type == 'H') or 
              (atom1_type == 'H' and atom2_type == 'C')) and bond_type == '1':
            h_idx = atom2_idx if atom2_type == 'H' else atom1_idx
            if h_idx in amide_h_indices:
                results['spurious_ch_bonds'] += 1
    
    return results

def parse_arguments():
    """Parse command line arguments using modern argparse."""
    parser = argparse.ArgumentParser(
        description="Build armchair carbon nanotube structures with amide (-CONH₂) functionalization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -n 8 -l 20.0 output.mol2
  %(prog)s -n 6 -l 15.0 --amide-angle 45.0 --no-stagger cnt_amide.mol2
          """)
    
    # Required positional argument
    parser.add_argument('output_file', 
                       help='Output file name (will add .mol2 extension if not present)')
    
    # CNT structure parameters
    parser.add_argument('-n', '--index', 
                       type=int, 
                       default=8,
                       help='Nanotube index n (default: %(default)s)')
    
    parser.add_argument('-l', '--length', 
                       type=float, 
                       default=20.0,
                       help='CNT backbone length in Angstroms BEFORE adding amide groups (default: %(default)s)')
    
    # Amide functionalization parameters
    parser.add_argument('--amide-angle', 
                       type=float, 
                       default=AMIDE_PARAMS['geometry']['default_tilt_angle'],
                       help='Tilt angle for amide bond axis in degrees (default: %(default)s)')
    
    parser.add_argument('--no-stagger', 
                       action='store_false', 
                       dest='stagger',
                       help='Disable staggered amide plane arrangement (default: staggered)')

    args = parser.parse_args()
    
    # Ensure output file has .mol2 extension
    if not args.output_file.endswith('.mol2'):
        args.output_file += '.mol2'
    
    # Validate arguments
    if args.index < 1:
        parser.error('Nanotube index n must be positive')
    
    if args.length <= 0:
        parser.error('CNT length must be positive')
    
    if not (0 <= args.amide_angle <= 90):
        parser.error('Amide angle must be between 0 and 90 degrees')
    
    return args

def create_cnt_backbone(n, backbone_length, ccbond):
    '''Build armchair carbon nanotube backbone structure'''
    atc = []
    circ1 = []
    circ2 = []
    
    # Simplified charges for CNT backbone (functionalization will add specific charges)
    c1charge = 0.0      # End carbon
    c1acharge = 0.0     # Adjacent to end carbon  
    c2charge = 0.0      # Regular carbon
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
    
    # Build the nanotube structure to achieve target backbone length
    while ycoord > -backbone_length:
        ycoord -= dy
        arc = 0.0
        
        # Assign charges based on position
        if ycoord == 0:  # Top end
            ccharge_circ1a = c1acharge
            ccharge_circ1 = c1charge
            ccharge_circ2a = c2charge
            ccharge_circ2 = c2charge
        elif ycoord < -backbone_length + dy:  # Bottom end
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
    
    actual_backbone_length = abs(ycoord)
    pbc_l = actual_backbone_length + dy
    
    return atc, natoms, pbc_l, len(atc), actual_backbone_length

def add_amide_groups(coords, natx, angle_deg=60.0, stagger=True):
    '''Add amide (-CONH₂) groups to every second carbon at both ends of armchair CNT
    
    Returns:
        tuple: (amide_extension, amide_h_indices) where amide_h_indices is a set of 
               indices for hydrogen atoms that belong to amide groups
    '''
    
    # Extract parameters from organized dictionaries
    bond_lengths = AMIDE_PARAMS['bond_lengths']
    charges = AMIDE_PARAMS['charges']
    
    # Bond lengths for amide group
    bond_c_cnt_c_amide = bond_lengths['c_cnt_c_amide']
    bond_c_o_double = bond_lengths['c_o_double']
    bond_c_n_amide = bond_lengths['c_n_amide']
    bond_n_h = bond_lengths['n_h']
    bond_c_h_passivate = bond_lengths['c_h_passivate']
    
    # Partial charges with descriptive names
    charge_C_amide = charges['c_amide']
    charge_O_carbonyl = charges['o_carbonyl']
    charge_N_amide = charges['n_amide']
    charge_H_amide = charges['h_amide']
    charge_C_passivate = charges['c_passivate']
    charge_H_passivate = charges['h_passivate']
    
    # Convert angle to radians
    theta = angle_deg * pi / 180.0
    
    # Calculate CNT radius for geometry
    ccbond = CNT_PARAMS['cc_bond_length']
    dx = ccbond * cos(120/2 * (pi/180.0))
    radius = (natx * (2*dx + ccbond) + natx * ccbond) / (2*pi)
    
    max_extension_top = 0.0
    max_extension_bottom = 0.0
    amide_h_indices = set()  # Track indices of hydrogen atoms in amide groups
    
    # Add functionalization to every second carbon at top end
    for i in range(natx):
        if modf(float(i)/2.0)[0] != 0:  # Unfunctionalized carbons get hydrogen atoms
            # Update carbon charge for C-H bonding
            coords[i][5] = charge_C_passivate
            
            # Add hydrogen atom
            c0_x, c0_y, c0_z = coords[i][1], coords[i][2], coords[i][3]
            
            # Calculate radial direction for H placement
            r_vec = array([c0_x, 0.0, c0_z])
            r_mag = linalg.norm(r_vec)
            if r_mag > 0:
                r_hat = r_vec / r_mag
            else:
                r_hat = array([1.0, 0.0, 0.0])
            
            # Place hydrogen radially outward with slight upward tilt
            y_hat = array([0.0, 1.0, 0.0])
            h_direction = 0.9 * r_hat + 0.1 * y_hat  # Mostly radial with slight upward tilt
            h_direction = h_direction / linalg.norm(h_direction)
            
            h_pos = array([c0_x, c0_y, c0_z]) + bond_c_h_passivate * h_direction
            
            tmpcoords = ['H'] + list(h_pos) + ['H', charge_H_passivate]
            coords.append(tmpcoords)
            continue  # Skip amide functionalization for this carbon
        
        # Get terminal carbon coordinates
        c0_x, c0_y, c0_z = coords[i][1], coords[i][2], coords[i][3]
        
        # Calculate local coordinate system
        # Radial unit vector (project to xz-plane and normalize)
        r_vec = array([c0_x, 0.0, c0_z])
        r_mag = linalg.norm(r_vec)
        if r_mag > 0:
            r_hat = r_vec / r_mag
        else:
            r_hat = array([1.0, 0.0, 0.0])
        
        # Axial unit vector (positive Y for top end)
        y_hat = array([0.0, 1.0, 0.0])
        
        # Tangent unit vector
        t_hat = cross(y_hat, r_hat)
        t_mag = linalg.norm(t_hat)
        if t_mag > 0:
            t_hat = t_hat / t_mag
        else:
            t_hat = array([0.0, 0.0, 1.0])
        
        # Apply staggering if enabled
        if stagger and i % 2 == 1:
            t_hat = -t_hat
        
        # Bond axis direction (tilted from radial)
        v_bond = cos(theta) * r_hat + sin(theta) * y_hat
        v_bond = v_bond / linalg.norm(v_bond)
        
        # Place amide carbon
        cam_pos = array([c0_x, c0_y, c0_z]) + bond_c_cnt_c_amide * v_bond
        
        # Define in-plane basis at amide carbon
        u_hat = -v_bond  # toward CNT
        w_hat = t_hat    # tangent direction
        
        # Place oxygen and nitrogen (120° separation in amide plane)
        angle_o = 120 * pi / 180
        angle_n = -120 * pi / 180
        
        dir_o = cos(angle_o) * u_hat + sin(angle_o) * w_hat
        dir_n = cos(angle_n) * u_hat + sin(angle_n) * w_hat
        
        o_pos = cam_pos + bond_c_o_double * dir_o
        n_pos = cam_pos + bond_c_n_amide * dir_n
        
        # Place hydrogens on nitrogen using robust orthonormal coordinate system
        # Build N-centric orthonormal frame to ensure proper H-N-H geometry
        y_axis_nh2 = (n_pos - cam_pos) / linalg.norm(n_pos - cam_pos)  # C→N direction
        
        # Orthonormalize tangent vector against C→N axis (Gram-Schmidt)
        w_dot_y = array([w_hat]).dot(y_axis_nh2)[0]  # dot product
        x_axis_nh2 = w_hat - w_dot_y * y_axis_nh2
        x_axis_norm = linalg.norm(x_axis_nh2)
        if x_axis_norm > 1e-6:
            x_axis_nh2 = x_axis_nh2 / x_axis_norm
        else:
            # Fallback if w_hat is parallel to C→N axis
            x_axis_nh2 = array([0.0, 0.0, 1.0]) if abs(y_axis_nh2[2]) < 0.9 else array([1.0, 0.0, 0.0])
        
        # Place hydrogens symmetrically at ±60° from C→N axis in the amide plane
        angle_h_half = (AMIDE_PARAMS['geometry']['amide_plane_angles']['n_h_h'] / 2.0) * (pi / 180.0)
        
        h1_vec = cos(angle_h_half) * y_axis_nh2 + sin(angle_h_half) * x_axis_nh2
        h2_vec = cos(angle_h_half) * y_axis_nh2 - sin(angle_h_half) * x_axis_nh2
        
        h1_pos = n_pos + bond_n_h * h1_vec
        h2_pos = n_pos + bond_n_h * h2_vec
        
        # Validate geometry: ensure H atoms are not too close to amide carbon
        min_h_c_distance = 1.3  # Minimum allowed H-C distance (Å)
        h1_c_dist = linalg.norm(h1_pos - cam_pos)
        h2_c_dist = linalg.norm(h2_pos - cam_pos)
        
        if h1_c_dist < min_h_c_distance or h2_c_dist < min_h_c_distance:
            # Flip x-axis and try again
            x_axis_nh2 = -x_axis_nh2
            h1_vec = cos(angle_h_half) * y_axis_nh2 + sin(angle_h_half) * x_axis_nh2
            h2_vec = cos(angle_h_half) * y_axis_nh2 - sin(angle_h_half) * x_axis_nh2
            h1_pos = n_pos + bond_n_h * h1_vec
            h2_pos = n_pos + bond_n_h * h2_vec
        
        # Add atoms in sequence: C(amide), O, N, H1, H2
        atoms_to_add = [
            (['C'] + list(cam_pos) + ['C.2', charge_C_amide]),
            (['O'] + list(o_pos) + ['O.2', charge_O_carbonyl]),
            (['N'] + list(n_pos) + ['N.am', charge_N_amide]),
            (['H'] + list(h1_pos) + ['H', charge_H_amide]),
            (['H'] + list(h2_pos) + ['H', charge_H_amide])
        ]
        
        for idx, atom_data in enumerate(atoms_to_add):
            coords.append(atom_data)
            # Track amide hydrogen indices (last two atoms in each amide group)
            if idx >= 3:  # H1 and H2 are indices 3 and 4
                amide_h_indices.add(len(coords) - 1)
        
        # Track maximum extension for top end
        for atom_data in atoms_to_add:
            extension = max(0, atom_data[2] - c0_y)
            max_extension_top = max(max_extension_top, extension)
    
    # Count atoms added at top end: amide groups + hydrogens
    functionalized_count_top = sum(1 for i in range(natx) if modf(float(i)/2.0)[0] == 0)
    hydrogen_count_top = sum(1 for i in range(natx) if modf(float(i)/2.0)[0] != 0)
    added_atoms = 5 * functionalized_count_top + hydrogen_count_top
    
    # Add amide groups to every second carbon at bottom end
    bottom_start_idx = len(coords) - added_atoms - natx
    for idx, i in enumerate(range(bottom_start_idx, bottom_start_idx + natx)):
        if modf(float(idx)/2.0)[0] != 0:  # Unfunctionalized carbons get hydrogen atoms
            # Update carbon charge for C-H bonding
            coords[i][5] = charge_C_passivate  # Update the charge of the carbon atom
            
            # Add hydrogen atom
            c0_x, c0_y, c0_z = coords[i][1], coords[i][2], coords[i][3]
            
            # Calculate radial direction for H placement
            r_vec = array([c0_x, 0.0, c0_z])
            r_mag = linalg.norm(r_vec)
            if r_mag > 0:
                r_hat = r_vec / r_mag
            else:
                r_hat = array([1.0, 0.0, 0.0])
            
            # Place hydrogen radially outward with slight downward tilt
            y_hat = array([0.0, -1.0, 0.0])  # Negative Y for bottom end
            h_direction = 0.9 * r_hat + 0.1 * y_hat  # Mostly radial with slight downward tilt
            h_direction = h_direction / linalg.norm(h_direction)
            
            h_pos = array([c0_x, c0_y, c0_z]) + bond_c_h_passivate * h_direction
            
            tmpcoords = ['H'] + list(h_pos) + ['H', charge_H_passivate]
            coords.append(tmpcoords)
            continue  # Skip amide functionalization for this carbon
        # Get terminal carbon coordinates
        c0_x, c0_y, c0_z = coords[i][1], coords[i][2], coords[i][3]
        
        # Calculate local coordinate system
        # Radial unit vector (project to xz-plane and normalize)
        r_vec = array([c0_x, 0.0, c0_z])
        r_mag = linalg.norm(r_vec)
        if r_mag > 0:
            r_hat = r_vec / r_mag
        else:
            r_hat = array([1.0, 0.0, 0.0])
        
        # Axial unit vector (negative Y for bottom end)
        y_hat = array([0.0, -1.0, 0.0])
        
        # Tangent unit vector
        t_hat = cross(y_hat, r_hat)
        t_mag = linalg.norm(t_hat)
        if t_mag > 0:
            t_hat = t_hat / t_mag
        else:
            t_hat = array([0.0, 0.0, 1.0])
        
        # Apply staggering if enabled
        if stagger and idx % 2 == 1:
            t_hat = -t_hat
        
        # Bond axis direction (tilted from radial)
        v_bond = cos(theta) * r_hat + sin(theta) * y_hat
        v_bond = v_bond / linalg.norm(v_bond)
        
        # Place amide carbon
        cam_pos = array([c0_x, c0_y, c0_z]) + bond_c_cnt_c_amide * v_bond
        
        # Define in-plane basis at amide carbon
        u_hat = -v_bond  # toward CNT
        w_hat = t_hat    # tangent direction
        
        # Place oxygen and nitrogen (120° separation in amide plane)
        angle_o = 120 * pi / 180
        angle_n = -120 * pi / 180
        
        dir_o = cos(angle_o) * u_hat + sin(angle_o) * w_hat
        dir_n = cos(angle_n) * u_hat + sin(angle_n) * w_hat
        
        o_pos = cam_pos + bond_c_o_double * dir_o
        n_pos = cam_pos + bond_c_n_amide * dir_n
        
        # Place hydrogens on nitrogen using robust orthonormal coordinate system
        # Build N-centric orthonormal frame to ensure proper H-N-H geometry
        y_axis_nh2 = (n_pos - cam_pos) / linalg.norm(n_pos - cam_pos)  # C→N direction
        
        # Orthonormalize tangent vector against C→N axis (Gram-Schmidt)
        w_dot_y = array([w_hat]).dot(y_axis_nh2)[0]  # dot product
        x_axis_nh2 = w_hat - w_dot_y * y_axis_nh2
        x_axis_norm = linalg.norm(x_axis_nh2)
        if x_axis_norm > 1e-6:
            x_axis_nh2 = x_axis_nh2 / x_axis_norm
        else:
            # Fallback if w_hat is parallel to C→N axis
            x_axis_nh2 = array([0.0, 0.0, 1.0]) if abs(y_axis_nh2[2]) < 0.9 else array([1.0, 0.0, 0.0])
        
        # Place hydrogens symmetrically at ±60° from C→N axis in the amide plane
        angle_h_half = (AMIDE_PARAMS['geometry']['amide_plane_angles']['n_h_h'] / 2.0) * (pi / 180.0)
        
        h1_vec = cos(angle_h_half) * y_axis_nh2 + sin(angle_h_half) * x_axis_nh2
        h2_vec = cos(angle_h_half) * y_axis_nh2 - sin(angle_h_half) * x_axis_nh2
        
        h1_pos = n_pos + bond_n_h * h1_vec
        h2_pos = n_pos + bond_n_h * h2_vec
        
        # Validate geometry: ensure H atoms are not too close to amide carbon
        min_h_c_distance = 1.3  # Minimum allowed H-C distance (Å)
        h1_c_dist = linalg.norm(h1_pos - cam_pos)
        h2_c_dist = linalg.norm(h2_pos - cam_pos)
        
        if h1_c_dist < min_h_c_distance or h2_c_dist < min_h_c_distance:
            # Flip x-axis and try again
            x_axis_nh2 = -x_axis_nh2
            h1_vec = cos(angle_h_half) * y_axis_nh2 + sin(angle_h_half) * x_axis_nh2
            h2_vec = cos(angle_h_half) * y_axis_nh2 - sin(angle_h_half) * x_axis_nh2
            h1_pos = n_pos + bond_n_h * h1_vec
            h2_pos = n_pos + bond_n_h * h2_vec
        
        # Add atoms in sequence: C(amide), O, N, H1, H2
        atoms_to_add = [
            (['C'] + list(cam_pos) + ['C.2', charge_C_amide]),
            (['O'] + list(o_pos) + ['O.2', charge_O_carbonyl]),
            (['N'] + list(n_pos) + ['N.am', charge_N_amide]),
            (['H'] + list(h1_pos) + ['H', charge_H_amide]),
            (['H'] + list(h2_pos) + ['H', charge_H_amide])
        ]
        
        for idx, atom_data in enumerate(atoms_to_add):
            coords.append(atom_data)
            # Track amide hydrogen indices (last two atoms in each amide group)
            if idx >= 3:  # H1 and H2 are indices 3 and 4
                amide_h_indices.add(len(coords) - 1)
        
        # Track maximum extension for bottom end
        for atom_data in atoms_to_add:
            extension = max(0, c0_y - atom_data[2])
            max_extension_bottom = max(max_extension_bottom, extension)
    
    # Return the maximum extension per end and amide hydrogen indices
    amide_extension = max(max_extension_top, max_extension_bottom)
    return amide_extension, amide_h_indices

def create_bonds(coords, natx, nohcoords, amide_h_indices=None):
    '''Build connectivity using distance-based method
    
    Args:
        coords: List of atom coordinates and properties
        natx: Number of atoms in CNT cross-section
        nohcoords: Number of CNT backbone atoms (without functional groups)
        amide_h_indices: Set of indices for hydrogen atoms in amide groups
    '''
    if amide_h_indices is None:
        amide_h_indices = set()
    # Extract parameters from organized dictionaries
    radii = COVALENT_RADII
    tolerances = BOND_TOLERANCES
    bond_lengths = AMIDE_PARAMS['bond_lengths']
    
    # Calculate bond tolerances
    btoll_cc_aromatic = (2 * radii['C']) * tolerances['cc_aromatic']
    btoll_ch = (radii['H'] + radii['C']) * tolerances['ch_single']
    btoll_oh = (radii['H'] + radii['O']) * tolerances['oh_single']
    btoll_co_single = (radii['C'] + radii['O']) * tolerances['co_single']
    btoll_cn_amide = bond_lengths['c_n_amide'] * tolerances['cn_amide']
    btoll_nh_amide = (radii['N'] + radii['H']) * tolerances['nh_amide']
    btoll_co_double = bond_lengths['c_o_double'] * tolerances['co_double']
    btoll_cc_amide = bond_lengths['c_cnt_c_amide'] * tolerances['cc_amide']
    
    # Bond distance ranges
    range_cc_aromatic = [2*radii['C'] - btoll_cc_aromatic, 2*radii['C'] + btoll_cc_aromatic]
    range_ch = [(radii['H'] + radii['C']) - btoll_ch, (radii['H'] + radii['C']) + btoll_ch]
    range_co_single = [(radii['O'] + radii['C']) - btoll_co_single, (radii['O'] + radii['C']) + btoll_co_single]
    range_oh = [(radii['H'] + radii['O']) - btoll_oh, (radii['H'] + radii['O']) + btoll_oh]
    range_cc_amide = [bond_lengths['c_cnt_c_amide'] - btoll_cc_amide, bond_lengths['c_cnt_c_amide'] + btoll_cc_amide]
    range_co_double = [bond_lengths['c_o_double'] - btoll_co_double, bond_lengths['c_o_double'] + btoll_co_double]
    range_cn_amide = [bond_lengths['c_n_amide'] - btoll_cn_amide, bond_lengths['c_n_amide'] + btoll_cn_amide]
    range_nh_amide = [(radii['N'] + radii['H']) - btoll_nh_amide, (radii['N'] + radii['H']) + btoll_nh_amide]
    
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
                if bond >= range_cc_aromatic[0] and bond < range_cc_aromatic[1]:
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
                bond_type = '1'  # default single bond
                
                if ((atom1_type == 'C' and atom2_type == 'H') or 
                    (atom1_type == 'H' and atom2_type == 'C')) and \
                   (bond >= range_ch[0] and bond < range_ch[1]):
                    # Prevent C-H bonds with amide hydrogens - they should only bond to N
                    h_idx = j if atom2_type == 'H' else i
                    if h_idx not in amide_h_indices:
                        valid_bond = True
                
                elif ((atom1_type == 'C' and atom2_type == 'O') or 
                      (atom1_type == 'O' and atom2_type == 'C')):
                    # Check for C=O double bond first
                    if bond >= range_co_double[0] and bond < range_co_double[1]:
                        valid_bond = True
                        bond_type = '2'  # double bond
                    # Then check for C-O single bond
                    elif bond >= range_co_single[0] and bond < range_co_single[1]:
                        valid_bond = True
                        bond_type = '1'  # single bond
                
                elif ((atom1_type == 'C' and atom2_type == 'N') or 
                      (atom1_type == 'N' and atom2_type == 'C')) and \
                     (bond >= range_cn_amide[0] and bond < range_cn_amide[1]):
                    valid_bond = True
                
                elif ((atom1_type == 'N' and atom2_type == 'H') or 
                      (atom1_type == 'H' and atom2_type == 'N')) and \
                     (bond >= range_nh_amide[0] and bond < range_nh_amide[1]):
                    valid_bond = True
                
                elif ((atom1_type == 'O' and atom2_type == 'H') or 
                      (atom1_type == 'H' and atom2_type == 'O')) and \
                     (bond >= range_oh[0] and bond < range_oh[1]):
                    valid_bond = True
                
                elif (atom1_type == 'C' and atom2_type == 'C'):
                    # Check for C(CNT)-C(amide) single bond first
                    if bond >= range_cc_amide[0] and bond < range_cc_amide[1]:
                        valid_bond = True
                        bond_type = '1'  # single bond
                    # Then check for CNT C-C aromatic bonds
                    elif bond >= range_cc_aromatic[0] and bond < range_cc_aromatic[1]:
                        valid_bond = True
                        bond_type = 'ar'  # aromatic bond (for CNT carbons)
                
                if valid_bond:
                    if i < j:
                        bondnumber += 1
                        # Use aromatic bond type for CNT C-C bonds, otherwise use determined bond_type
                        if (atom1_type == 'C' and atom2_type == 'C' and 
                            bond >= range_cc_aromatic[0] and bond < range_cc_aromatic[1] and
                            i < nohcoords and j < nohcoords):
                            tmpbond = [bondnumber, i+1, j+1, 'ar']
                        else:
                            tmpbond = [bondnumber, i+1, j+1, bond_type]
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
    '''Write mol2 file for amide-functionalized CNT'''
    file.write("@<TRIPOS>MOLECULE\n")
    file.write("AMIDE_CNT\n")
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
    """Main function to build amide-functionalized CNT structures."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Extract CNT structural parameters
    ccbond = CNT_PARAMS['cc_bond_length']
    
    # Create CNT structure with specified backbone length
    coords, natx, pbc_l, nohcoords, actual_backbone_length = create_cnt_backbone(
        args.index, args.length, ccbond)
    
    # Add amide groups to every second carbon at both ends
    amide_extension, amide_h_indices = add_amide_groups(coords, natx, args.amide_angle, args.stagger)
    
    # Calculate total length including amide groups
    total_length = actual_backbone_length + (2 * amide_extension)
    
    # Calculate radius and diameter
    dx = ccbond * cos(120/2 * (pi/180.0))
    radius = (args.index * (2*dx + ccbond) + args.index * ccbond) / (2*pi)
    diameter = 2 * radius
    
    # Display results
    print('\n' + '='*60)
    print('AMIDE-FUNCTIONALIZED ARMCHAIR CNT STRUCTURE')
    print('='*60)
    print(f'  Nanotube index (n)           = {args.index}')
    print(f'  CNT backbone length (Å)      = {actual_backbone_length:.3f}')
    print(f'  Amide extension per end (Å)  = {amide_extension:.3f}')
    print(f'  Total length with amide (Å)  = {total_length:.3f}')
    print(f'  Diameter (Å)                 = {diameter:.3f}')
    print(f'  Periodicity (Å)              = {pbc_l:.3f}')
    print(f'  Amide tilt angle (degrees)   = {args.amide_angle:.1f}')
    print(f'  Staggered amide planes       = {args.stagger}')
    print('='*60)
    print(f'  Total atoms                  = {len(coords)}')
    
    # Create bond connectivity
    conn, bondlist = create_bonds(coords, natx, nohcoords, amide_h_indices)
    print(f'  Total bonds                  = {len(bondlist)}')
    
    # Validate amide geometry and bonding
    geometry_validation = validate_amide_geometry(coords, amide_h_indices)
    bonding_validation = validate_amide_bonding(coords, bondlist, amide_h_indices)
    
    print(f'  Amide hydrogens              = {geometry_validation["total_amide_h"]}')
    print(f'  Min H-C distance (Å)         = {geometry_validation["min_h_c_distance"]:.3f}')
    print(f'  C-N amide bonds              = {bonding_validation["cn_bonds"]}')
    print(f'  C=O double bonds             = {bonding_validation["co_double_bonds"]}')
    print(f'  N-H bonds                    = {bonding_validation["nh_bonds"]}')
    
    if geometry_validation['problematic_h_indices']:
        print(f'  WARNING: {len(geometry_validation["problematic_h_indices"])} H atoms too close to C!')
    if bonding_validation['spurious_ch_bonds'] > 0:
        print(f'  WARNING: {bonding_validation["spurious_ch_bonds"]} spurious C-H bonds with amide H!')
    
    # Calculate expected values for validation
    functionalized_carbons = sum(1 for i in range(args.index) if i % 2 == 0) * 2  # Both ends
    expected_amide_groups = functionalized_carbons
    print(f'  Expected amide groups        = {expected_amide_groups}')
    
    print('='*60)
    
    # Write mol2 file
    with open(args.output_file, 'w') as outfile:
        write_mol2(outfile, coords, bondlist)
    
    print(f'\nAmide-functionalized CNT saved to: {args.output_file}')
if __name__ == "__main__":
    main()
