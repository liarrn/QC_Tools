def xyz2xsf(xyz_filename, xsf_prefix):
    with open(xyz_filename) as fp:
        xyz = fp.readlines()
    natoms = int(xyz[0])
    num_step = len(xyz)/(natoms+2)
    if len(xyz[2].split()) == 4:
        is_force = False
    elif len(xyz[2].split()) == 7:
        is_force = True
    for i in range(num_step):
        energy = float(xyz[(natoms+2)*i+1].split()[1])
        header = '# total energy = %.9f eV\n\n'%energy
        if is_force == True:
            coor = xyz[(natoms+2)*i+2: (natoms+2)*(i+1)]
        elif is_force == False:
            coor = map(lambda s: s.rstrip()+'%15.9f%15.9f%15.9f\n'%(0.0, 0.0, 0.0), xyz[(natoms+2)*i+2: (natoms+2)*(i+1)])
        coor = reduce(lambda s1, s2: s1+s2, coor)
        with open(xsf_prefix+'%05d'%i+'.xsf', 'w') as fp:
            fp.write(header)
            fp.write(coor)
			
if __name__ == '__main__':
    xyz2xsf('xyz-force.xyz', 'au10_')