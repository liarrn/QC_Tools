
        
if __name__ == '__main__':
    xyzfilename = './au20-xyz-force-123.xyz'
    ns = 0
    with open(xyzfilename, 'r') as fp:
        while True:
            natoms = fp.readline()
            if not natoms:
                break
            ns += 1
            natoms = int(natoms)
            energy = float(fp.readline().split()[1])
            coor = ''
            for i in range(natoms):
                coor += fp.readline()
            xsf = '# total energy = %.9f eV\n\nATOMS\n'%energy
            xsf += coor
            with open('structure%06d.xsf'%ns, 'w') as fp2:
                fp2.write(xsf)
             
                
    
    
    