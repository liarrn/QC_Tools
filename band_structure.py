import xml.etree.ElementTree as ET

#for spin polarized calculation !!!
#output format:
#			band1 band2 band3 ...
#kpoint1
#kpoint2
#kpoint3
#...
#...

dir1 = 'band_alpha.dat'
dir2 = 'band_beta.dat'

tree = ET.parse('vasprun.xml')
root = tree.getroot()

bands1 = []
bands2 = []
efermi = float(root.find('calculation/dos/i').text.strip())
spin1 = root.find('calculation/eigenvalues/array/set').findall('set')[0]
spin2 = root.find('calculation/eigenvalues/array/set').findall('set')[1]
for kpoint in spin1:
	k_tmp = []
	for eigen in kpoint:
	[energy, fill] = eigen.text.split()
	[energy, fill] = [float(energy)-efermi, int(float(fill))] #band structure with fermi energy correction
	k_tmp.append({energy: fill})
	bands1.append(k_tmp)
with open(dir1, mode = 'w') as fp:
	for kpoint in bands1:
		for band in kpoint:
			fp.write(str(band.keys()[0]) + ' ')
		fp.write('\n')
for kpoint in spin2:
	k_tmp = []
	for eigen in kpoint:
	[energy, fill] = eigen.text.split()
	[energy, fill] = [float(energy)-efermi, int(float(fill))] #band structure with fermi energy correction
	k_tmp.append({energy: fill})
	bands2.append(k_tmp)
with open(dir2, mode = 'w') as fp:
	for kpoint in bands2:
		for band in kpoint:
			fp.write(str(band.keys()[0]) + ' ')
		fp.write('\n')