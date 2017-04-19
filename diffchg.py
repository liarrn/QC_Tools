#differential charge density -- CHGCAR1+CHGCAR2-CHGCAR3
#positive--lose electron after adsorption
#negative--gain electron after adsorption

CHGCAR1 = 'CHGCAR_fa'
CHGCAR2 = 'CHGCAR_substrate'
CHGCAR3 = 'CHGCAR_all'
CHGCAR_diff = "CHGCAR_diff"

i = j = k = 0
na = 0
ng = 1
chg1 = []
chg2 = []
chg3 = []
with open(CHGCAR1, mode = 'r') as fchg1:
	with open(CHGCAR2, mode = 'r') as fchg2:
		with open(CHGCAR3, mode = 'r') as fchg3:
			with open(CHGCAR_diff, mode = 'w') as fpo:
				for i in range(6):
					fchg1.readline()
					fchg2.readline()
					fpo.write(fchg3.readline())
				for i in fchg1.readline().split():
					na = na + int(i)
				for i in range(na+2):
					fchg1.readline()
				na = 0
				for i in fchg2.readline().split():
					na = na + int(i)
				for i in range(na+2):
					fchg2.readline()
				na = 0
				tmp = fchg3.readline()
				fpo.write(tmp)
				for i in tmp.split():
					na = na + int(i)
				for i in range(na+2):
					fpo.write(fchg3.readline())
				for i in fchg1.readline().split():
					ng = ng * int(i)
				fchg2.readline()
				fpo.write(fchg3.readline())
				for i in range ( (ng + 4) / 5):
					print i
					chg1= fchg1.readline().split()
					chg2= fchg2.readline().split()
					chg3= fchg3.readline().split()
					for j in range(5):
						tmp = float(chg1[j])+float(chg2[j])-float(chg3[j])
						s = '%0.11e'%tmp
						fpo.write(s)
						fpo.write(" ")
					fpo.write("\n")
						