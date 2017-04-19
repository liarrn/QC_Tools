#only for movies shorter than 60s

movTime = 17
trjTime = 2780

unitTime = movTime/trjTime*1000*100
for i in range(int(trjTime/100)+1):
	print i+1
	print "00:00:"+"%02d"%(i*unitTime/1000)+",%03d"%(i*unitTime%1000)+ " --> ""00:00:"+"%02d"%((i+1)*unitTime/1000)+",%03d"%((i+1)*unitTime%1000)
	print "Time: %.1f"%(float(i)/10)+" <i>ps</i>"
	print "\n"
	