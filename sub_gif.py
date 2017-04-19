from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import StringIO

frame = 961 #time of simulation
tpf = 15.0 #simualation time per frame, fs

for i in range(0, frame):
	a= StringIO.StringIO()
	a.write("%05d"%i)
	img = Image.open('untitled.'  +a.getvalue()+'.bmp')
	draw = ImageDraw.Draw(img)
	font = ImageFont.truetype("calibri.ttf", 40)
	time = i * tpf /1000
	draw.text((300, 590), "%0.1f ps"%time, (0, 0, 0), font = font)
	img.save('sub.'+a.getvalue()+'.bmp')
	a.close()