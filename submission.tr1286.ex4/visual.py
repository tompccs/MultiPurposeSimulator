#!/usr/bin/env python 

from graphics import *
from sys import argv, stdin
from random import randint

fileout = None

if len(argv) > 1:
	scale = float(argv[1])
	if len(argv) > 2:
		fileout_path = argv[2]
		fileout = open(fileout_path, "w")

noflush = False
if "--noflush" in argv:
	noflush = True

width = 600
height = 600

flush_wait = 5
flush_counter = 0

colour = 0x000000
timestep = 0

timetemp = 0

win = GraphWin("Simulation", width, height, autoflush=False)
c2b = Text(Point(width/2, height/2), "Click to begin")
c2b.draw(win)
win.getMouse()
c2b.undraw()
if noflush:
	c2b.setText("Drawing...")
	c2b.draw(win)
win.flush()

key_line = stdin.readline()
if fileout:
	fileout.write(key_line)

# determine the number of bodies
bodies = []
xpos_indecies = []
index = 0
keys = key_line.split(",")
for key in keys:
	split_key = key.split(".")
	if key == "xpos" or (len(split_key) > 1 and split_key[1] == "xpos"):
		bodies.append(Circle(Point(0, 0), 5))
		xpos_indecies.append(index)
	index += 1

index = 0
for body in bodies:
	colour = color_rgb(randint(0, 255), randint(0, 255), randint(0, 255))
	body.setFill(colour)
	keytext = Text(Point(width - 100, 10 * (index + 1)), "body "+str(index))
	keytext.setTextColor(colour)
	body.draw(win)
	keytext.draw(win)
	index += 1

timebox = Rectangle(Point(3, 3), Point(100, 20))
timebox.setFill('yellow')
timebox.draw(win)
timetext = Text(Point(50, 10), 'time = ')
timetext.draw(win)

lastpos = {}

for line in stdin:
	if fileout:
		fileout.write(line)
	sline = line.split(",")
	time = float(sline[0])

	if timetemp == 0:
		timetemp = time
	else:
		timestep = time - timetemp
		colour += int(timetemp)
		if colour >= 0xFFFFFF:
			colour = 0x000000

	timetext.setText('time = '+str(time))
	
	xposs = []
	yposs = []
	body_index = 0
	for i in xpos_indecies:
		xpos = float(sline[i])
		ypos = float(sline[i + 1])
		xpos_scaled = xpos / scale * width + width/2
		ypos_scaled = ypos / scale * height + height/2
		if not noflush:
			print keys[i], "=", xpos_scaled
			print keys[i + 1], "=", ypos_scaled
		if (xpos_scaled, ypos_scaled) not in lastpos:
			try:
				win.plotPixel(xpos_scaled, ypos_scaled, color_rgb(colour & 0xFF0000, colour & 0x00FF00, colour & 0x0000FF))
			except Exception:
				colour = 0x000000
			if not noflush:
				bodies[body_index].move(xpos_scaled - bodies[body_index].getCenter().getX(), ypos_scaled - bodies[body_index].getCenter().getY())
			body_index += 1
			lastpos[i] = (xpos_scaled, ypos_scaled)

	if not noflush and flush_counter >= flush_wait:
		win.flush()
		flush_counter = 0
	else:
		flush_counter += 1

if noflush:
	c2b.undraw()
win.flush()
if fileout:
	fileout.close()
win.getMouse()
