import matplotlib.pyplot as plt
import numpy as np

inputFilename = '../stationary.dat'
outputFilename = 'stationaryWater2D.case'
params = []
inputLines = []

with open(inputFilename) as f:
    inputLines = f.readlines()
    params = [line.strip().split() for line in inputLines]
    params = [line[0] for line in params if line]

fluidParticleSize, boundaryParticleSize, smoothingLength = float(params[4]), float(params[5]), float(params[6])
restDensity = float(params[1])
boundaryWidth = int(np.ceil(smoothingLength*2 / boundaryParticleSize))

## Generate boundary particles for fluid container
# start with bottom
bottom = []
startx, endx = boundaryParticleSize / 2, 1
xs = np.arange(startx, endx, boundaryParticleSize)
ys = [i * -boundaryParticleSize - (boundaryParticleSize / 2) for i in range(boundaryWidth)]
for x in xs:
    for y in ys:
        bottom.append((round(x, 5), round(y, 5), 1, 0, 0))

# do the sides
leftside, rightside = [], []
starty, endy = bottom[boundaryWidth-1][1], 0.5
ys = np.arange(starty, endy, boundaryParticleSize)
lxs = [i * -boundaryParticleSize - (boundaryParticleSize / 2) for i in range(boundaryWidth)]
rxs = [(i+1) * boundaryParticleSize + (bottom[-1][0]) for i in range(boundaryWidth)]
for lx, rx in zip(lxs, rxs):
    for y in ys:
        leftside.append((round(lx, 5), round(y, 5), 1, 0, 0))
        rightside.append((round(rx, 5), round(y, 5), 1, 0, 0))

## Now fill in the fluid
fluid = []
startx, endx = fluidParticleSize / 2, 1
starty, endy = fluidParticleSize / 2, 0.2
xs = np.arange(startx, endx, fluidParticleSize)
ys = np.arange(starty, endy, fluidParticleSize)
for x in xs:
    for y in ys:
        hydrostaticPressure = restDensity * 9.81 * (0.2 - y)
        fluid.append((round(x, 5), round(y, 5), 0, hydrostaticPressure, 0))

boundaryParticles = leftside[::-1] + bottom + rightside

with open(outputFilename, 'w') as f:
    f.write(f'{len(boundaryParticles)} {len(fluid)}\n')
    f.write(f'-0.2 -0.2 0 1.2 1 0\n')
    for p in boundaryParticles + fluid:
        f.write(' '.join(list(map(str, p))) + '\n')

bxs = [p[0] for p in boundaryParticles]
bys = [p[1] for p in boundaryParticles]
fxs = [p[0] for p in fluid]
fys = [p[1] for p in fluid]
plt.plot(bxs, bys, 'go')
plt.plot(fxs, fys, 'bo')
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 0.6)
plt.grid(True)
plt.show()
