import matplotlib.pyplot as plt
import numpy as np

inputFilename = '../sph.dat'
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
x = boundaryParticleSize / 2
while x < 1:
    y = -boundaryParticleSize / 2
    for w in range(boundaryWidth):
        bottom.append((round(x, 5), round(y, 5), 1))
        y -= boundaryParticleSize
    x += boundaryParticleSize

# do the left side
leftside = []
y = bottom[boundaryWidth-1][1]
while y < 0.5:
    x = -boundaryParticleSize / 2 
    for w in range(boundaryWidth):
        leftside.append((round(x, 5), round(y, 5), 1))
        x -= boundaryParticleSize
    y += boundaryParticleSize

# and the right
rightside = []
y = bottom[boundaryWidth-1][1]
while y < 0.5:
    x = bottom[-1][0] + boundaryParticleSize
    for w in range(boundaryWidth):
        leftside.append((round(x, 5), round(y, 5), 1))
        x += boundaryParticleSize
    y += boundaryParticleSize

## Now fill in the fluid
fluid = []
y = fluidParticleSize / 2
while y < 0.2:
    x = fluidParticleSize / 2
    while x < 1:
        hydrostaticPressure = restDensity * 9.81 * (0.2 - y) 
        fluid.append((round(x, 5), round(y, 5), 0, hydrostaticPressure))
        x += fluidParticleSize
    y += fluidParticleSize

boundaryParticles = leftside[::-1] + bottom + rightside

with open(outputFilename, 'w') as f:
    f.write(f'{len(boundaryParticles)} {len(fluid)}\n')
    for p in boundaryParticles + fluid:
        f.write(' '.join(list(map(str, p))) + '\n')

outputLines = inputLines[:-1] + [outputFilename + '\n']
with open(inputFilename, 'w') as f:
    f.writelines(outputLines)

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

