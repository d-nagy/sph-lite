import matplotlib.pyplot as plt
import numpy as np

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

inputFilename = '../sph.dat'
outputFilename = 'sedovBlast2D.case'
params = []
inputLines = []

with open(inputFilename) as f:
    inputLines = f.readlines()
    params = [line.strip().split() for line in inputLines]
    params = [line[0] for line in params if line]

fluidParticleSize = float(params[7])
restDensity = float(params[1])
energySpike = 1.0
spikeRadius = 0.01
bgdU = 1e-9

start = (fluidParticleSize/2) - 0.5
end = 0.5
ns = np.arange(start, end, fluidParticleSize)
N = len(ns)**2
pmass = restDensity / N

fluid = []
spikeParticles = []
for i, x in enumerate(ns):
    for j, y in enumerate(ns):
        u = bgdU
        if np.sqrt(x**2 + y**2) < spikeRadius:
            spikeParticles.append(i*len(ns) + j)
        fluid.append([round(x, 5), round(y, 5), 0, 0, u])

uSpike = energySpike / (pmass * len(spikeParticles))

for sp in spikeParticles:
    fluid[sp][4] = uSpike

with open(outputFilename, 'w') as f:
    f.write(f'0 {len(fluid)}\n')
    for p in fluid:
        f.write(' '.join(list(map(str, p))) + '\n')

outputLines = inputLines[:-1] + [outputFilename + '\n']
with open(inputFilename, 'w') as f:
    f.writelines(outputLines)

xs = [p[0] for p in fluid]
ys = [p[1] for p in fluid]
us = [p[4] for p in fluid]

plt.scatter(xs, ys, c=us, s=5)
plt.xlim(-0.6, 0.6)
plt.ylim(-0.6, 0.6)
plt.grid(True)
plt.show()
