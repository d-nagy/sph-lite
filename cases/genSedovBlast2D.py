import matplotlib.pyplot as plt
import numpy as np

def gaussianKernel(x, y, sig):
    return np.exp(-(x**2 + y**2) / (2 * sig**2)) / (2 * np.pi * sig**2)

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
h = float(params[9])
kernel = params[10]
gamma = float(params[4])
energySpike = 0.1
spikeRadius = 0.08 # fluidParticleSize * 1.2
bgdU = 1e-5 / (restDensity * (gamma - 1))

if kernel == 'cubicspline':
    supportRadius = 2 * h

start = (fluidParticleSize/2) - 0.5
end = int(np.ceil(1 / supportRadius)) * (supportRadius / 2)
ns_ = np.arange(0, end, fluidParticleSize)
ns = np.concatenate((-ns_[::-1], ns_[1:]))
N = len(ns)**2
pmass = restDensity / N

fluid = []
spikeParticles = []
for i, x in enumerate(ns):
    for j, y in enumerate(ns):
        u = bgdU
        dist = np.sqrt(x**2 + y**2)
        if dist < spikeRadius:
            spikeParticles.append(i*len(ns) + j)
        fluid.append([round(x, 5), round(y, 5), 0, 0, u])

# uSpike = energySpike / (pmass * len(spikeParticles))

for sp in spikeParticles:
    fluid[sp][4] = round(energySpike * gaussianKernel(fluid[sp][0], fluid[sp][1], spikeRadius/3), 5)

with open(outputFilename, 'w') as f:
    f.write(f'0 {len(fluid)}\n')
    domainMin = ns[0] - fluidParticleSize/2
    domainMax = ns[-1] + fluidParticleSize/2
    f.write(f'{domainMin} {domainMin} 0 {domainMax} {domainMax} 0\n')
    for p in fluid:
        f.write(' '.join(list(map(str, p))) + '\n')

outputLines = inputLines[:-1] + [outputFilename + '\n']
with open(inputFilename, 'w') as f:
    f.writelines(outputLines)

xs = [p[0] for p in fluid]
ys = [p[1] for p in fluid]
us = [p[4] for p in fluid]

usArr = np.array(us)
print(sum(pmass * usArr[usArr > bgdU]))

plt.scatter(xs, ys, c=us, s=5)
plt.xlim(-0.6, 0.6)
plt.ylim(-0.6, 0.6)
plt.grid(True)
plt.show()
