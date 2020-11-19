from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

inputFilename = '../sedov.dat'
outputFilename = 'sedovBlast3D.case'
params = []
inputLines = []

with open(inputFilename) as f:
    inputLines = f.readlines()
    params = [line.strip().split() for line in inputLines]
    params = [line[0] for line in params if line]

fluidParticleSize = float(params[4])
restDensity = float(params[1])
h = float(params[6])
kernel = params[7]
gamma = float(params[10])
energySpike = 0.1
spikeRadius = 0.08 # fluidParticleSize * 1.2
bgdU = 1e-5 / (restDensity * (gamma - 1))

if kernel == 'cubicspline':
    supportRadius = 2 * h

start = (fluidParticleSize/2) - 0.5
end = int(np.ceil(1 / supportRadius)) * (supportRadius / 2)
ns_ = np.arange(0, end, fluidParticleSize)
ns = np.concatenate((-ns_[::-1], ns_[1:]))
N = len(ns)**3
pmass = restDensity / N

fluid = []
spikeParticles = []
for i, x in enumerate(ns):
    for j, y in enumerate(ns):
        for k, z in enumerate(ns):
            u = bgdU
            if np.sqrt(x**2 + y**2 + z**2) < spikeRadius:
                spikeParticles.append(i*(len(ns)**2) + j*len(ns) + k)
            fluid.append([round(x, 5), round(y, 5), round(z, 5), 0, 0, u])

uSpike = energySpike / (pmass * len(spikeParticles))

for sp in spikeParticles:
    fluid[sp][5] = uSpike

with open(outputFilename, 'w') as f:
    f.write(f'0 {len(fluid)}\n')
    domainMin = ns[0] - fluidParticleSize/2
    domainMax = ns[-1] + fluidParticleSize/2
    f.write(f'{domainMin} {domainMin} {domainMin} {domainMax} {domainMax} {domainMax}\n')
    for p in fluid:
        f.write(' '.join(list(map(str, p))) + '\n')

xs = [p[0] for p in fluid]
ys = [p[1] for p in fluid]
zs = [p[2] for p in fluid]
us = [p[5] for p in fluid]

print(len(fluid))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zs, c=us, s=1)
plt.grid(True)
plt.show()
