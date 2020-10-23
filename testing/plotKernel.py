import matplotlib.pyplot as plt

filename = 'cubicKernelPlot.txt'
lines = []

with open(filename) as f:
    lines = [line.strip() for line in f.readlines()]

h = float(lines[0].split('=')[-1])

qs, Ws, dWs = [], [], []

for line in lines[2:]:
    q, W, dW = map(float, line.split(','))
    qs.append(q)
    Ws.append(W)
    dWs.append(dW)

plt.plot(qs, Ws, label="W")
plt.plot(qs, dWs, label="dW")
plt.grid(True)
plt.xlabel('q')
plt.legend()
plt.title(f'Cubic spline kernel, h={h}')
plt.show()

