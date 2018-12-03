# http://matplotlib.org/


import numpy as np
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import scipy.interpolate

plt.close("all")

filename = "output.out"
test = False

timeUNIT = 1
distUNIT = 0.0529
voltageUNIT = 27.21138
chargeUNIT = 1
currentUNIT = 1

lc = 0
lt = 0
variables = []
x = []
y = []
types = []
numofcolor = 0
colortypelist = []
numoffields = 0
u = []
q = []
varlist = []

maxx = 0
maxy = 0
maxu = 0
maxq = 0

with open(filename, "r") as f:
    for line in f:
        # print lt

        if (lt == 0):
            variables = [float(item) for item in line.split()[:-1]]
            varlist += [variables]

        if (lt == 1): maxx = max(maxx, max([float(item) for item in line.split()]))
        if (lt == 2): maxy = max(maxy, max([float(item) for item in line.split()]))
        if (lt == 4): numofcolor = int(line)
        if (lt == numofcolor + 5): numoffields = int(line)
        if (lt == numofcolor + 5 + 1): maxu = max(maxu, max([abs(float(item)) for item in line.split()]))
        if (lt == numofcolor + 5 + 2): maxq = max(maxq, max([abs(float(item)) for item in line.split()]))

        if (lt == numoffields + numofcolor + 6): lt = -1
        lt += 1

maxx *= distUNIT
maxy *= distUNIT
maxu *= voltageUNIT
maxq *= chargeUNIT
print maxx, maxy, maxu, maxq

fig = matplotlib.pyplot.figure(figsize=(6, 5))
timeSERIES = [(v[0]) * timeUNIT for v in varlist]
biasSERIES = [(v[2] - v[1]) * voltageUNIT for v in varlist]
currentSERIES = [(v[3]) * currentUNIT for v in varlist]
tfig1 = fig.add_subplot(2, 1, 1)
tfig1.set_xlim(min(timeSERIES), max(timeSERIES))
tfig1.set_xlabel("t/t0")
tfig1.set_ylim(min(biasSERIES), max(biasSERIES))
tfig1.set_ylabel("Ub/V")
tfig1.plot(timeSERIES, biasSERIES)
tfig1 = fig.add_subplot(2, 1, 2)
tfig1.set_xlim(min(timeSERIES), max(timeSERIES))
tfig1.set_xlabel("t/t0")
tfig1.set_ylim(min(currentSERIES), max(currentSERIES))
tfig1.set_ylabel("I/I0")
tfig1.plot(timeSERIES, currentSERIES)
fig.tight_layout()
fig.savefig("pngs/result_timevs.png", dpi=60)
plt.close()

lc = 0
lt = 0
with open(filename, "r") as f:
    for line in f:
        # print lt

        if (lt == 0): variables = [float(item) for item in line.split()[:-1]]
        if (lt == 1): x = [float(item) for item in line.split()]
        if (lt == 2): y = [float(item) for item in line.split()]
        if (lt == 3): types = [int(item) for item in line.split()]
        if (lt == 4): numofcolor = int(line)
        if (5 <= lt < numofcolor + 5):
            items = line.split()
            colortypelist += [(items[0], [int(item) for item in items[1:]])]
        if (lt == numofcolor + 5): numoffields = int(line)
        if (lt == numofcolor + 5 + 1): u = [float(item) for item in line.split()]
        if (lt == numofcolor + 5 + 2): q = [float(item) for item in line.split()]
        if (lt == numofcolor + 5 + 3): IFx = [float(item) for item in line.split()]
        if (lt == numofcolor + 5 + 4): IFy = [float(item) for item in line.split()]
        if (lt == numofcolor + 5 + 5): EFx = [float(item) for item in line.split()]
        if (lt == numofcolor + 5 + 6): EFy = [float(item) for item in line.split()]

        if (lt == numoffields + numofcolor + 6):
            # print variables,x,y,types,numofcolor,colortypelist,numoffields
            lt = -1

            x = np.array(x) * distUNIT
            y = np.array(y) * distUNIT
            u = np.array(u) * voltageUNIT
            q = np.array(q) * chargeUNIT

            npx = np.array(x)
            npy = np.array(y)
            # npu = np.array(u)
            # npq = np.array(q)

            npxmin = 0  # npx.min()
            npxmax = maxx  # npx.max()
            npymin = 0  # npy.min()
            npymax = maxy  # npy.max()

            xi, yi = np.linspace(npxmin, npxmax, 300), np.linspace(npymin, npymax, 300)
            xi, yi = np.meshgrid(xi, yi)

            time = variables[0] * timeUNIT
            # bias = (variables[2] - variables[1])*voltageUNIT
            # current = variables[3]*currentUNIT

            fig = matplotlib.pyplot.figure(figsize=(14, 8))
            # fig.suptitle("t/t0 = % 3.3e Ub/V = % 3.3e I/I0 = % 3.3e  " % (time,bias,current))
            fig.suptitle("t/t0 = % g " % (time))
            subfig1 = fig.add_subplot(1, 2, 1)
            subfig1.set_xlim(npxmin, npxmax)
            subfig1.set_xlabel("x/a")
            subfig1.set_ylim(npymin, npymax)
            subfig1.set_ylabel("y/a")
            subfig1.set_aspect('equal')
            subfig1.set_title('surface charge')
            subfig1ax = make_axes_locatable(subfig1)
            subfig1cax = subfig1ax.append_axes("right", size="10%", pad=0.05)
            qz = scipy.interpolate.griddata((npx, npy), q, (xi, yi), method='linear')
            subfig1qzimshow = subfig1.imshow(qz, cmap='seismic', vmin=-maxq, vmax=maxq, origin='left',
                                             extent=[npxmin, npxmax, npymin, npymax])
            # subfig1qzimshow = subfig1.imshow(q, cmap='RdYlBu_r', vmin=-maxq, vmax=maxq, origin='left',  extent=[npxmin, npxmax, npymin, npymax],interpolation="nearest" )
            cbar1 = fig.colorbar(subfig1qzimshow, cax=subfig1cax)
            cbar1.set_label("q/e")
            # subfig1.quiver(x,y,EFx,EFy,pivot='mid',scale_units='xy', angles='xy', scale=0.1,color='g',width=0.004)

            for (c, t) in colortypelist: subfig1.plot(npx[t], npy[t], '.', color=c, markeredgecolor='none',
                                                      markersize=3)

            subfig2 = fig.add_subplot(1, 2, 2)
            subfig2.set_xlim(npxmin, npxmax)
            subfig2.set_xlabel("x/a")
            subfig2.set_ylim(npymin, npymax)
            subfig2.set_aspect('equal')
            subfig2.set_title('potential')
            subfig2ax = make_axes_locatable(subfig2)
            subfig2cax = subfig2ax.append_axes("right", size="10%", pad=0.05)
            uz = scipy.interpolate.griddata((npx, npy), u, (xi, yi), method='linear')
            subfig2uzimshow = subfig2.imshow(uz, cmap='seismic', vmin=-maxu, vmax=maxu, origin='left',
                                             extent=[npxmin, npxmax, npymin,
                                                     npymax])  # http://matplotlib.org/users/colormaps.html
            cbar2 = fig.colorbar(subfig2uzimshow, cax=subfig2cax)
            cbar2.set_label('U/V')
            # subfig2.quiver(x,y,IFx,IFy,pivot='mid',scale_units='xy', angles='xy', scale=0.1,color='m',width=0.004)

            for (c, t) in colortypelist: subfig2.plot(npx[t], npy[t], '.', color=c, markeredgecolor='none',
                                                      markersize=3)

            fig.tight_layout()
            fig.savefig("pngs/result%d.png" % (lc + 100000), dpi=260)
            plt.close()
            print '.',

        lt += 1
        lc += 1
        if test:
            if lc > 30: break
