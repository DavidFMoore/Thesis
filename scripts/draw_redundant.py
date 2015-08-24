#! /usr/bin/env python

from pylab import *

fig = figure(figsize=(8,4))

ax = fig.add_subplot(111)

for col in range(8):
    for row in range(4):
        plot(col, row,'ko')

#Three types
plot((0,1), (2,3), 'b', lw=2)
text(0.0, 2.6, "(1,1)", size=15)
plot((0,1), (2,2), 'b', lw=2)
text(0.5, 2.1, "(1,0)", size=15)
plot((0,1), (2,1), 'b', lw=2)
text(0.6, 1.5, "(1,-1)", size=15)

#closed triangle
plot((4,5), (0,0), 'c', lw=2)
plot((5,6), (0,0), 'c', lw=2)
plot((4,5), (0,1), 'c', lw=2)
plot((5,6), (1,0), 'c', lw=2)

axis('off')
xticks([])
yticks([])
xlim([-0.2,7.2])
ylim([-0.2,3.2])
savefig('Tools/figures/AntennaGrid.eps', fmt='eps')
show()
