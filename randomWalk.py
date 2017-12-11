# Kip Thorne <Modern Classical Physics>
# ex-6.3 random walk example
# Shakuna Vimana
# 2017-12-8
import random
import numpy as np
import matplotlib.pyplot as plt
y0 = 0
yj = []
for i in range(10000):
 y = random.choice([-1,1])
 y0 = y0+y
 yj.append(y0)

hist, bin_edges = np.histogram(yj, density=True)
plt.hist(yj,bins="auto")
plt.show()
