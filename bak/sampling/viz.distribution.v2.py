

########################################################################

# Read count coverage. till this part everyhting is the same
#3 when snps are actually being called is when the negative binomial enters

# identifies true variants POS, ALT, REF -
# filter indivisual sequence/s by variable POS
# get TRUE genotype per position
# rc-sample
## sample coverage for position
## if 1 seq
### compute error
## else if 2 seq
### sample binomial 50% each strand (coverage)
### compute error
### DP=coverage

## when there is already the RC i can compute GL,AD
QUAL????? Hownto calculate it?
GL quality- is it just ohred score?


# get genotype per position


import numpy as np
from scipy.stats import nbinom, poisson
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)
n, p = 100, 0.4
mean, var, skew, kurt = nbinom.stats(n, p, moments='mvsk')
x = np.arange(nbinom.ppf(0.01, n, p),nbinom.ppf(0.99, n, p))
ax.plot(x, nbinom.pmf(x, n, p), 'bo', ms=8, label='nbinom pmf')
plt.show()

fig, ax = plt.subplots(1, 1) # mean is half n
n, p = 100, 0.5
mean, var, skew, kurt = nbinom.stats(n, p, moments='mvsk')
x = np.arange(nbinom.ppf(0.01, n, p),nbinom.ppf(0.99, n, p))
ax.plot(x, nbinom.pmf(x, n, p), 'bo', ms=8, label='nbinom pmf')
plt.show()


fig, ax = plt.subplots(2, 1)
n, p = 100/2, 0.5
mean, var, skew, kurt = nbinom.stats(n, p, moments='mvsk')
x = np.arange(nbinom.ppf(0.01, n, p),nbinom.ppf(0.99, n, p))
ax[0].plot(x, nbinom.pmf(x, n, p), 'bo', ms=8, label='nbinom pmf')

n, p = 100/2, 0.5
mean, var, skew, kurt = poisson.stats(n, p, moments='mvsk')
x = np.arange(poisson.ppf(0.01, n, p),poisson.ppf(0.99, n, p))
ax[1].plot(x, poisson.pmf(x, n, p), 'bo', ms=8, label='poisson pmf')
plt.show()



import numpy as np
from scipy.stats import nbinom, poisson
import matplotlib.pyplot as plt
plt.figure(1)
plt.subplot(211)
n, p = 50,
mean, var, skew, kurt = nbinom.stats(n, p, moments='mvsk')
x = np.arange(nbinom.ppf(0.01, n, p),nbinom.ppf(0.99, n, p))
plt.plot(x, nbinom.pmf(x, n, p), 'bo', ms=8, label='nbinom pmf')
plt.vlines(mean, 0, nbinom.pmf(mean,n,p), colors='k', linestyles='-', lw=1,label='frozen pmf')
plt.xlim(0,100)
plt.ylim(0, .08)

plt.subplot(212)
n, p = 50, 0.5
mean, var, skew, kurt = poisson.stats(n, p, moments='mvsk')
x = np.arange(poisson.ppf(0.01, n, p),poisson.ppf(0.99, n, p))
plt.plot(x, poisson.pmf(x, n, p), 'bo', ms=8, label='poisson pmf')
plt.vlines(mean, 0, poisson.pmf(mean,n,p), colors='k', linestyles='-', lw=1,label='frozen pmf')
plt.xlim(0,100)
plt.ylim(0, .08)

plt.show()
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
plt.figure(1)
mu=100 # coverage
p=norm(mu)
mean, var, skew, kurt = p.stats(moments='mvsk')
x = np.arange(p.ppf(0.01),p.ppf(0.99))
plt.plot(x, p.pdf(x), 'bo', ms=8, label='p pdf')
plt.vlines(mean, 0, p.pdf(mean), colors='k', linestyles='-', lw=1,label='frozen pmf')
plt.xlim(0,200)
plt.ylim(0, .08)
