
import numpy as np
from scipy.stats import nbinom, poisson
import matplotlib.pyplot as plt
plt.figure(1)
plt.subplot(211)
mu=100 # coverage
poissonDistribution=poisson(mu)
mean, var, skew, kurt = poissonDistribution.stats(moments='mvsk')
x = np.arange(poissonDistribution.ppf(0.01),poissonDistribution.ppf(0.99))
plt.plot(x, poissonDistribution.pmf(x), 'bo', ms=8, label='poissonDistribution pmf')
plt.vlines(mean, 0, poissonDistribution.pmf(mean), colors='k', linestyles='-', lw=1,label='frozen pmf')
plt.xlim(0,200)
plt.ylim(0, .08)

plt.subplot(212)
r=90 # controls the deviation from the poisson
# This makes the negative binomial distribution suitable as a robust alternative to the Poisson,
# which approaches the Poisson for large r, but which has larger variance than the Poisson for small r.
p=((mu*1.0)/((mu+r)*1.0))
nbinomDistribution=nbinom(r, p)
mean, var, skew, kurt = nbinomDistribution.stats(moments='mvsk')
x = np.arange(nbinomDistribution.ppf(0.01),nbinomDistribution.ppf(0.99))
plt.plot(x, nbinomDistribution.pmf(x), 'bo', ms=8, label='nbinomDistribution pmf')
plt.vlines(mean,0,nbinomDistribution.pmf(int(mean)), colors='k', linestyles='-', lw=1,label='frozen pmf')
plt.xlim(0,200)
plt.ylim(0, .08)



plt.show()
