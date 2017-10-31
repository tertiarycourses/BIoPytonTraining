# selecting from a normal distribution using random.normal, which we then #show as a histogram:

from matplotlib import pyplot
from numpy import random
mean = 0.0
stdDev = 1.0

for nPoints in (10, 100, 1000, 10000,100000):sample = random.normal(mean, stdDev, nPoints)


pyplot.hist(sample, bins=20, range=(-4,4), normed=True)
pyplot.show()



# mode, median, mean

values = [1,2,2,3,2,1,4,2,3,1,0]

from scipy import stats
from numpy import array
valArray = array(values, float)
mode, count = stats.mode(valArray)
print('Mode:', mode)
print('Count:', count)


from numpy import median
med = median(valArray)
print('Median:', med)

## variance and std dev

from numpy import array
valArray = array(values)
variance = valArray.var() # Biased estimate
print('Var 1', variance) # Result is 1.1736
variance = valArray.var(ddof=1) # Unbiased estimate
print('Var 2', variance)


from numpy import std, sqrt
stdDev = sqrt(variance)
stdDev = std(valArray)  # biased estimate
stdDev = valArray.std(ddof=1)  # unbiased estimate
print('Std:', stdDev)



##############  Statistical tests  ################

#The two-tailed binomial test is available as a function in SciPy:

from scipy.stats import binom_test
count, nTrials, pEvent = 530, 1000, 0.5
result = binom_test(count, nTrials, pEvent)
print('Binomial two tail', result)



#We can easily calculate a Z-score in Python,
#here taking parameters from the human height example:

from numpy import abs
mean = 1.76
stdDev = 0.075
values = array([1.8, 1.9, 2.0])
zScores = abs(values - mean)/stdDev
print('Z scores', zScores)


#Note that the T-statistic as well as the two-tailed probability
#are passed back by the function:

from scipy.stats import ttest_1samp
trueMean = 1.76
samples = array([1.752, 1.818, 1.597, 1.697, 1.644, 1.593, 1.878, 1.648, 1.819, 1.794, 1.745, 1.827])


# the two-sample T-test,
from scipy.stats import ttest_ind
samples1 = array([1.752, 1.818, 1.597, 1.697, 1.644, 1.593])
samples2 = array([1.878, 1.648, 1.819, 1.794, 1.745, 1.827])
tStat, twoTailProb = ttest_ind(samples1, samples2)
print('tStat', tStat)
print('twoTailProb', twoTailProb)


tStat, twoTailProb = ttest_1samp(samples, trueMean)
print('tStat', tStat)
print('mean', trueMean)

# Chisquare

from scipy.stats import chisquare
from scipy.stats import norm
from numpy import array
bins = array([1.65, 1.7, 1.75, 1.8, 1.85,])
obsd = array([ 14, 15, 33, 22, 8,])
mean = 1.76
std = 0.075
nObs = obsd.sum()
expd = norm.pdf(bins, mean, std)
expd *= nObs / expd.sum()
# Expected counts: 9.196, 19.576, 26.720, 23.385, 13.122

chSqStat, pValue = chisquare(obsd, expd)
print('Chi square A', chSqStat, pValue)



## Correlations

from numpy import random, cov

xVals = random.normal(0.0, 1.0, 100)
yVals1 = random.normal(0.0, 1.0, 100) # Random, independent of xVals
deltas = random.normal(0.0, 0.75, 100)
yVals2 = 0.5 + 2.0 * (xVals - deltas) # Derived from xVals

from numpy import corrcoef

r1 = corrcoef(xVals, yVals1)[0, 1] 
r2 = corrcoef(xVals, yVals2)[0, 1] 

r1
r2

### Simple linear regression

from numpy import cov, var, mean, random
xVals = random.normal(0.0, 1.0, 100)
yVals = 2.0 + -0.7 * xVals + random.normal(0.0, 0.2, 100)
grad = cov(xVals, yVals)/var(xVals, ddof=1)
yInt = mean(yVals) - grad * mean(xVals)
print('LR 1:', grad, yInt) 


from scipy.stats import linregress
from matplotlib import pyplot
grad, yInt, corrCoeff, pValue, stdErr = linregress(xVals, yVals)
print('LR 2:', grad, yInt, corrCoeff, pValue, stdErr)
xValsFit = [xVals.min(),xVals.max()]
yValsFit = [yInt + x*grad for x in xValsFit]
pyplot.plot(xVals, yVals, 'o')
pyplot.plot(xValsFit, yValsFit)
pyplot.show()

