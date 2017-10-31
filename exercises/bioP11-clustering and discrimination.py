
#### Kmeans Clustering


from numpy import array, random
from random import sample

def kMeans(data, k, centers=None):
    if centers is None:
        centers = array( sample(list(data), k) ) 
        change = 1.0
    while change > 1e-8:
        clusters = [[] for x in range(k)]
        for vector in data:
            diffs = centers - vector
            dists = (diffs * diffs).sum(axis=1)
            closest = dists.argmin()
            clusters[closest].append(vector)
            change = 0
        for i, cluster in enumerate(clusters):
            cluster = array(cluster)
            center = cluster.sum(axis=0)/len(cluster)
            diff = center - centers[i]
            change += (diff * diff).sum()
            centers[i] = center
    return centers, clusters

testDataA = random.random((1000,2)) # No clumps
centers, clusters = kMeans(testDataA, 3)

print('centers', centers)
print('clusters', clusters)

from numpy import random, vstack

testDataB1 = random.normal(0.0, 2.0, (100,2))
testDataB2 = random.normal(7.0, 2.0, (100,2))
testDataB = vstack([testDataB1, testDataB2]) # Two clumps
centers, clusters = kMeans(testDataB, 2)

colors = ['#FF0000','#00FF00','#0000FF','#FFFF00','#00FFFF','#FF00FF']
for i, cluster in enumerate(clusters):
    x, y = zip(*cluster)
    color = colors[i % len(colors)]
    pyplot.scatter(x, y, c=color, marker='o')
    x, y = zip(*centers)


#pyplot.scatter(x, y, s=40, c='black', marker='o')
#pyplot.show()



####### PCA  - Pricipal component analysis

from numpy import cov, linalg, sqrt, zeros, ones, diag

def principalComponentAnalysis(data, n=2):
    samples, features = data.shape
    meanVec = data.mean(axis=0)
    dataC = (data - meanVec).T
    covar = cov(dataC)
    evals, evecs = linalg.eig(covar)
    indices = evals.argsort()[::-1]
    evecs = evecs[:,indices]
    basis = evecs[:,:n]
    energy = evals[:n].sum()
# norm wrt to variance
#sd = sqrt(diag(covar))
#zscores = dataC.T / sd
    return basis, energy


from numpy import random, dot, array, outer

def extractPrincipalComponent(data, precision=1e-9):
    samples, features = data.shape
    meanVec = data.mean(axis=0)
    dataC = data - meanVec
    pc1 = random.random(features)
    pc0 = pc1 - 1.0
    while abs((pc0-pc1).sum()) > precision:
            t = zeros(features)
            for datum in dataC:
                    t += dot(datum, pc1) * datum
                    pc1 = t / sqrt(dot(t,t))
                    pc0 = pc1
            return pc1


testData = random.normal(0.0, 2.0, (100,2))
shear = array([[2,1],[1,0]])
testData = dot(testData, shear)


pc1 = extractPrincipalComponent(testData)
print('Quick PC1:', pc1)

basis, energy = principalComponentAnalysis(testData, n=2)
print('Full PCA:', basis, energy)



from matplotlib import pyplot
x,y = zip(*testData)
pyplot.scatter(x, y, s=20, c='#F0F0F0', marker='o')
x,y = zip(-10*pc1, 10*pc1)
pyplot.plot(x, y)


transformed = dot(testData, basis)
x,y = zip(*transformed)
pyplot.scatter(x, y, s=10, c='#000000', marker='^')
pyplot.show()



############ LDA - linear discriminant analysis


def twoClassLda(dataA, dataB):
    meanA = dataA.mean(axis=0)
    meanB = dataB.mean(axis=0)
    covA = cov(dataA.T)
    covB = cov(dataB.T)
    nA = len(dataA)-1.0
    nB = len(dataB)-1.0
    scatterWithin = nA * covA + nB * covB
    scatterBetween = meanA - meanB
    discrim = dot(linalg.inv(scatterWithin),scatterBetween)
    transfA = dot(dataA, discrim.T)
    transfB = dot(dataB, discrim.T)
    divide = dot(discrim,(meanA+meanB))/2.0
    return transfA, transfB, divide


testData1 = random.normal(0.0, 2.0, (100,2)) + array([-10.0,5.0])
testData2 = random.normal(0.0, 6.0, (100,2))

from matplotlib import pyplot
x, y = zip(*testData1)
pyplot.scatter(x, y, s=25, c='#404040', marker='o')
x, y = zip(*testData2)
pyplot.scatter(x, y, s=25, c='#7fffd4', marker='^')
pyplot.show()



proj1, proj2, div = twoClassLda(testData1, testData2)
print(div)
x = proj1
y = [0.5] * len(x)
pyplot.scatter(x, y, s=35, c='#404040', marker='o')
x = proj2
y = [-0.5] * len(x)
pyplot.scatter(x, y, s=35, c='#7fffd4', marker='^')
pyplot.plot((div, div), (1.0, -1.0))
pyplot.show()