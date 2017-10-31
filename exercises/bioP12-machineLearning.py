 ############# kNN  - k Nearest neighbours

def getFeatureDistance(vector1, vector2):
    distance = 0.0
    for a, b in zip(vector1, vector2):
        delta = a-b
        distance += delta * delta
        return distance



def kNearestNeighbour(knowns, query, k=7):
    if k >= len(knowns):
        raise Exception('Length of training data must be larger than k')
    dists = []
    for vector, cat in knowns[:k]:
        dist = getFeatureDistance(vector, query)
        dists.append( (dist, cat) ) # Vector could be included
        dists.sort()
        closest = dists[:k]
        counts = {}
    for dist, cat in closest:
        counts[cat] = counts.get(cat, 0) + 1
        bestCount = max(counts.values())
        bestCats = [cat for cat in counts if counts[cat] == bestCount]
    for dist, cat in closest:
        if cat in bestCats:
            return cat


knownClasses = [((1.0, 0.0, 0.0), 'warm'), # red
((0.0, 1.0, 0.0), 'cool'), # green
((0.0, 0.0, 1.0), 'cool'), # blue
((0.0, 1.0, 1.0), 'cool'), # cyan
((1.0, 1.0, 0.0), 'warm'), # yellow
((1.0, 0.0, 1.0), 'warm'), # magenta
((0.0, 0.0, 0.0), 'cool'), # black
((0.5, 0.5, 0.5), 'cool'), # grey
((1.0, 1.0, 1.0), 'cool'), # white
((1.0, 1.0, 0.5), 'warm'), # light yellow
((0.5, 0.0, 0.0), 'warm'), # maroon
((1.0, 0.5, 0.5), 'warm'), # pink
]

result = kNearestNeighbour(knownClasses, (0.7,0.7,0.2), k=3)

print('Colour class:', result)