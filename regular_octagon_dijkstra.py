

# Dijkstra's algorithm on a weighted graph (the edges are weighted).
# Find the least cost path from the initial node to the final node.

nodes = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I')
distances = {
    "A": {"B": 13, "I": 12, "H": 6},
    "B": {"A": 13, "I": 4, "C": 5}, 
    "C": {"B": 5, "I": 9, "D": 14},
    "D": {"C": 14, "I": 1, "E": 8},
    "E": {"D": 8, "I": 10, "F": 15},
    "F": {"E": 15, "I": 2, "G": 7},
    "G": {"F": 7, "I": 11, "H": 16},
    "H": {"G": 16, "I": 3, "A": 6},
    "I": {"A": 12, "B": 4, "C": 9, "D": 1, "E": 10, "F": 2, "G": 11, "H": 3} 
    }

unvisited = {node: None for node in nodes}   # a dictionary with unvisited vertices  
visited = {}
loc = "C"                                    # initial starting point 
loc_dist = 0                                 # initial distance 
unvisited[loc] = loc_dist

while True:
    for i, x in distances[loc].items():
        if i not in unvisited:
            continue
        new_dist = loc_dist + x
        if unvisited[i] is None or unvisited[i] > new_dist:
            unvisited[i] = new_dist
    visited[loc] = loc_dist                  # mark the visited location with the minimal distance from loc
    del unvisited[loc]
    if not unvisited:
        break
    candidates = [node for node in unvisited.items() if node[1]]
    loc, loc_dist = sorted(candidates, key = lambda x: x[1])[0]

print(visited)
 



