import graphviz

def read_tree_file(fname):
    f = open(fname, "r")
    header = f.readline().rstrip()
    assert header == "vertex_id\tparent_id\tpoint_id\tradius"
    vertices = []
    for line in f.readlines():
        tokens = line.rstrip().split("\t")
        vtx, par, pt = (int(v) for v in tokens[:3])
        rad = float(tokens[3])
        vertices.append((vtx,par,pt,rad))
    f.close()

    adjlist = [list() for i in range(len(vertices))]
    for vertex, parent, point, radius in vertices:
        if parent == -1: continue
        adjlist[parent].append(vertex)

    return vertices, adjlist


vertices, adjlist = read_tree_file("tree.tsv")

g = graphviz.Digraph()
#  g.attr(rankdir="LR")

for vertex, parent, point, radius in vertices:
    label = f"{radius:.3f}" if adjlist[vertex] else f"{point}"
    color = "lightblue" if adjlist[vertex] else "green"
    #  label = f"ball({point}, {radius:.3f})" if adjlist[vertex] else f"leaf({point})"
    g.node(str(vertex), label=label, color=color)

for source in range(len(vertices)):
    neighbors = adjlist[source]
    for target in neighbors:
        color = None if vertices[target][2] != vertices[source][2] else "red"
        g.edge(str(source), str(target), color=color)

g.render(filename="output")
