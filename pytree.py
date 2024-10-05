import numpy as np
from collections import namedtuple

def create_points(n, d, seed=-1):
    if seed >= 0: np.random.seed(seed)
    points = 2*np.random.random((n,d)).astype(np.float64)-1
    return points

def distance(p, q):
    return np.linalg.norm(p-q)

HubPoint = namedtuple("HubPoint", ["id", "leader", "dist"])

class Hub(object):

    def __init__(self, leaders, hub_points, hub_size, representative, candidate, candidate_point, hub_parent, hub_radius, hub_sep):
        self.leaders = leaders
        self.hub_points = hub_points
        self.hub_size = hub_size
        self.representative = representative
        self.candidate = candidate
        self.candidate_point = candidate_point
        self.hub_parent = hub_parent
        self.hub_radius = hub_radius
        self.hub_sep = hub_sep
        self.active = True

        self.new_hubs = []
        self.leaves = []

    @classmethod
    def from_protected(cls, points, repr_pt):

        leaders = [0]
        hub_size = len(points)
        representative = 0
        candidate = 0
        hub_parent = -1
        hub_radius = 0.
        hub_sep = 0.

        hub_points = []

        for i in range(hub_size):
            hub_points.append(HubPoint(i, 0, distance(repr_pt, points[i])))
            if hub_points[i].dist > hub_points[candidate].dist:
                candidate = i

        hub_radius = hub_sep = hub_points[candidate].dist
        candidate_point = points[candidate]

        return cls(leaders, hub_points, hub_size, representative, candidate, candidate_point, hub_parent, hub_radius, hub_sep)

    @classmethod
    def from_points(cls, points):
        return cls.from_protected(points, points[0])

    @classmethod
    def from_hub(cls, hub_points, representative, candidate, candidate_point, parent, radius):
        return cls([representative], hub_points, len(hub_points), representative, candidate, candidate_point, parent, radius, radius)

    def size(self): return self.hub_size
    def repr(self): return self.representative
    def cand(self): return self.candidate
    def parent(self): return self.hub_parent
    def radius(self): return self.hub_radius
    def sep(self): return self.hub_sep

    def getids(self):
        for hpt in self.hub_points:
            yield hpt.id

    def is_split(self, split_ratio):
        return self.sep() <= split_ratio*self.radius()

    def get_cand_ball(self):
        return (candidate_point, candidate, hub_sep)

    def add_new_leader(self, points):
        n = self.hub_size
        new_leader = self.candidate
        self.leaders.append(new_leader)
        self.hub_sep = 0.

        for i in range(n):
            #  hpt = self.hub_points[i]
            idx, leader, dist = self.hub_points[i]
            new_leader_dist = distance(self.candidate_point, points[idx])
            if new_leader_dist < dist:
                dist = new_leader_dist
                leader = new_leader
            if dist > self.hub_sep:
                self.candidate = idx
                self.hub_sep = dist
            self.hub_points[i] = HubPoint(idx, leader, dist)

        self.candidate_point = points[self.candidate]

    def split_leaders(self, points):
        n = self.hub_size
        for leader in self.leaders:
            relcand = 0
            new_hub_points = []
            for i in range(n):
                if self.hub_points[i].leader == leader:
                    new_hub_points.append(self.hub_points[i])
                    if new_hub_points[-1].dist > new_hub_points[relcand].dist:
                        relcand = len(new_hub_points)-1
            new_radius = new_hub_points[relcand].dist
            new_candidate = new_hub_points[relcand].id
            new_candidate_point = points[new_hub_points[relcand].id]
            self.new_hubs.append(Hub.from_hub(new_hub_points, leader, new_candidate, new_candidate_point, -1, new_radius))

    def find_leaves(self, min_hub_size):
        updated_new_hubs = []
        for new_hub in self.new_hubs:
            if new_hub.size() <= min_hub_size:
                for p in new_hub.getids():
                    self.leaves.append(p)
            else:
                updated_new_hubs.append(new_hub)
        self.new_hubs[:] = updated_new_hubs

    def add_hub_vertex(self, tree):
        assert self.active
        self.hub_vertex = len(tree)
        vtx = {'repr' : self.repr(), 'radius' : self.radius(), 'parent' : self.parent(), 'children' : []}
        tree[self.hub_vertex] = vtx
        if self.parent() >= 0: tree[self.parent()]['children'].append(self.hub_vertex)
        return self.hub_vertex

    def update_tree(self, tree, next_hubs):
        assert self.active
        self.add_hub_vertex(tree)
        self.active = False
        for new_hub in self.new_hubs:
            new_hub.hub_parent = self.hub_vertex
            next_hubs.append(new_hub)
        for leaf in self.leaves:
            vtx = len(tree)
            tree[vtx] = {'repr' : leaf, 'radius' : 0., 'parent' : self.hub_vertex, 'children' : []}
        return len(self.leaves)

class CoverTree(object):

    def __init__(self, points, globids=None):
        self.points = points
        self.globids = globids if globids else list(range(len(points)))
        self.tree = dict()

    def build(self, split_ratio, min_hub_size):
        size = len(self.points)
        hubs = [Hub.from_points(self.points)]
        it = 1
        leaf_count = 0

        while True:
            next_hubs = []
            num_hubs = len(hubs)
            for i in range(num_hubs):
                hubs[i].add_new_leader(self.points)
                if hubs[i].is_split(split_ratio):
                    hubs[i].split_leaders(points)
                    hubs[i].find_leaves(min_hub_size)
                    leaf_count += hubs[i].update_tree(self.tree, next_hubs)
                else:
                    next_hubs.append(hubs[i])
            hubs, next_hubs = next_hubs, hubs
            avg_hub_size = (size-leaf_count)/len(hubs) if len(hubs) else float('inf')
            leaf_percent = (100.*leaf_count)/size
            print(f"[msg::build] {leaf_percent:.2f} percent leaves reached [iter={it},vertices={len(self.tree)},hubs={len(hubs)},avg_hub_size={avg_hub_size:.3f}]")
            it += 1
            if leaf_count >= size:
                break

#  points = create_points(1000, 2)
fv = np.fromfile("points", dtype="int32")
d = fv.view(np.int32)[0]
new = fv.reshape(-1,d+1)[:,1:]
points = new.view(np.float32)
ctree = CoverTree(points)
ctree.build(0.5, 2)
