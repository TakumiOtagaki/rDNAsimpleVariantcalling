# define the segment tree

# add +1 i,j --> array[i:j]++
# query i,j --> sum(array[i:j])
class SegmentTree:
    def __init__(self, n):
        self.n = n
        self.tree = [0] * (2*n)
    def add(self, i, j):
        i += self.n
        j += self.n
        while i < j:
            if i % 2 == 1:
                self.tree[i] += 1
                i += 1
            if j % 2 == 1:
                j -= 1
                self.tree[j] += 1
            i //= 2
            j //= 2
    def query(self, i):
        i += self.n
        res = 0
        while i > 0:
            res += self.tree[i]
            i //= 2
        return res
    def print(self):
        print(self.tree)
        
