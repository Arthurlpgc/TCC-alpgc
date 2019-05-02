from sarray import Sarray
from random import randint
class minimizer:
  def __init__(self, pat, k_max, get_mnm):
    self.k_max = k_max
    self.sze = len(pat)
    self.pat = pat
    self.get_mnm = get_mnm
    self.minimizer_st = [[x for x in range(self.sze)]]
    bucket = 1
    while bucket <= self.sze:
      next_lvl = [-1] * self.sze
      ind = 0
      while ind + bucket < self.sze:
        next_lvl[ind] = get_mnm(self.minimizer_st[-1][ind], self.minimizer_st[-1][ind + bucket])
        ind += 1
      self.minimizer_st.append(next_lvl)
      bucket += bucket

  def get_minimizer_pos(self, k, w, pos):
    if k + w < self.k_max:
      return -1
    bucket = 1
    ind = 0
    while bucket <= w:
      bucket += bucket
      ind += 1
    bucket //= 2
    ind -= 1
    pos = self.get_mnm(self.minimizer_st[ind][pos], self.minimizer_st[ind][pos + w - 1 - bucket])
    return pos
  def get_minimizer(self, k, w, pos):
    pos = self.get_minimizer_pos(k, w, pos)
    return self.pat[pos:(pos+k)]
    

def comp(x, y, pat):
  pat += pat
  if pat[x:] < pat[y:]:
    return x
  else:
    return y
    
def comp2(x, y, inv_sarray):
  if inv_sarray[x] < inv_sarray[y]:
    return x
  else:
    return y

patt = ""
dics = {}
i = 0
while i < 100000:
  i += 1
  patt += "actg"[randint(0,3)]
print(patt)
srr = Sarray(patt)
k = 8
mmn = minimizer(patt, k, lambda x, y: comp2(x,y,srr.rank))
while True:
  # [k, w, pos] = list(map(lambda x:eval(x), input().split(' ')))
  w = 100
  for pos in range(0, 100000 - k - w):
    dics[mmn.get_minimizer_pos(k, w, pos)] = True
  break
print(len(dics.keys()))