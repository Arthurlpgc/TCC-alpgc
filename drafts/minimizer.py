from sarray import Sarray
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

fle = open('dna.200MB', 'r')
patt = fle.read()
sze = 10000000
srr = Sarray(patt[0:sze])
mmn = minimizer(patt[0:sze], 3, lambda x, y: comp2(x,y,srr.rank))

for wpot in range(2,10):
  w = 10 * (2 ** wpot)
  for kpot in range(4,7):
    indx = {}
    k = 2 ** kpot
    for pos in range(sze - k - w):
      minimizer_pos = mmn.get_minimizer_pos(k, w, pos)
      key = patt[minimizer_pos:(minimizer_pos + k)]
      if key not in indx:
        indx[key] = {}
      indx[patt[minimizer_pos:(minimizer_pos + k)]][minimizer_pos] = True
    cnt = 0
    for key in indx:
      for key2 in indx[key]:
        cnt += 1
    print(k, w, cnt, sep=', ')

cnt = 0
for key in indx:
  for key2 in indx[key]:
    cnt += 1
print(cnt)