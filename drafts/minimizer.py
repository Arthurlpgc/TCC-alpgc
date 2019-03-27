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

  def get_minimizer(self, k, w, pos):
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
    return self.pat[pos:(pos+k)]
    

def comp(x, y, pat):
  pat += pat
  if pat[x:] < pat[y:]:
    return x
  else:
    return y
    
patt = "abacabadabacaba"
mmn = minimizer(patt, 3, lambda x, y: comp(x,y,patt))
while True:
  [k, w, pos] = list(map(lambda x:eval(x), input().split(' ')))
  print(patt[pos:(pos + w)])
  print(mmn.get_minimizer(k, w, pos))