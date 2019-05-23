from sarray import Sarray
class minimizer:
  def __init__(self, pat, get_mnm):
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
patt = patt[0:sze]
sze = 10000000
srr = Sarray(patt[0:sze])
mmn = minimizer(patt[0:sze], lambda x, y: comp2(x,y,srr.rank))

w = 100
Ks = []
indx = {}
for kpot in range(4,7):
  k = 2 ** kpot
  Ks.append(k)
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

import random
runs = 1000
stats = {}
prob_del = 0.001
prob_ins = 0.001
prob_mut = 0.001
offset = 100

def mutate(patt, prob_del, prob_ins, prob_mut):
  i = 0
  ret = ''
  while i < len(patt):
    p = random.uniform(0, 1)
    if p < prob_del:
      i += 1
    elif p < prob_del + prob_ins:
      ret += 'ACTG'[random.randint(0, 3)]
    elif p < prob_del + prob_ins + prob_mut:
      ret += 'ACTG'[random.randint(0, 3)]
      i += 1
    else:
      ret += patt[i]
      i += 1
  return ret

while runs > 0:
  runs -= 1
  stt = random.randint(0, sze - w - Ks[-1] - offset)
  sub_patt = patt[stt:(stt + w + Ks[-1] + offset)]
  sub_patt = mutate(sub_patt, prob_del, prob_ins, prob_mut)
  srr = Sarray(sub_patt)
  mmn = minimizer(sub_patt, lambda x, y: comp2(x,y,srr.rank))
  for k in Ks:
    mnm = mmn.get_minimizer(k, w, stt)
    if mnm in indx: # Positive
      print('pos')
    else: # Negative
      print('neg')