from sarray import Sarray
from logger import logger
import sys
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

fle = open('dna.50MB', 'r')
patt = fle.read()
sze = 10000000
patt = patt[0:sze]
srr = Sarray(patt[0:sze], True)
mmn = minimizer(patt[0:sze], lambda x, y: comp2(x,y,srr.rank))

w = 100
Ks = []
indx = {}
min_kpot = 4
max_kpot = 9
lgg = logger(max_kpot - min_kpot)
for kpot in range(min_kpot,max_kpot):
  lgg.log('Building K index: ', kpot)
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
runs = 10000
stats = {}
prob_del = 0.001
prob_ins = 0.001
prob_mut = 0.001
offset = 10

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
cnt = {}
for k in Ks:
  cnt[k] = {'FP':0, 'TP':0, 'FN':0}
lgg = logger(runs / 100)
while runs > 0:
  runs -= 1
  if runs % 100 == 0:
    lgg.log('Remaining runs:', runs)
  stt = random.randint(0, sze - w - Ks[-1] - offset)
  sub_patt = patt[stt:(stt + w + Ks[-1] + offset)]
  sub_patt = mutate(sub_patt, prob_del, prob_ins, prob_mut)
  if len(sub_patt) < (w + Ks[-1]):
    continue
  srr = Sarray(sub_patt)
  mmn = minimizer(sub_patt, lambda x, y: comp2(x,y,srr.rank))
  for k in Ks:
    for stt2 in range(0, len(sub_patt) - k - w):
      mnm = mmn.get_minimizer(k, w, stt2)
      if mnm in indx: # Positive
        tr = False
        for mnm2 in indx[mnm]:
          if mnm2 >= stt and mnm2 <= stt + k + w:
            tr = True
        if tr:
          cnt[k]['TP'] += 1
        else:
          cnt[k]['FP'] += 1
      else: # Negative
        cnt[k]['FN'] += 1
print(cnt)