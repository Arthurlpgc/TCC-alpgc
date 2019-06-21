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

test_run = 1

fle = open('dna.50MB', 'r')
pattern = fle.read()
pattern = pattern.replace('\n', '')
sze = 10000000 // (1000 ** test_run)
pattern = pattern[0:sze]
srr = Sarray(pattern, True)
mmn = minimizer(pattern, lambda x, y: comp2(x,y,srr.rank))

w = 200
Ks = [16,32,64,128]

def get_index(w, Ks, patt, mmn):
  indx = {}
  lgg = logger(len(Ks))
  for k in Ks:
    lgg.log('Building K index: ', k)
    for pos in range(sze - w):
      minimizer_pos = mmn.get_minimizer_pos(k, w - k, pos)
      key = patt[minimizer_pos:(minimizer_pos + k)]
      if key not in indx:
        indx[key] = {}
      indx[patt[minimizer_pos:(minimizer_pos + k)]][minimizer_pos] = True
  return indx

indx = get_index(w, Ks, pattern, mmn)
# print(indx)
import random
runs = 10000 // (100 ** test_run)
stats = {}
prob_del = 0.01
prob_ins = 0.01
prob_mut = 0.01
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
cnt = { -1: {'SZ':0, 'FNDSZE':0, 'FND':0, 'MIS':0, 'NOT': 0} }
for k in Ks:
  cnt[k] = {'SZ':0, 'FNDSZE':0, 'FND':0, 'MIS':0, 'NOT': 0}
lgg = logger(runs / 100)

# patts = open('files/m130928_232712_42213_c100518541910000001823079209281310_s1_p0.1.subreads.fastq', 'r').read()
# patts = patts.split('\n+\n')
# patts = patts[:-1]
# for patt_long in patts:
#   patt_long = patt_long.split('\n')
#   patt = patt_long[-1]
#   patt_data = patt_long[-2]
#   print('[',len(patt), patt_data,']')

def process_overlap(overlap):
  if len(overlap) == 0:
    return {'size': 0, 'pos': 0}
  covered = 0
  size = 0
  for minim in overlap:
    size += minim[2]
    size -= min(max(covered - minim[1], 0), minim[2])
    covered = max(covered, minim[1] + minim[2])
  return {'size': size, 'pos': overlap[len(overlap) // 2][0]}

def get_overlap(finger_print, cluster_threshold = 10):
  finger_print = list(set(finger_print))
  finger_print.sort()
  best_overlap = {'size': 0, 'pos': 0}
  overlap = []
  for fp in finger_print:
    if len(overlap) == 0:
      overlap.append(fp)
    else:
      if fp[0] - overlap[-1][0] > cluster_threshold:
        overlap = process_overlap(overlap)
        if overlap['size'] > best_overlap['size']:
          best_overlap = overlap
        overlap = [fp]
      else:
        overlap.append(fp)
  overlap = process_overlap(overlap)
  if overlap['size'] > best_overlap['size']:
    best_overlap = overlap
  return best_overlap
  
patt_size = 1000
while runs > 0:
  runs -= 1
  if runs % 100 == 0:
    lgg.log('Remaining runs:', runs)
  stt = random.randint(0, sze - patt_size)
  sub_patt = pattern[stt:(stt + patt_size)]
  sub_patt = mutate(sub_patt, prob_del, prob_ins, prob_mut)
  if len(sub_patt) < (w):
    continue
  srr = Sarray(sub_patt)
  mmn = minimizer(sub_patt, lambda x, y: comp2(x,y,srr.rank))
  for k in Ks:
    finger_print = []
    for stt2 in range(0, len(sub_patt) - w):
      mnm_pos = mmn.get_minimizer_pos(k, w - k, stt2)
      mnm = mmn.pat[mnm_pos:(mnm_pos+k)]
      if mnm in indx: # Positive
        for mnm2 in indx[mnm]:
          finger_print.append((mnm2 - mnm_pos, mnm2, k))
          assert(mnm == pattern[mnm2:(mnm2+len(mnm))])
    overlap = get_overlap(finger_print)
    cnt[k]['SZ'] += overlap['size']
    if overlap['size'] == 0:
      cnt[k]['NOT'] += 1
    else:
      if abs(overlap['pos'] - stt) < 10:
        cnt[k]['FND'] += 1
        cnt[k]['FNDSZE'] += overlap['size']
      else:
        cnt[k]['MIS'] += 1

  finger_print = []
  for stt2 in range(0, len(sub_patt) - w):
    ind = len(Ks) - 1
    found = False
    while (ind >= 0) and not found:
      k = Ks[ind]
      mnm_pos = mmn.get_minimizer_pos(k, w - k, stt2)
      mnm = mmn.pat[mnm_pos:(mnm_pos+k)]
      if mnm in indx: # Positive
        for mnm2 in indx[mnm]:
          finger_print.append((mnm2 - mnm_pos, mnm2, k))
          assert(mnm == pattern[mnm2:(mnm2+len(mnm))])
        found = True
      ind -= 1  
  k=-1  
  overlap = get_overlap(finger_print)
  cnt[k]['SZ'] += overlap['size']
  if overlap['size'] == 0:
    cnt[k]['NOT'] += 1
  else:
    if abs(overlap['pos'] - stt) < 10:
      cnt[k]['FND'] += 1
    else:
      cnt[k]['MIS'] += 1
  

print(cnt)
for key in cnt:
  print(key, end=":\t")
  for key2 in cnt[key]:
    for key3 in cnt[key]:
      if key2 == key3:
        continue
      print("{}/{}".format(key2,key3), '{:.3f}'.format(cnt[key][key2] / max(cnt[key][key3], 1)), end=' | ')
  print()

