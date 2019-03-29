import copy 

class Sarray:
  def __init__(self, pat):
    self.sarray = []
    sze = len(pat)
    for i in range(sze):
      self.sarray.append([[ord(pat[i]),0],i])
    bucket = 1
    while bucket <= sze:
      rank = [0]*sze
      for i in range(sze):
        cr = self.sarray[i]
        rank[cr[1]] = cr[0][0]
      for i in range(sze):
        ind = self.sarray[i][1]
        nr = -1
        if ind + bucket < sze:
          nr = rank[ind + bucket]
        self.sarray[i][0][1] = nr
      self.sarray.sort()
      fr = None
      ind = 0
      for i in range(sze):
        if self.sarray[i][0] == fr:
          self.sarray[i][0][0] = ind
        else:
          fr = copy.copy(self.sarray[i][0])
          ind += 1
          self.sarray[i][0][0] = ind
      bucket += bucket
    self.rank = [0]*sze
    for i in range(sze):
      cr = self.sarray[i]
      self.rank[cr[1]] = cr[0][0]

          