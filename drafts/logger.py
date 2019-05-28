import sys
import datetime
class logger:
  def __init__(self, log_steps):
    self.tm = datetime.datetime.now()
    self.log_steps = log_steps
    self.cur_step = 0

  def log(self, *args):
    self.cur_step += 1
    cur_tm = datetime.datetime.now()
    print(*args, '\033[0;31mPAST:', (cur_tm - self.tm), '\033[0;32mETA:', (cur_tm - self.tm) * (self.log_steps - float(self.cur_step)) / self.cur_step, '\033[0;0m',file=sys.stderr)