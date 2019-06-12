import sys
import datetime
import time, threading



class logger:
  def __init__(self, log_steps):
    self.tm = datetime.datetime.now()
    self.log_steps = log_steps
    self.cur_step = 0
    self.current_time = True
    # self.thr = threading.Thread(target=self.current_time_func)
    # self.thr.start()
  # def current_time_func(self):
  #   while True:
  #     if self.current_time:
  #       cur_tm = datetime.datetime.now()
  #       print('\033[38;2;200;200;0m\033[100;0H', (cur_tm - self.tm),'\033[0m\033[1A\033[1000D', end='', file=sys.stderr)
  #       sys.stderr.flush()
  #     else:
  #       break
  #     time.sleep(1)
  def log(self, *args):
    self.cur_step += 1
    if self.cur_step == self.log_steps:
      self.current_time = False
    cur_tm = datetime.datetime.now()
    print(*args, '\033[0;31mPAST:', (cur_tm - self.tm), '\033[0;32mETA:', (cur_tm - self.tm) * (self.log_steps - float(self.cur_step)) / self.cur_step, '\033[0;0m',file=sys.stderr)