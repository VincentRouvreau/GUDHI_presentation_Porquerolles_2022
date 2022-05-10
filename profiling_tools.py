# Some time profiling tools
from timeit import default_timer as timer
from datetime import timedelta

# Some memory profiling tools
import os, psutil

class Chrono(object):
    start = 0

    def __init__(self):
        self.start = timer()
    
    def end_timer(self):
        return timedelta(seconds = timer() - self.start) / timedelta(seconds=1)
  
class Ram(object):
    start = 0

    @staticmethod
    def __get_ram_used():
        return psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2

    def __init__(self):
        self.start = Ram._Ram__get_ram_used()
    
    def get(self):
        return Ram._Ram__get_ram_used() - self.start
