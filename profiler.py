# profiler.py
import time
import resource
from collections import OrderedDict

class StepProfiler:
    def __init__(self):
        self.results = OrderedDict()
        self._t0 = time.perf_counter()

    def mark(self, label):
        """Record time and memory at this point."""
        t_now = time.perf_counter()
        elapsed = t_now - self._t0
        mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        self.results[label] = {
            "time_sec": elapsed,
            "mem_mb": mem_mb,
        }
        self._t0 = t_now

    def summary(self, sort_by="time_sec"):
        print("\n==================== PROFILING SUMMARY ====================")
        print(f"{'Step':40} | {'Time (s)':>10} | {'Memory (MB)':>10}")
        print("-" * 70)
        for step, d in self.results.items():
            print(f"{step:40} | {d['time_sec']:10.2f} | {d['mem_mb']:10.1f}")
        print("===========================================================\n")

    def to_csv(self, path):
        import csv
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["step", "time_sec", "mem_mb"])
            for step, d in self.results.items():
                w.writerow([step, d["time_sec"], d["mem_mb"]])
        print(f"Profiler saved → {path}")
