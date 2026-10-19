#!/usr/bin/env python3
"""Run a command and append its peak RSS (ru_maxrss) to the same log.

/usr/bin/time is absent on the benchmark host, so peak-RSS measurement needs
this wrapper. Declared as an input by the fastmaps rule so editing it forces a
rerun rather than reusing a stale result.
"""
import os
import sys
import time


def main():
    argv = sys.argv[1:]
    if not argv:
        print("usage: rss_wrap.py CMD [ARGS...]", file=sys.stderr)
        return 2
    t0 = time.time()
    pid = os.fork()
    if pid == 0:
        os.execvp(argv[0], argv)
    _, status, ru = os.wait4(pid, 0)
    dt = time.time() - t0
    name = os.path.basename(argv[0])
    print(f"rss_wrap: {name} exit={status} wall={dt:.1f}s "
          f"peak_rss={ru.ru_maxrss / 1024.0:.1f} MB", file=sys.stderr)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
