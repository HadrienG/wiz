#!/usr/bin/env python

import os
import platform
import subprocess


def max_cpus():
    """return max cpu on the machine
    """
    return os.cpu_count()


def max_mem():
    """return max memory on the machine
    """
    if platform.system() == 'Darwin':
        hostinfo = subprocess.Popen(["hostinfo"], stdout=subprocess.PIPE)
        mem = subprocess.check_output(
            ('grep', 'memory'), stdin=hostinfo.stdout)
        hostinfo.wait()
        max_mem = float(mem.split()[3].decode())
        return max_mem
    else if platform.system() == 'Linux':
        mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        max_mem = mem_bytes/(1024.**3)
        return max_mem
