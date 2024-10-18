#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 05:23:33 2018

@author: rob
"""

import sys
try:
    import wrap
    import ffea_script
    import ffea_rod
except ImportError:
    from ffeatools import wrap
    from ffeatools import ffea_script
    from ffeatools import ffea_rod

def main():

    try:
        #wrap.wrap_process("ffeadev", ["symmetry_test_stretch_only.ffea"])
        wrap.wrap_process("../../../../src/ffea", ["symmetry_test_bend_only.ffea"])
        wrap.wrap_process("../../../../src/ffea", ["symmetry_test_stretch_only.ffea"])
    except OSError:
        raise
    
    bend_script = ffea_script.ffea_script("symmetry_test_bend_only.ffea")
    bend_rod = bend_script.rod[0]
    bend_rod.set_avg_energies()
    bend_analysis = ffea_rod.anal_rod(bend_rod)
    bend_test_result = bend_analysis.do_bend_symmetry_test()
    
    stretch_script = ffea_script.ffea_script("symmetry_test_stretch_only.ffea")
    stretch_rod = stretch_script.rod[0]
    stretch_rod.set_avg_energies()
    stretch_analysis = ffea_rod.anal_rod(stretch_rod)
    stretch_test_result = stretch_analysis.do_stretch_symmetry_test()
    
    if bend_test_result == False or stretch_test_result == False:
        return 1

    return 0



if __name__ == "__main__":
    error_status = main()
    sys.exit(error_status)
