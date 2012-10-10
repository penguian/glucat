# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorials.py: This file implements a simple method to run
#                        PyClical tutorials.
#
#    copyright            : (C) 2012 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

from pyclical_tutorial_utils import *

class tutorial_module:
    def __init__(self, tutorial_title, tutorial_module_name):
        self.title       = tutorial_title
        self.module_name = tutorial_module_name

tutorial_dict = dict( [
    ("0.0", tutorial_module("0.0 Notation",
                            "pyclical_tutorial_0_0_notation")),
                            
    ("0.1", tutorial_module("0.1 Index sets",
                            "pyclical_tutorial_0_1_index_sets")),
                            
    ("0.2", tutorial_module("0.2 Operations",
                            "pyclical_tutorial_0_2_operations")),
                            
    ("0.3", tutorial_module("0.3 Algebraic functions",
                            "pyclical_tutorial_0_3_functions")),
                            
    ("0.4", tutorial_module("0.4 Square root and transcendental functions.",
                            "pyclical_tutorial_0_4_transcendental")),
                            
    ("1.0", tutorial_module("1.0 Plane geometry",
                            "pyclical_tutorial_1_0_plane")),
                            
    ("1.1", tutorial_module("1.1 Complex numbers",
                            "pyclical_tutorial_1_1_complex")),
                            
    ("1.2", tutorial_module("1.2 Space geometry and vector algebra",
                            "pyclical_tutorial_1_2_space")),
                            
    ("1.3", tutorial_module("1.3 Electromagnetism and Lorentz transformations",
                            "pyclical_tutorial_1_3_lorentz")),
                            
    ("1.4", tutorial_module("1.4 The fourth dimension",
                            "pyclical_tutorial_1_4_fourth")),

    ("1.5", tutorial_module("1.5 Conformal Geometric Algebra.",
                            "pyclical_tutorial_1_5_conformal")),
                      ])

def tutorial():
    input_str = ""
    while True:
        print ""
        print_fill("Currently available PyClical tutorials:")
        print ""
        for key, tut in sorted(tutorial_dict.iteritems()):
            print_fill(tut.title)
            
        print ""
        input_str = raw_input(fill("Enter the number for the tutorial you want to run, or Q to quit:") + " ")
        if (len(input_str) > 0) and (input_str[0].upper() == 'Q'):
            break
            
        try:
            tut_nbr_str = "{:3.1f}".format(float(input_str))
            tut = tutorial_dict[tut_nbr_str]
            tut_module = __import__(tut.module_name)
            tut_module.tutorial()
        except KeyboardInterrupt:
            print_fill("The tutorial was interrupted.")
            continue
        except:
            print_fill("Please enter a valid tutorial number.")

if __name__ == "__main__":
    tutorial()
    