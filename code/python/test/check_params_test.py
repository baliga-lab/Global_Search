#!/usr/bin/env python3

"""
check_params_test.py - Unit tests for the globalsearch.control.gs_prepare module
"""

import unittest
import os, sys
import fs

from globalsearch.control.gs_prepare import check_star_options

class CheckParamsTest(unittest.TestCase):
    def test_check_star_options_out_samattrs_special_success(self):
        """check the special attributes of outSAMattributes"""
        check_star_options({'outSAMattributes': ['Standard']})
        check_star_options({'outSAMattributes': ['None']})
        check_star_options({'outSAMattributes': ['All']})

    def test_check_star_options_out_samattrs_single_fail(self):
        """check error on non existing value"""
        self.assertRaises(ValueError, check_star_options, {'outSAMattributes': ['notexists']})

    def test_check_star_options_out_samattrs_multi_fail(self):
        """check error when special attribute is used in multi context"""
        self.assertRaises(ValueError, check_star_options, {'outSAMattributes': ['All', 'None']})

    def test_check_star_options_out_samattrs_multi_success(self):
        """check successful multi context"""
        check_star_options({'outSAMattributes': ['NH', 'HI']})

    def test_check_star_options_out_samattrs_wrong_type(self):
        """check wrong type for outSAMattributes"""
        self.assertRaises(TypeError, check_star_options, {'outSAMattributes': 'NH HI'})

if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(CheckParamsTest))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
