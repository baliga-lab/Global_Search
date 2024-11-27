#!/usr/bin/env python3

"""
trim_galore_test.py - Unit tests for the globalsearch.rnaseq.trim_galore module
"""

import unittest
import os, sys
import fs
import globalsearch.rnaseq.trim_galore as trim_galore


class TrimGaloreTest(unittest.TestCase):

    def __make_paired_folder(self, num_samples, pattern="R%d_%d_val_%d.fq.gz"):
        for i in range(1, num_samples + 1):
            dir = "inputdata/R%d" % i
            self.mem_fs.makedir(dir)
            for j in range(2):
                fname = pattern % (i, j + 1, j + 1)
                self.mem_fs.touch(fs.path.combine(dir, fname))

    def __make_single_folder(self, num_samples=1, pattern="R%d_%d_trimmed.fq.gz", readnum=1):
        for i in range(1, num_samples + 1):
            dir = "inputdata/R%d" % i
            self.mem_fs.makedir(dir)
            fname = pattern % (i, readnum)
            self.mem_fs.touch(fs.path.combine(dir, fname))


    def setUp(self):
        self.mem_fs = fs.open_fs("mem://")
        self.mem_fs.makedir("inputdata")

    def tearDown(self):
        self.mem_fs.close()

    def test_collect_trimmed_data_paired_end(self):
        """setup: 1 first pair file, 1 second pair file"""
        self.__make_paired_folder(2)
        first_pair_group, second_pair_group = trim_galore.collect_trimmed_data(
            "/inputdata/R1", "gz", is_paired_end=True, rootfs=self.mem_fs)
        self.assertEqual("/inputdata/R1/R1_1_val_1.fq.gz", first_pair_group)
        self.assertEqual("/inputdata/R1/R1_2_val_2.fq.gz", second_pair_group)

    def test_collect_trimmed_data_single_end(self):
        """setup: 1 first pair file"""
        self.__make_single_folder()
        first_pair_group, second_pair_group = trim_galore.collect_trimmed_data(
            "/inputdata/R1", "gz", is_paired_end=False, rootfs=self.mem_fs)
        self.assertEqual("/inputdata/R1/R1_1_trimmed.fq.gz", first_pair_group)
        self.assertEqual("", second_pair_group)

    def test_collect_trimmed_data_single_end_non_std(self):
        """setup: 1 first pair file"""
        self.__make_single_folder(readnum=3)
        first_pair_group, second_pair_group = trim_galore.collect_trimmed_data(
            "/inputdata/R1", "gz", is_paired_end=False, rootfs=self.mem_fs)
        self.assertEqual("/inputdata/R1/R1_3_trimmed.fq.gz", first_pair_group)
        self.assertEqual("", second_pair_group)


if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(TrimGaloreTest))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
