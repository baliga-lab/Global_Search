#!/usr/bin/env python3

"""
find_files_test.py - Unit tests for the globalsearch.rnaseq.find_files module
"""

import unittest
import os, sys
import fs
import globalsearch.rnaseq.find_files as find_files


class FindFilesTest(unittest.TestCase):

    def __make_input_folder(self, num_samples, dir_levels=0,
                            pattern="R%d_%d.fq.gz"):
        parent_dir = "inputdata"
        for level in range(dir_levels):
            parent_dir = os.path.join(parent_dir, "d%d" % level)

        for i in range(1, num_samples + 1):
            dir = "%s/R%d" % (parent_dir, i)
            self.mem_fs.makedirs(dir)
            for j in range(2):
                fname = pattern % (i, j + 1)
                self.mem_fs.touch(fs.path.combine(dir, fname))

    def __make_include_file(self, path, includes):
        self.mem_fs.writetext(path, '\n'.join(includes))


    def setUp(self):
        self.mem_fs = fs.open_fs("mem://")
        self.mem_fs.makedir("inputdata")

    def tearDown(self):
        self.mem_fs.close()

    def test_rnaseq_data_folder_list_scandir(self):
        self.__make_input_folder(2, dir_levels=2)
        config = {'input_dir': '/inputdata'}
        result = find_files.rnaseq_data_folder_list(config, filesys=self.mem_fs)
        self.assertEqual(result, ['d0/d1/R1', 'd0/d1/R2'])

    def test_rnaseq_data_folder_list_includes(self):
        self.__make_input_folder(2)
        config = {'input_dir': '/inputdata', 'includes': ['R1', 'R2']}
        result = find_files.rnaseq_data_folder_list(config, filesys=self.mem_fs)
        self.assertEqual(result, ['R1', 'R2'])

    def test_rnaseq_data_folder_list_include_file(self):
        self.__make_input_folder(2)
        self.__make_include_file("/my_include.txt", ['R1', 'R2'])
        config = {'input_dir': '/inputdata', 'include_file': '/my_include.txt'}
        result = find_files.rnaseq_data_folder_list(config, filesys=self.mem_fs)
        self.assertEqual(result, ['R1', 'R2'])

    def test_rnaseq_data_folder_list_includes_and_include_file(self):
        self.__make_input_folder(2)
        self.__make_include_file("/my_include.txt", ['R1'])
        config = {
            'input_dir': '/inputdata',
            'includes': ['R2'],
            'include_file': '/my_include.txt'
        }
        result = find_files.rnaseq_data_folder_list(config, filesys=self.mem_fs)
        self.assertEqual(sorted(result), ['R1', 'R2'])

    def test_find_fastq_files(self):
        self.__make_input_folder(2)
        result = find_files.find_fastq_files("/inputdata/R1", ["*_{{readnum}}.fq.*"], filesys=self.mem_fs)
        self.assertEqual(result, [('/inputdata/R1/R1_1.fq.gz', '/inputdata/R1/R1_2.fq.gz')])

    def test_find_fastq_files_multi_pattern(self):
        self.__make_input_folder(2)
        result = find_files.find_fastq_files("/inputdata/R1",
                                             ["*_{{readnum}}.fq.*", "*_{{readnum}}.fastq.*"],
                                             filesys=self.mem_fs)
        self.assertEqual(result, [('/inputdata/R1/R1_1.fq.gz', '/inputdata/R1/R1_2.fq.gz')])

    #
    # TEST FOR SINGLE READ SETUPS
    #
    def __make_input_folder_single(self, folder_name, samples, pattern="%s_%s_1.fastq.gz"):
        for sample in samples:
            dir = "inputdata/%s" % folder_name
            self.mem_fs.makedir(dir)
            fname = pattern % (folder_name, sample)
            self.mem_fs.touch(fs.path.combine(dir, fname))

    def test_find_fastq_files_single(self):
        # Find single reads
        self.__make_input_folder_single("SRR401853", ["GSM818457_Control"])
        result = find_files.find_fastq_files("/inputdata/SRR401853", ["*_{{readnum}}.fastq.gz"], filesys=self.mem_fs)
        self.assertEqual(result, [('/inputdata/SRR401853/SRR401853_GSM818457_Control_1.fastq.gz', None)])

    def test_find_fastq_files_single_no_wildcard(self):
        # Find single reads
        self.__make_input_folder_single("SRR401853", ["GSM818457_Control"])
        result = find_files.find_fastq_files("/inputdata/SRR401853", ["*.fastq.gz"], filesys=self.mem_fs)
        self.assertEqual(result, [('/inputdata/SRR401853/SRR401853_GSM818457_Control_1.fastq.gz', None)])

if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(FindFilesTest))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
