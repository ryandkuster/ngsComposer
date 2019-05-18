import unittest
import sys
import os
sys.path.append("../tools/")
import crinoid
import rotifer
import scallop


'''
usage:
python3 -m unittest test_tools.py 
'''


class TestCrinoid(unittest.TestCase):

    def test_dict_maker(self):
        max_len = 1
        score_ref_dt, score_dt, base_ref_dt, base_dt = crinoid.dict_maker(max_len)
        expected = {'!': '0', '"': '1', '#': '2', '$': '3', '%': '4',
                  '&': '5', "'": '6', '(': '7', ')': '8', '*': '9',
                  '+': '10', ',': '11', '-': '12', '.': '13', '/': '14',
                  '0': '15', '1': '16', '2': '17', '3': '18', '4': '19',
                  '5': '20', '6': '21', '7': '22', '8': '23', '9': '24',
                  ':': '25', ';': '26', '<': '27', '=': '28', '>': '29',
                  '?': '30', '@': '31', 'A': '32', 'B': '33', 'C': '34',
                  'D': '35', 'E': '36', 'F': '37', 'G': '38', 'H': '39',
                  'I': '40'}
        self.assertEqual(score_ref_dt, expected)

    def test_output_prep(self):
        dt1 = {0: {'!': 5, '"': 5, '#': 5}, 1: {'!': 3, '"': 3, '#': 3}}
        mx = [[0, 0, 0], [0, 0, 0]]
        dt2 = {'!': '0', '"': '1', '#': '2'}
        col = 0
        result = crinoid.output_prep(dt1, mx, dt2, col)
        expected = [[5, 5, 5], [3, 3, 3]]
        self.assertEqual(result, expected)


class TestRotifer(unittest.TestCase):

    def test_rotifer1(self):
        in1 = './data/input/R1_1000.fastq'
        in2 = './data/input/R2_1000.fastq'
        R1_bases_ls = ['AA', 'CC', 'GG', 'TT']
        R2_bases_ls = None
        pe_1 = './data/output/pe.R1_single.fastq'
        pe_2 = './data/output/pe.R2_single.fastq'
        se_1 = './data/output/se.R1_single.fastq'
        se_2 = './data/output/se.R2_single.fastq'
        rotifer.rotifer(R1_bases_ls, R2_bases_ls, in1, in2, pe_1, pe_2, se_1, se_2)
        with open(pe_1) as f1, open('./data/expected/test_rotifer1.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove(pe_1)
        os.remove(pe_2)
        os.remove(se_1)
        os.remove(se_2)

    def test_rotifer2(self):
        in1 = './data/input/R1_1000.fastq'
        in2 = './data/input/R2_1000.fastq'
        R1_bases_ls = ['AA', 'CC', 'GG', 'TT']
        R2_bases_ls = ['AA', 'CC', 'GG', 'TT']
        pe_1 = './data/output/pe.R1_single.fastq'
        pe_2 = './data/output/pe.R2_single.fastq'
        se_1 = './data/output/se.R1_single.fastq'
        se_2 = './data/output/se.R2_single.fastq'
        rotifer.rotifer(R1_bases_ls, R2_bases_ls, in1, in2, pe_1, pe_2, se_1, se_2)
        with open(pe_2) as f1, open('./data/expected/test_rotifer2.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove(pe_1)
        os.remove(pe_2)
        os.remove(se_1)
        os.remove(se_2)

    def test_rotifer_single(self):
        in1 = './data/input/R1_1000.fastq'
        se_1 = './data/output/se.R1_single.fastq'
        R1_bases_ls = ['GG', 'TT']
        rotifer.rotifer_single(R1_bases_ls, in1, se_1)
        with open(se_1) as f1, open('./data/expected/test_rotifer_single.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove(se_1)

    def test_rotifer_test1(self):
        line = 'ACGT'
        bases_ls = ['AG', 'AC']
        result = rotifer.rotifer_test(line, bases_ls)
        self.assertEqual(result, False)

    def test_rotifer_test2(self):
        line = 'ACGT'
        bases_ls = ['GG', 'TC']
        result = rotifer.rotifer_test(line, bases_ls)
        self.assertEqual(result, True)


class TestScallop(unittest.TestCase):

    def test_scallop1(self):
        in1 = './data/input/R1_single.fastq'
        front_trim = 20
        back_trim = None
        out1 = './data/output/trimmed.R1_single.fastq'
        scallop.scallop(in1, front_trim, back_trim, out1)
        with open(out1) as f1, open('./data/expected/test_scallop1.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove(out1)

    def test_scallop2(self):
        in1 = './data/input/R1_single.fastq'
        front_trim = 3
        back_trim = 6
        out1 = './data/output/trimmed.R1_single.fastq'
        scallop.scallop(in1, front_trim, back_trim, out1)
        with open(out1) as f1, open('./data/expected/test_scallop2.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove(out1)


if __name__ == "__main__":
    unittest.main()
