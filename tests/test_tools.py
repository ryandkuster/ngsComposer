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


# class TestAnemone(unittest.TestCase):


# class TestCrinoid(unittest.TestCase):

#     def test_crinoid_main(self):
#         pass

#     def test_crinoid_comp(self):
#         pass
#         crinoid.crinoid_comp(curr, all_qc, p64, in1)

#     def test_crinoid_open(self):
#         pass
#         crinoid.crinoid_open(in1, out1, out2, p64)

#     def test_crinoid(self):
#         pass
#         crinoid.crinoid(f, out1, out2, p64)

#     def test_bespoke_matrix(self):
#         pass
#         crinoid.bespoke_matrix(mx, max_len)

#     def test_matrix_print(self):
#         pass
#         crinoid.matrix_print(mx, outfile)

#     def test_visualizer(self):
#         pass
#         crinoid.visualizer(out1, out2)

#     def test_combine_matrix(self):
#         pass
#         crinoid.combine_matrix(in_ls, out)

#     def test_matrix_mash(self):
#         pass
#         crinoid.matrix_mash(in_ls, a1, a2)


# class TestKrill(unittest.TestCase):


# class TestPorifera(unittest.TestCase):


class TestRotifer(unittest.TestCase):

    # def test_rotifer_main(self):
    #     pass

    # def test_rotifer_comp(self):
    #     pass

    # def test_rotifer_open(self):
    #     pass
    #     rotifer.rotifer_open(R1_bases_ls, R2_bases_ls, in1, in2, pe_1, pe_2,
    #                          se_1, se_2, trim)
 
    # def test_rotifer_trim(self):
    #     pass
    #     rotifer.rotifer_trim(R1_bases_ls, R2_bases_ls, f1, f2, pe_o1, pe_o2,
    #                          se_o1, se_o2, trim):

    def test_rotifer1(self):
        '''
        test pe files with R1 motifs
        '''
        R1_bases_ls = ['AA', 'CC', 'GG', 'TT']
        R2_bases_ls = None
        with open('./data/input/R1_1000.fastq') as f1, \
            open('./data/input/R2_1000.fastq') as f2, \
            open('./data/output/pe.R1_single.fastq', 'w') as pe_o1, \
            open('./data/output/pe.R2_single.fastq', 'w') as pe_o2, \
            open('./data/output/se.R1_single.fastq', 'w') as se_o1, \
            open('./data/output/se.R2_single.fastq', 'w') as se_o2:
            rotifer.rotifer(R1_bases_ls, R2_bases_ls, f1, f2, pe_o1, pe_o2, se_o1, se_o2)
        with open('./data/output/pe.R1_single.fastq') as f1, \
            open('./data/expected/test_rotifer1.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove('./data/output/pe.R1_single.fastq')
        os.remove('./data/output/pe.R2_single.fastq')
        os.remove('./data/output/se.R1_single.fastq')
        os.remove('./data/output/se.R2_single.fastq')

    def test_rotifer2(self):
        '''
        test pe files with R1 and R2 motifs
        '''
        R1_bases_ls = ['AA', 'CC', 'GG', 'TT']
        R2_bases_ls = ['AA', 'CC', 'GG', 'TT']
        with open('./data/input/R1_1000.fastq') as f1, \
            open('./data/input/R2_1000.fastq') as f2, \
            open('./data/output/pe.R1_single.fastq', 'w') as pe_o1, \
            open('./data/output/pe.R2_single.fastq', 'w') as pe_o2, \
            open('./data/output/se.R1_single.fastq', 'w') as se_o1, \
            open('./data/output/se.R2_single.fastq', 'w') as se_o2:
            rotifer.rotifer(R1_bases_ls, R2_bases_ls, f1, f2, pe_o1, pe_o2, se_o1, se_o2)
        with open('./data/output/pe.R2_single.fastq') as f1, \
            open('./data/expected/test_rotifer2.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove('./data/output/pe.R1_single.fastq')
        os.remove('./data/output/pe.R2_single.fastq')
        os.remove('./data/output/se.R1_single.fastq')
        os.remove('./data/output/se.R2_single.fastq')

    # def test_rotifer_single_open(self):
    #     pass
    #     rotifer.single_open(R1_bases_ls, in1, se_1, trim)

    # def test_rotifer_single_trim(self):
    #     pass
    #     rotifer.single_trim(R1_bases_ls, f1, se_o1, trim)

    def test_rotifer_single(self):

        R1_bases_ls = ['GG', 'TT']
        with open('./data/input/R1_1000.fastq') as f1, \
            open('./data/output/se.R1_single.fastq', 'w') as se_o1:
            rotifer.rotifer_single(R1_bases_ls, f1, se_o1)
        with open('./data/output/se.R1_single.fastq') as f1, \
            open('./data/expected/test_rotifer_single.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove('./data/output/se.R1_single.fastq')

    def test_rotifer_test1(self):
        '''
        test that motif test returns False
        '''
        line = 'ACGT'
        bases_ls = ['AG', 'AC']
        result = rotifer.rotifer_test(line, bases_ls)
        self.assertEqual(result, False)

    def test_rotifer_test2(self):
        '''
        test that motif test returns True
        '''
        line = 'ACGT'
        bases_ls = ['GG', 'TC']
        result = rotifer.rotifer_test(line, bases_ls)
        self.assertEqual(result, True)


class TestScallop(unittest.TestCase):

    # def test_scallop_main(self):
    #     pass
    #     scallop.scallop_main()

    # def test_scallop_comp(self):
    #     pass
    #     scallop.scallop_comp(front_trim, end_trim, curr, in1)

    # def test_scallop_open(self):
    #     pass
    #     scallop_open(in1, front_trim, end_trim, out1)

    def test_scallop1(self):
        '''
        test front_trim behavior
        '''
        front_trim = 20
        end_trim = None
        with open('./data/input/R1_single.fastq') as f, \
            open('./data/output/trimmed.R1_single.fastq', 'w') as o:
            scallop.scallop(front_trim, end_trim, f, o)
        with open('./data/output/trimmed.R1_single.fastq') as f1, \
            open('./data/expected/test_scallop1.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove('./data/output/trimmed.R1_single.fastq')

    def test_scallop2(self):
        '''
        test front_trim and end_trim behavior
        '''
        front_trim = 3
        end_trim = 6
        with open('./data/input/R1_single.fastq') as f, \
            open('./data/output/trimmed.R1_single.fastq', 'w') as o:
            scallop.scallop(front_trim, end_trim, f, o)
        with open('./data/output/trimmed.R1_single.fastq') as f1, \
            open('./data/expected/test_scallop2.fastq') as f2:
            self.assertEqual(f1.readlines(), f2.readlines())
        os.remove('./data/output/trimmed.R1_single.fastq')

    # def test_scallop_end(self):
    #     pass
    #     scallop.scallop_end(curr, auto_trim, trim_mode, in1)

    # def scallop_stats(self):
    #     pass
    #     scallop.scallop_stats(target_index, pos):


if __name__ == "__main__":
    unittest.main()
