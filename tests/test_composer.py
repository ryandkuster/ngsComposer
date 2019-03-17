import unittest
import sys
import os
sys.path.append("..")
import composer


class TestComposer(unittest.TestCase):

    def test_initialize(self):
        result = composer.initialize(__file__)
        self.assertEqual(result, os.path.abspath(__file__))

#    def test_index_reader(self):
#        result = index_reader(bcs_index)


if __name__ == "__main__":
    unittest.main()

