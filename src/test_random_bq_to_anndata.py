import unittest
import numpy

from random_bq_to_anndata import generate_sparse_matrix

class TestRandomBqToAnnData(unittest.TestCase):
    def test_one_cell_matrix(self):
        #       ENS1
        # AAA   1

        cell_ids = ["AAA"]
        feature_ids = ["ENS1"]

        matrix = numpy.full([len(cell_ids), len(feature_ids)], None)
        matrix[0][0] = 1

        expected_index_ptr = [0, 1]
        expected_indices   = [0]
        expected_data      = [1]

        self.run_test_on_matrix(cell_ids, feature_ids, matrix, expected_index_ptr, expected_indices, expected_data)

    def test_really_sparse_matrix(self):

        #       ENS1 ENS2 ENS3 ENS4
        # AAA   n/a  n/a  n/a  n/a
        # AAC   n/a  n/a  n/a  n/a
        # AAG   n/a  n/a  88   n/a
        # AAT   n/a  n/a  n/a  n/a
        #
        # we represent it as:
        # index_ptr = [0, 1]
        # indices   = [2]
        # data      = [88]

        cell_ids = ["AAA", "AAC", "AAG", "AAT"]
        feature_ids = ["ENS1", "ENS2", "ENS3", "ENS4"]

        matrix = numpy.full([len(cell_ids), len(feature_ids)], None)
        matrix[2][2] = 88

        # Note that the expected index_ptr, indices, and data contain the '0' data point that `generate_sparse_matrix`
        # will add in order to ensure that there is data for the maximum possible feature id / column
        expected_index_ptr = [0, 2]
        expected_indices   = [2, 3]
        expected_data      = [88, 0]

        self.run_test_on_matrix(cell_ids, feature_ids, matrix, expected_index_ptr, expected_indices, expected_data)


    def test_another_sparse_matrix(self):

        #       ENS1 ENS2 ENS3 ENS4
        # AAA   n/a  n/a  n/a  n/a
        # AAC   45   n/a  n/a  n/a
        # AAG   n/a  n/a  88   n/a
        # AAT   n/a  n/a  n/a  n/a
        #
        # we represent it as:
        # index_ptr = [0, 1, 2]
        # indices   = [0, 2, 3]
        # data      = [45, 88, 0]

        cell_ids = ["AAA", "AAC", "AAG", "AAT"]
        feature_ids = ["ENS1", "ENS2", "ENS3", "ENS4"]

        matrix = numpy.full([len(cell_ids), len(feature_ids)], None)
        matrix[1][0] = 45
        matrix[2][2] = 88

        # Note that the expected index_ptr, indices, and data contain the '0' data point that `generate_sparse_matrix`
        # will add in order to ensure that there is data for the maximum possible feature id / column
        expected_index_ptr = [0, 1, 3]
        expected_indices   = [0, 2, 3]
        expected_data      = [45, 88, 0]

        self.run_test_on_matrix(cell_ids, feature_ids, matrix, expected_index_ptr, expected_indices, expected_data)

    def test_sparse_matrix(self):

        #       ENS1 ENS2 ENS3 ENS4
        # AAA   1    7     2   n/a
        # AAC   0    n/a   4   n/a
        # AAG   3    n/a   0   1
        # AAT   n/a  n/a   n/a 9
        #
        # we represent it as:
        # index_ptr = [0, 3, 5, 8, 9]
        # indices   = [0, 1, 2, 0, 2, 0, 2, 3, 3]
        # data      = [1, 7, 2, 0, 4, 3, 0, 1, 9]

        cell_ids = ["AAA", "AAC", "AAG", "AAT"]
        feature_ids = ["ENS1", "ENS2", "ENS3", "ENS4"]

        matrix = numpy.full([len(cell_ids), len(feature_ids)], None)
        matrix[0][0] = 1
        matrix[0][1] = 7
        matrix[0][2] = 2

        matrix[1][0] = 0
        matrix[1][2] = 4

        matrix[2][0] = 3
        matrix[2][2] = 0
        matrix[2][3] = 1

        matrix[3][3] = 9

        expected_index_ptr = [0, 3, 5, 8, 9]
        expected_indices   = [0, 1, 2, 0, 2, 0, 2, 3, 3]
        expected_data      = [1, 7, 2, 0, 4, 3, 0, 1, 9]

        self.run_test_on_matrix(cell_ids, feature_ids, matrix, expected_index_ptr, expected_indices, expected_data)


    def run_test_on_matrix(self, cell_ids, feature_ids, matrix, expected_index_ptr, expected_indices, expected_data):
        cas_rows = []

        for i in range (0, len(cell_ids)):
            cell_id = cell_ids[i]
            for j in range (0, len(feature_ids)):
                if (matrix[i][j] is not None):
                    cas_row = {}
                    cas_row["cas_cell_index"] = cell_id
                    cas_row["cas_feature_index"] = j
                    cas_row["count"] = matrix[i][j]
                    cas_rows.append(cas_row)

        (index_ptr, indices, data) = generate_sparse_matrix(cas_rows, 0, len(feature_ids) - 1)

        self.assertEqual(index_ptr, expected_index_ptr)
        self.assertEqual(indices, expected_indices)
        self.assertEqual(data, expected_data)


def main():
    unittest.main()

if __name__ == "__main__":
    main()