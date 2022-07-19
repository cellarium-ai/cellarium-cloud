import unittest
import numpy

from random_bq_to_anndata import generate_sparse_matrix

class TestRandomBqToAnndata(unittest.TestCase):
    def test_one_cell_sparse_matrix(self):
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

    # def test_really_sparse_matrix(self):
    #
    #     #       ENS1 ENS2 ENS3 ENS4
    #     # AAA   n/a  n/a  n/a  n/a
    #     # AAC   n/a  n/a  n/a  n/a
    #     # AAG   n/a  n/a  88   n/a
    #     # AAT   n/a  n/a  n/a  n/a
    #     #
    #     # we represent it as:
    #     # index_ptr = [0, 1]
    #     # indices   = [0]
    #     # data      = [88]
    #
    #     cell_ids = ["AAA", "AAC", "AAG", "AAT"]
    #     feature_ids = ["ENS1", "ENS2", "ENS3", "ENS4"]
    #
    #     matrix = numpy.full([len(cell_ids), len(feature_ids)], None)
    #     matrix[2][2] = 88
    #
    #     expected_index_ptr = [0, 1]
    #     expected_indices   = [0]
    #     expected_data      = [88]
    #
    #     self.run_test_on_matrix(cell_ids, feature_ids, matrix, expected_index_ptr, expected_indices, expected_data)


    def test_generate_sparse_matrix(self):

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

        expected_cell_types = []
        expected_feature_names = []
        for i in range (0, len(cell_ids)):
            cell_id = cell_ids[i]
            cell_type = "CT_" + cell_id
            expected_cell_types.append(cell_type)
            for j in range (0, len(feature_ids)):
                feature_name = "FT_" + feature_ids[j]
                if (len(expected_feature_names) != len(feature_ids)):
                    expected_feature_names.append(feature_name)
                if (matrix[i][j] is not None):
                    cas_row = {}
                    cas_row["original_cell_id"] = cell_id
                    cas_row["cell_type"] = cell_type
                    cas_row["original_feature_id"] = feature_ids[j]
                    cas_row["feature_name"] = feature_name
                    cas_row["count"] = matrix[i][j]
                    cas_rows.append(cas_row)

        (index_ptr, indices, data, cell_names, cell_types, feature_ids, feature_names) = generate_sparse_matrix(cas_rows)

        self.assertEqual(index_ptr, expected_index_ptr)
        self.assertEqual(indices, expected_indices)
        self.assertEqual(data, expected_data)
        self.assertEqual(cell_names, cell_ids)
        self.assertEqual(cell_types, expected_cell_types)
        self.assertEqual(feature_ids, feature_ids)
        self.assertEqual(feature_names, expected_feature_names)

def main():
    unittest.main()

if __name__ == "__main__":
    main()