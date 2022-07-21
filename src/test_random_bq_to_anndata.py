import unittest
import numpy

from random_bq_to_anndata import convert_matrix_data_to_coo_matrix_input_format

class TestRandomBqToAnnData(unittest.TestCase):
    cell_ids = ["AAA", "AAC", "AAG", "AAT"]
    feature_ids = ["ENS1", "ENS2", "ENS3", "ENS4"]

    def test_one_cell_matrix(self):
        #       ENS1
        # AAA   1

        cell_ids = ["AAA"]
        feature_ids = ["ENS1"]

        matrix = numpy.full([len(cell_ids), len(feature_ids)], None)
        matrix[0][0] = 1

        expected_rows = [0]
        expected_columns = [0]
        expected_data = [1]

        self.run_test_on_matrix(cell_ids, feature_ids, matrix, expected_rows, expected_columns, expected_data)

    def test_really_sparse_matrix(self):
        #       ENS1 ENS2 ENS3 ENS4
        # AAA   n/a  n/a  n/a  n/a
        # AAC   n/a  n/a  n/a  n/a
        # AAG   n/a  n/a  88   n/a
        # AAT   n/a  n/a  n/a  n/a

        matrix = numpy.full([len(TestRandomBqToAnnData.cell_ids), len(TestRandomBqToAnnData.feature_ids)], None)
        matrix[2][2] = 88

        expected_rows    = [2]
        expected_columns = [2]
        expected_data    = [88]

        self.run_test_on_matrix(TestRandomBqToAnnData.cell_ids, TestRandomBqToAnnData.feature_ids, matrix, expected_rows, expected_columns, expected_data)


    def test_another_sparse_matrix(self):
        #       ENS1 ENS2 ENS3 ENS4
        # AAA   n/a  n/a  n/a  n/a
        # AAC   45   n/a  n/a  n/a
        # AAG   n/a  n/a  88   n/a
        # AAT   n/a  n/a  n/a  n/a

        matrix = numpy.full([len(TestRandomBqToAnnData.cell_ids), len(TestRandomBqToAnnData.feature_ids)], None)
        matrix[1][0] = 45
        matrix[2][2] = 88

        expected_rows    = [1, 2]
        expected_columns = [0, 2]
        expected_data    = [45, 88]

        self.run_test_on_matrix(TestRandomBqToAnnData.cell_ids, TestRandomBqToAnnData.feature_ids, matrix, expected_rows, expected_columns, expected_data)

    def test_sparse_matrix(self):
        #       ENS1 ENS2 ENS3 ENS4
        # AAA   1    7     2   n/a
        # AAC   0    n/a   4   n/a
        # AAG   3    n/a   0   1
        # AAT   n/a  n/a   n/a 9

        matrix = numpy.full([len(TestRandomBqToAnnData.cell_ids), len(TestRandomBqToAnnData.feature_ids)], None)
        matrix[0][0] = 1
        matrix[0][1] = 7
        matrix[0][2] = 2

        matrix[1][0] = 0
        matrix[1][2] = 4

        matrix[2][0] = 3
        matrix[2][2] = 0
        matrix[2][3] = 1

        matrix[3][3] = 9

        expected_rows    = [0, 0, 0, 1, 1, 2, 2, 2, 3]
        expected_columns = [0, 1, 2, 0, 2, 0, 2, 3, 3]
        expected_data    = [1, 7, 2, 0, 4, 3, 0, 1, 9]

        self.run_test_on_matrix(TestRandomBqToAnnData.cell_ids, TestRandomBqToAnnData.feature_ids, matrix, expected_rows, expected_columns, expected_data)

    def run_test_on_matrix(self, cell_ids, feature_ids, matrix, expected_rows, expected_columns, expected_data):
        cell_ids_to_row_num = {}
        for i in range (0, len(cell_ids)):
            cell_ids_to_row_num[cell_ids[i]] = i

        feature_ids_to_col_num = {}
        for i in range (0, len(feature_ids)):
            feature_ids_to_col_num[feature_ids[i]] = i

        cas_matrix_rows = []
        for i in range (0, len(cell_ids)):
            for j in range (0, len(feature_ids)):
                if (matrix[i][j] is not None):
                    cas_matrix_row = {}
                    cas_matrix_row["cas_cell_index"] = cell_ids[i]
                    cas_matrix_row["cas_feature_index"] = feature_ids[j]
                    cas_matrix_row["count"] = matrix[i][j]
                    cas_matrix_rows.append(cas_matrix_row)

        (rows, columns, data) = convert_matrix_data_to_coo_matrix_input_format(cas_matrix_rows, cell_ids_to_row_num, feature_ids_to_col_num)

        self.assertEqual(rows, expected_rows)
        self.assertEqual(columns, expected_columns)
        self.assertEqual(data, expected_data)


def main():
    unittest.main()

if __name__ == "__main__":
    main()