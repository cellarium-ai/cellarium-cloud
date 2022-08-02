"""
Tests of core matrix conversion functionality in `random_bq_to_anndata` as part of BQ extract.
"""

import unittest

import numpy

from casp import random_bq_to_anndata as rand_bq


class Expected:
    """
    Encapsulates an expectation for this test class.
    """

    def __init__(self, rows, columns, data):
        self.rows = rows
        self.columns = columns
        self.data = data


class TestRandomBqToAnnData(unittest.TestCase):
    """
    Various tests of matrix conversion as part of BQ extract.
    """

    cell_ids = ["AAA", "AAC", "AAG", "AAT"]
    feature_ids = ["ENS1", "ENS2", "ENS3", "ENS4"]

    def test_one_cell_matrix(self):
        """
        Single cell matrix test:

              ENS1
        AAA   1

        """

        cell_ids = ["AAA"]
        feature_ids = ["ENS1"]

        matrix = numpy.full([len(cell_ids), len(feature_ids)], None)
        matrix[0][0] = 1

        expected = Expected([0], [0], [1])

        self.run_test_on_matrix(cell_ids, feature_ids, matrix, expected)

    def test_really_sparse_matrix(self):
        """
        Really sparse matrix test:

              ENS1 ENS2 ENS3 ENS4
        AAA   n/a  n/a  n/a  n/a
        AAC   n/a  n/a  n/a  n/a
        AAG   n/a  n/a  88   n/a
        AAT   n/a  n/a  n/a  n/a

        """

        matrix = numpy.full([len(TestRandomBqToAnnData.cell_ids), len(TestRandomBqToAnnData.feature_ids)], None)
        matrix[2][2] = 88

        expected = Expected([2], [2], [88])

        self.run_test_on_matrix(
            TestRandomBqToAnnData.cell_ids,
            TestRandomBqToAnnData.feature_ids,
            matrix,
            expected,
        )

    def test_another_sparse_matrix(self):
        """
        Another sparse matrix test:

              ENS1 ENS2 ENS3 ENS4
        AAA   n/a  n/a  n/a  n/a
        AAC   45   n/a  n/a  n/a
        AAG   n/a  n/a  88   n/a
        AAT   n/a  n/a  n/a  n/a

        """

        matrix = numpy.full([len(TestRandomBqToAnnData.cell_ids), len(TestRandomBqToAnnData.feature_ids)], None)
        matrix[1][0] = 45
        matrix[2][2] = 88

        expected = Expected(rows=[1, 2], columns=[0, 2], data=[45, 88])

        self.run_test_on_matrix(
            TestRandomBqToAnnData.cell_ids,
            TestRandomBqToAnnData.feature_ids,
            matrix,
            expected,
        )

    def test_sparse_matrix(self):
        """
        And another sparse matrix test:

              ENS1 ENS2 ENS3 ENS4
        AAA   1    7     2   n/a
        AAC   0    n/a   4   n/a
        AAG   3    n/a   0   1
        AAT   n/a  n/a   n/a 9

        """

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

        expected = Expected(
            rows=[0, 0, 0, 1, 1, 2, 2, 2, 3], columns=[0, 1, 2, 0, 2, 0, 2, 3, 3], data=[1, 7, 2, 0, 4, 3, 0, 1, 9]
        )

        self.run_test_on_matrix(
            TestRandomBqToAnnData.cell_ids,
            TestRandomBqToAnnData.feature_ids,
            matrix,
            expected,
        )

    def run_test_on_matrix(self, cell_ids, feature_ids, matrix, expected):
        """
        Worker method for running matrix tests.
        """
        cell_ids_to_row_num = {}
        for i, cell_id in enumerate(cell_ids):
            cell_ids_to_row_num[cell_id] = i

        feature_ids_to_col_num = {}
        for i, feature_id in enumerate(feature_ids):
            feature_ids_to_col_num[feature_id] = i

        cas_matrix_rows = []
        for i, cell_id in enumerate(cell_ids):
            for j, feature_id in enumerate(feature_ids):
                if matrix[i][j] is not None:
                    cas_matrix_row = {"cas_cell_index": cell_id, "cas_feature_index": feature_id, "count": matrix[i][j]}
                    cas_matrix_rows.append(cas_matrix_row)

        (rows, columns, data) = rand_bq.convert_matrix_data_to_coo_matrix_input_format(
            cas_matrix_rows, cell_ids_to_row_num, feature_ids_to_col_num
        )

        self.assertEqual(rows, expected.rows)
        self.assertEqual(columns, expected.columns)
        self.assertEqual(data, expected.data)


def main():
    """
    Standard entry point that forwards to unit test main.
    """
    unittest.main()


if __name__ == "__main__":
    main()
