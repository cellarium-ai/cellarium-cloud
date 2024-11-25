from unittest.mock import patch

import pandas as pd
import pytest

from casp.scripts.benchmarking.calculate_metrics import (
    calculate_f1,
    calculate_metrics_for_cas_output,
    calculate_metrics_for_query_cell,
    calculate_precision,
    calculate_tps_and_fps,
    split_into_batches,
)


@pytest.fixture
def mock_co_resource():
    return {
        "CL_0000197": {
            "all_ancestors": {"CL_0000000", "CL_0000003", "CL_0000197"},
            "all_descendants": {"CL_0000197"},
            "hop_0": {
                "nodes": {"CL_0000197"},
                "all_ancestors": {"CL_0000000", "CL_0000003", "CL_0000197"},
                "all_descendants": {"CL_0000197"},
            },
            "hop_1": {
                "nodes": {"CL_0000003", "CL_0000197"},
                "all_ancestors": {"CL_0000000", "CL_0000003", "CL_0000197"},
                "all_descendants": {"CL_0000003", "CL_0000197"},
            },
        },
        "CL_0000003": {
            "all_ancestors": {"CL_0000000", "CL_0000003"},
            "all_descendants": {"CL_0000003", "CL_0000197"},
            "hop_0": {
                "nodes": {"CL_0000003"},
                "all_ancestors": {"CL_0000000", "CL_0000003"},
                "all_descendants": {"CL_0000003", "CL_0000197"},
            },
            "hop_1": {
                "nodes": {"CL_0000000", "CL_0000003", "CL_0000197"},
                "all_ancestors": {"CL_0000000", "CL_0000003", "CL_0000197"},
                "all_descendants": {"CL_0000000", "CL_0000003", "CL_0000197"},
            },
        },
        "CL_0000000": {
            "all_ancestors": {"CL_0000000"},
            "all_descendants": {"CL_0000000", "CL_0000003", "CL_0000197"},
            "hop_0": {
                "nodes": {"CL_0000000"},
                "all_ancestors": {"CL_0000000"},
                "all_descendants": {"CL_0000000", "CL_0000003", "CL_0000197"},
            },
            "hop_1": {
                "nodes": {"CL_0000000", "CL_0000003"},
                "all_ancestors": {"CL_0000000", "CL_0000003"},
                "all_descendants": {"CL_0000000", "CL_0000003", "CL_0000197"},
            },
        },
    }


def test_calculate_precision():
    assert calculate_precision(tp=10, fp=5) == pytest.approx(expected=0.66666667, rel=1e-6)
    assert calculate_precision(tp=10, fp=0) == pytest.approx(expected=1.0, rel=1e-6)
    assert calculate_precision(tp=0, fp=10) == pytest.approx(expected=0.0, rel=1e-6)
    assert calculate_precision(tp=0, fp=0) == pytest.approx(expected=0.0, rel=1e-6)


def test_calculate_f1():
    assert calculate_f1(precision=0.66666667, recall=0.5) == pytest.approx(expected=0.57142857, rel=1e-6)
    assert calculate_f1(precision=1.0, recall=1.0) == pytest.approx(expected=1.0, rel=1e-6)
    assert calculate_f1(precision=0.0, recall=1.0) == pytest.approx(expected=0.0, rel=1e-6)
    assert calculate_f1(precision=1.0, recall=0.0) == pytest.approx(expected=0.0, rel=1e-6)
    assert calculate_f1(precision=0.0, recall=0.0) == pytest.approx(expected=0.0, rel=1e-6)


def test_calculate_tps_and_fps(mock_co_resource):
    query_cell_obj = {
        "matches": [
            {"cell_type_ontology_term_id": "CL_0000000", "score": 0.9},
            {"cell_type_ontology_term_id": "CL_0000003", "score": 0.7},
            {"cell_type_ontology_term_id": "CL_0000197", "score": 0.4},
        ]
    }
    true_positives, false_positives = calculate_tps_and_fps(
        query_cell_obj=query_cell_obj,
        ground_truth_cl_name="CL_0000000",
        num_hops=1,
        co_resource=mock_co_resource,
    )
    assert true_positives == pytest.approx(expected=[0.9, 0.9], rel=1e-6)
    assert false_positives == pytest.approx(expected=[0.0, 0.0], rel=1e-6)


def test_calculate_metrics_for_query_cell(mock_co_resource):
    query_cell_obj = {
        "query_cell_id": "cell_1",
        "matches": [
            {"cell_type_ontology_term_id": "CL_0000000", "score": 0.9},
            {"cell_type_ontology_term_id": "CL_0000003", "score": 0.7},
            {"cell_type_ontology_term_id": "CL_0000197", "score": 0.4},
        ],
    }
    metrics = calculate_metrics_for_query_cell(
        query_cell_obj=query_cell_obj,
        ground_truth_cl_name="CL_0000000",
        co_resource=mock_co_resource,
        num_hops=1,
    )
    expected_metrics = {
        "query_cell_id": "cell_1",
        "detail": "",
        "hop_0_sensitivity": 0.9,
        "hop_0_specificity": 1.0,
        "hop_0_f1_score": 0.9473684210526316,
        "hop_1_sensitivity": 0.9,
        "hop_1_specificity": 1.0,
        "hop_1_f1_score": 0.9473684210526316,
    }
    assert metrics == pytest.approx(expected=expected_metrics, rel=1e-6)


@patch("casp.scripts.benchmarking.calculate_metrics.calculate_metrics_for_query_cell")
def test_calculate_metrics_for_cas_output(mock_calculate_metrics_for_query_cell, mock_co_resource):
    mock_calculate_metrics_for_query_cell.side_effect = [
        {"query_cell_id": "cell_1", "hop_0_f1_score": 0.9, "hop_1_f1_score": 0.8},
        {"query_cell_id": "cell_2", "hop_0_f1_score": 0.8, "hop_1_f1_score": 0.7},
    ]
    ground_truth_cl_names = ["CL_0000000", "CL_0000197"]
    cas_result = [{"query_cell_id": "cell_1"}, {"query_cell_id": "cell_2"}]
    df_result = calculate_metrics_for_cas_output(
        ground_truth_cl_names=ground_truth_cl_names,
        cas_result=cas_result,
        co_resource=mock_co_resource,
        num_hops=1,
    )
    expected_df = pd.DataFrame(
        [
            {"query_cell_id": "cell_1", "hop_0_f1_score": 0.9, "hop_1_f1_score": 0.8},
            {"query_cell_id": "cell_2", "hop_0_f1_score": 0.8, "hop_1_f1_score": 0.7},
        ]
    ).set_index("query_cell_id")
    pd.testing.assert_frame_equal(df_result, expected_df)


def test_split_into_batches():
    data_list = list(range(10))
    batch_size = 3
    batches = split_into_batches(data_list=data_list, batch_size=batch_size)
    expected_batches = [
        [0, 1, 2],
        [3, 4, 5],
        [6, 7, 8],
        [9],
    ]
    assert batches == expected_batches
