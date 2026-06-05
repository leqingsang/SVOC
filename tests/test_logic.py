import pytest
from svoc import calculate_score_category 

def test_oncogenic_threshold():
    assert calculate_score_category(10) == "Oncogenic"
    assert calculate_score_category(12) == "Oncogenic"

def test_vus_range():
    assert calculate_score_category(0) == "VUS"
    assert calculate_score_category(5) == "VUS"

def test_benign_threshold():
    assert calculate_score_category(-7) == "Benign"
