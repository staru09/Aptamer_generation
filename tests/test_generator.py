import pytest
from generator import generate_aptamers, evaluate_aptamer

def test_generate_aptamers_count():
    """Test if the correct number of aptamers are generated."""
    aptamers = generate_aptamers(count=5)
    assert len(aptamers) == 5

def test_generate_aptamers_length():
    """Test if aptamers have the correct length."""
    target_length = 25
    aptamers = generate_aptamers(length=target_length, count=3)
    for aptamer in aptamers:
        assert len(aptamer["sequence"]) == target_length

def test_gc_content_in_range():
    """Test if generated aptamers have GC content in the specified range."""
    gc_min, gc_max = 0.4, 0.6
    aptamers = generate_aptamers(count=10, gc_content_range=(gc_min, gc_max))
    
    for aptamer in aptamers:
        gc_content = aptamer["gc_content"]
        assert gc_min <= gc_content <= gc_max

def test_seed_sequence_variation():
    """Test if variation from seed sequence works."""
    seed = "ATGCATGCATGCATGCATGCATGCATGCAT"
    aptamers = generate_aptamers(count=5, seed_sequence=seed)
    
    # At least one aptamer should be different from the seed
    different_found = False
    for aptamer in aptamers:
        if aptamer["sequence"] != seed:
            different_found = True
            break
    
    assert different_found

def test_evaluate_aptamer():
    """Test if aptamer evaluation returns a valid score."""
    sequence = "ATGCATGCATGCATGCATGC"
    target_properties = {
        "target_length": 20,
        "target_gc_content": 0.5
    }
    
    score = evaluate_aptamer(sequence, target_properties)
    assert 0.0 <= score <= 1.0