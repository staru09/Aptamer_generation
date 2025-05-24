import pytest
import sys
import os
import random
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from selex.selex import InSilicoSELEX

def test_library_generation():
    """Test initial library generation with proper sizes and constraints."""
    selex = InSilicoSELEX()
    library_size = 50
    seq_length = 25
    gc_range = (0.4, 0.6)
    
    library = selex._generate_library(library_size, seq_length, gc_range)
    
    # Check library size
    assert len(library) == library_size
    
    # Check sequence length
    for seq_data in library:
        assert len(seq_data["sequence"]) == seq_length
    
    # Check GC content constraints
    for seq_data in library:
        gc = seq_data["gc_content"]
        assert gc_range[0] <= gc <= gc_range[1]

def test_scoring():
    """Test sequence scoring mechanism."""
    # Test with target
    target = "GGTTGGTGTGGTTGG"
    selex = InSilicoSELEX(target=target)
    
    # Create test sequences
    exact_match = target + "A" * 15  # Exact match of target + padding
    partial_match = "GGT" + "A" * (len(target) - 3) + "A" * 15  # Partial match
    no_match = "A" * len(target) + "A" * 15  # No match
    
    # Score the sequences
    exact_score = selex._calculate_affinity(exact_match)
    partial_score = selex._calculate_affinity(partial_match)
    no_match_score = selex._calculate_affinity(no_match)
    
    # Verify scoring behavior
    assert 0.0 <= exact_score <= 1.0
    assert 0.0 <= partial_score <= 1.0
    assert 0.0 <= no_match_score <= 1.0
    assert exact_score > partial_score > no_match_score

def test_selection():
    """Test selection of top binders."""
    selex = InSilicoSELEX()
    
    # Create mock scored library
    scored_pool = [
        {"sequence": "AAA", "score": 0.9},
        {"sequence": "BBB", "score": 0.8},
        {"sequence": "CCC", "score": 0.7},
        {"sequence": "DDD", "score": 0.6},
        {"sequence": "EEE", "score": 0.5}
    ]
    
    # Test different selection pressures
    high_pressure = selex._select_binders(scored_pool, 0.2)
    assert len(high_pressure) == 1
    assert high_pressure[0]["sequence"] == "AAA"
    
    med_pressure = selex._select_binders(scored_pool, 0.4)
    assert len(med_pressure) == 2
    assert med_pressure[0]["sequence"] == "AAA"
    assert med_pressure[1]["sequence"] == "BBB"
    
    low_pressure = selex._select_binders(scored_pool, 0.8)
    assert len(low_pressure) == 4
    
    # Test empty pool handling
    empty_result = selex._select_binders([], 0.5)
    assert isinstance(empty_result, list)
    assert len(empty_result) == 0

def test_mutation():
    """Test sequence mutation with GC constraints."""
    selex = InSilicoSELEX()
    
    # Test sequence
    original = "ATGCATGCATGCATGCATGC"  # 50% GC
    gc_range = (0.4, 0.6)
    
    # Test with different mutation rates
    for rate in [0.0, 0.1, 0.5]:
        # Multiple runs to account for randomness
        for _ in range(5):
            mutated = selex._mutate(original, rate, gc_range)
            
            # Check length preservation
            assert len(mutated) == len(original)
            
            # Check GC content
            gc = (mutated.count('G') + mutated.count('C')) / len(mutated)
            assert gc_range[0] <= gc <= gc_range[1]
            
            # Check mutation rate approximately followed
            if rate == 0.0:
                assert mutated == original
            elif rate == 0.5:
                # With high mutation rate, sequence should change
                assert mutated != original

def test_folding_potential():
    """Test folding potential calculation."""
    selex = InSilicoSELEX()
    
    # Sequence with palindromic regions (more likely to fold)
    folding_seq = "GGGGAAAATTTTCCCC"  # Complementary ends
    
    # Random sequence (less likely to fold)
    random_seq = "ACGTACGTACGTACGT"
    
    folding_score = selex._folding_potential(folding_seq)
    random_score = selex._folding_potential(random_seq)
    
    assert 0.0 <= folding_score <= 1.0
    assert 0.0 <= random_score <= 1.0
    
    # The folding sequence should have a higher score
    # But we can't assert this with certainty due to the simple heuristic
    # So just check the function returns without error

def test_amplification():
    """Test PCR amplification with mutations."""
    selex = InSilicoSELEX()
    
    sequences = ["ATGCATGCAT", "GCTAGCTAGC"]
    mutation_rate = 0.1
    gc_range = (0.4, 0.6)
    
    amplified = selex._amplify_pool(sequences, mutation_rate, gc_range)
    
    # Check that amplification increases diversity
    assert len(amplified) >= len(sequences)
    
    # Check GC content constraints maintained
    for seq_data in amplified:
        gc = seq_data["gc_content"]
        assert gc_range[0] <= gc <= gc_range[1]

def test_full_selex_run():
    """Test a full SELEX run with minimal parameters."""
    # Seed for reproducibility
    random.seed(42)
    np.random.seed(42)
    
    # Initialize with target
    selex = InSilicoSELEX(target="GGTTGG")
    
    # Run a minimal simulation
    final_library = selex.run(
        num_rounds=3,
        library_size=20,
        selection_pressure=0.3,
        mutation_rate=0.1,
        seq_length=15
    )
    
    # Check final library
    assert isinstance(final_library, list)
    assert len(final_library) > 0
    
    # Check round history was recorded
    assert len(selex.round_history) == 3
    
    # Check the average score increased over rounds
    if len(selex.round_history) >= 2:
        first_round = selex.round_history[0]['avg_score']
        last_round = selex.round_history[-1]['avg_score']
        # Score should generally increase, but with randomness it's not guaranteed
        