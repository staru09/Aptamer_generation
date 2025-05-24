import random
import numpy as np
from Bio.Seq import Seq
from typing import List, Dict, Union, Tuple, Optional

def generate_aptamers(
    length: int = 30,
    count: int = 10, 
    gc_content_range: Tuple[float, float] = (0.4, 0.6),
    seed_sequence: Optional[str] = None,
    constraints: Optional[Dict] = None
) -> List[Dict[str, Union[str, float]]]:
    """
    Generate candidate aptamer sequences.
    
    Args:
        length: Length of aptamer sequences to generate
        count: Number of aptamer candidates to generate
        gc_content_range: Desired GC content range as (min, max)
        seed_sequence: Optional seed sequence to base variations on
        constraints: Optional constraints for sequence generation
        
    Returns:
        List of dictionaries containing aptamer sequences and properties
    """
    nucleotides = ['A', 'T', 'G', 'C']
    aptamers = []
    
    constraints = constraints or {}
    
    for _ in range(count):
        if seed_sequence:
            # Generate variations based on seed sequence
            aptamer_seq = _generate_variation(seed_sequence, mutation_rate=0.2)
        else:
            # Generate new random sequence with specified GC content
            while True:
                seq = ''.join(random.choices(nucleotides, k=length))
                gc_count = seq.count('G') + seq.count('C')
                gc_content = gc_count / length
                
                if gc_content_range[0] <= gc_content <= gc_content_range[1]:
                    aptamer_seq = seq
                    break
        
        # Calculate properties
        properties = {
            "sequence": aptamer_seq,
            "length": len(aptamer_seq),
            "gc_content": (aptamer_seq.count('G') + aptamer_seq.count('C')) / len(aptamer_seq),
            "folding_score": _estimate_folding_potential(aptamer_seq),
        }
        
        aptamers.append(properties)
    
    return aptamers

def _generate_variation(seed_sequence: str, mutation_rate: float = 0.2) -> str:
    """Generate a variation of the seed sequence with random mutations."""
    nucleotides = ['A', 'T', 'G', 'C']
    variant = list(seed_sequence)
    
    for i in range(len(variant)):
        if random.random() < mutation_rate:
            # Replace with a different nucleotide
            current = variant[i]
            options = [n for n in nucleotides if n != current]
            variant[i] = random.choice(options)
    
    return ''.join(variant)

def _estimate_folding_potential(sequence: str) -> float:
    """
    Estimate the folding potential of a sequence using a simple heuristic.
    Real implementation would use proper RNA folding algorithms.
    """
    # This is a simplified approximation, real implementation would use
    # libraries like ViennaRNA or Nupack for proper secondary structure prediction
    # Count potential base pairs
    pairs = {
        'A': 'T',
        'T': 'A',
        'G': 'C', 
        'C': 'G'
    }
    
    score = 0.0
    seq_len = len(sequence)
    
    # Simple sliding window to detect potential hairpins
    for i in range(seq_len - 3):
        for j in range(i + 3, min(i + 15, seq_len)):
            if pairs.get(sequence[i]) == sequence[j]:
                # Add score for potential base pairing
                score += 1.0 / (j - i)
    
    return score

def evaluate_aptamer(sequence: str, target_properties: Dict) -> float:
    """
    Evaluate an aptamer sequence against desired target properties.
    
    Args:
        sequence: The aptamer sequence to evaluate
        target_properties: Dictionary of desired properties
        
    Returns:
        Score indicating how well the aptamer matches desired properties (0-1)
    """
    scores = []
    
    # Check length
    if "target_length" in target_properties:
        target_length = target_properties["target_length"]
        length_score = 1.0 - abs(len(sequence) - target_length) / target_length
        scores.append(length_score)
    
    # Check GC content
    if "target_gc_content" in target_properties:
        target_gc = target_properties["target_gc_content"]
        actual_gc = (sequence.count('G') + sequence.count('C')) / len(sequence)
        gc_score = 1.0 - abs(actual_gc - target_gc)
        scores.append(gc_score)
    
    # Check secondary structure potential
    folding_score = min(_estimate_folding_potential(sequence), 1.0)  # Cap at 1.0
    scores.append(folding_score)
    
    # Average all scores
    return sum(scores) / len(scores) if scores else 0.0