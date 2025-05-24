import random
from typing import List, Dict, Tuple, Optional
import numpy as np

class InSilicoSELEX:
    def __init__(self, target: Optional[str] = None):
        self.target = target
        self.round_history = []
        
    def run(
        self,
        num_rounds: int = 8,
        library_size: int = 1000,
        selection_pressure: float = 0.1,
        mutation_rate: float = 0.05,
        gc_range: Tuple[float, float] = (0.4, 0.6),
        seq_length: int = 30
    ) -> List[Dict]:
        """
        Execute full SELEX protocol
        
        Args:
            num_rounds: Number of selection rounds
            library_size: Initial pool size
            selection_pressure: Fraction of pool to keep each round
            mutation_rate: Probability of mutation during amplification
            gc_range: Acceptable GC content range
            seq_length: Oligonucleotide length
            
        Returns:
            Final enriched library with metadata
        """
        library = self._generate_library(library_size, seq_length, gc_range)
        
        for round_num in range(1, num_rounds + 1):
            library = self._selection_round(
                library,
                round_num,
                selection_pressure,
                mutation_rate,
                gc_range
            )
            
        return library

    
    def _selection_round(
        self,
        library: List[Dict],
        round_num: int,
        selection: float,
        mutation_rate: float,
        gc_range: Tuple[float, float]
    ) -> List[Dict]:
        """Execute one complete selection round"""
        # 1. Binding Selection
        scored_pool = self._score_library(library)
        selected = self._select_binders(scored_pool, selection)

        if not selected:
            print(f"Warning: No binders found in round {round_num}. Using top sequence.")
            if scored_pool:
                selected = scored_pool[0]
            else:
                new_seq = self._generate_library(1, seq_length=30, gc_range=gc_range)[0]
                selected = [{'sequence': new_seq['sequence'], 'score':0.1}]
        
        # 2. Amplification
        amplified = self._amplify_pool(
            [seq['sequence'] for seq in selected],
            mutation_rate,
            gc_range
        )
        
        # 3. Analysis and tracking
        round_data = {
            'round': round_num,
            'avg_score': np.mean([s['score'] for s in selected]),
            'top_seq': selected[0]['sequence'],
            'diversity': len(set([s['sequence'] for s in amplified]))
        }
        self.round_history.append(round_data)
        
        self._print_round_summary(round_data)
        return amplified

    def _score_library(self, library: List[Dict]) -> List[Dict]:
        """Calculate binding scores for all sequences"""
        return [{
            'sequence': seq['sequence'],
            'score': self._calculate_affinity(seq['sequence']),
            'gc': seq['gc_content'],
            'folding': seq['folding_score']
        } for seq in library]

    def _calculate_affinity(self, sequence: str) -> float:
        """Calculate binding affinity score (0-1 scale)"""
        score = 0.0
        
        # Target complementarity component
        if self.target:
            overlap = sum(1 for s, t in zip(sequence, self.target) if s == t)
            score += 0.6 * (overlap / len(self.target))
        
        # Structural component
        score += 0.25 * self._folding_potential(sequence)
        
        # Stability component
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        score += 0.15 * (1 - abs(gc_content - 0.5))
        
        return min(score, 1.0)  # Cap at 1.0

    def _select_binders(
        self,
        scored_pool: List[Dict],
        selection: float
    ) -> List[Dict]:
        """Select top performing sequences"""
        if not scored_pool:
            return []
        
        sorted_pool = sorted(scored_pool, key=lambda x: x['score'], reverse=True)
        num_to_select = max(1, int(len(sorted_pool) * selection))
        return sorted_pool[:num_to_select]

    def _amplify_pool(
        self,
        sequences: List[str],
        mutation_rate: float,
        gc_range: Tuple[float, float]
    ) -> List[Dict]:
        """PCR amplification with mutagenesis"""
        # Mutate sequences
        mutated = [self._mutate(seq, mutation_rate, gc_range) for seq in sequences]
        
        # Duplicate while maintaining original library size
        amplified = random.choices(mutated + sequences, k=len(sequences))
        
        return self._generate_library(len(amplified), gc_range=gc_range, pool=amplified)


    def _generate_library(
        self,
        size: int,
        length: int = 30,
        gc_range: Tuple[float, float] = (0.4, 0.6),
        pool: Optional[List[str]] = None
    ) -> List[Dict]:
        """Generate sequence library with GC constraints"""
        if pool:  # Amplification from existing pool
            valid_seqs = [seq for seq in pool if self._check_gc(seq, gc_range)]
            return [self._seq_metadata(s) for s in random.choices(valid_seqs, k=size)]
        
        # Generate new library
        library = []
        while len(library) < size:
            seq = ''.join(random.choices('ATGC', k=length))
            if self._check_gc(seq, gc_range):
                library.append(self._seq_metadata(seq))
        return library

    def _mutate(
        self,
        sequence: str,
        rate: float,
        gc_range: Tuple[float, float]
    ) -> str:
        """Introduce mutations while maintaining GC constraints"""
        nt = ['A', 'T', 'G', 'C']
        seq = list(sequence)
        for i in range(len(seq)):
            if random.random() < rate:
                original = seq[i]
                candidates = [n for n in nt if n != original]
                random.shuffle(candidates)
                
                # Find mutation that preserves GC constraints
                for candidate in candidates:
                    new_seq = seq[:i] + [candidate] + seq[i+1:]
                    if self._check_gc(new_seq, gc_range):
                        seq = new_seq
                        break
        return ''.join(seq)

    def _folding_potential(self, sequence: str) -> float:
        """Estimate secondary structure potential"""
        score = 0.0
        pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        for i in range(len(sequence) - 4):
            for j in range(i+4, min(i+20, len(sequence))):
                if pairs.get(sequence[i]) == sequence[j]:
                    score += 1/(j-i)  # Closer pairs get higher scores
        return score / 10  # Normalized heuristic

    def _check_gc(self, sequence: str, gc_range: Tuple[float, float]) -> bool:
        gc = (sequence.count('G') + sequence.count('C')) / len(sequence)
        return gc_range[0] <= gc <= gc_range[1]

    def _seq_metadata(self, sequence: str) -> Dict:
        return {
            'sequence': sequence,
            'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence),
            'folding_score': self._folding_potential(sequence)
        }

    def _print_round_summary(self, data: Dict):
        print(f"\nRound {data['round']} Results:")
        print(f"  Average Score: {data['avg_score']:.3f}")
        print(f"  Sequence Diversity: {data['diversity']} unique sequences")
        print(f"  Top Sequence: {data['top_seq']}")

    def plot_progress(self):
        """Generate evolution plot (requires matplotlib)"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Matplotlib required for plotting")
            return

        rounds = [r['round'] for r in self.round_history]
        scores = [r['avg_score'] for r in self.round_history]
        diversity = [r['diversity'] for r in self.round_history]

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        ax1.plot(rounds, scores, 'b-', marker='o', label='Binding Score')
        ax2.plot(rounds, diversity, 'r--', marker='s', label='Diversity')

        ax1.set_xlabel('SELEX Round')
        ax1.set_ylabel('Binding Score', color='b')
        ax2.set_ylabel('Sequence Diversity', color='r')
        plt.title('SELEX Progression')
        fig.legend(loc='upper right')
        plt.show()

