# example_selex.py
import json
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from selex.selex import InSilicoSELEX

def main():
    # Initialize SELEX with target motif
    selex = InSilicoSELEX(target="GGTTGGTGTGGTTGG")  # Thrombin binding aptamer motif
    
    # Run simulation
    final_library = selex.run(
        num_rounds=8,
        library_size=1500,
        selection_pressure=0.12,
        mutation_rate=0.04,
        gc_range=(0.45, 0.55),
        seq_length=30
    )
    
    # Show top 5 aptamers
    print("\nFinal Top Aptamers:")
    scored_library = selex._score_library(final_library)
    top_aptamers = sorted(scored_library, key=lambda x: x['score'], reverse=True)[:5]
    
    for idx, apt in enumerate(top_aptamers):
        print(f"{idx+1}. {apt['sequence']} | Score: {apt['score']:.3f}")
        print(f"   GC: {apt['gc']:.2%} | Folding: {apt['folding']:.2f}")
    
    # Save results to JSON
    output_file = os.path.join(os.path.dirname(__file__), "selex_results.json")
    with open(output_file, "w") as f:
        json.dump({
            "top_aptamers": top_aptamers,
            "round_history": selex.round_history
        }, f, indent=2)
    
    print(f"\nResults saved to {output_file}")
    
    # Plot progression
    selex.plot_progress()

if __name__ == "__main__":
    main()