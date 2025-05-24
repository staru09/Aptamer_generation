import json
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generator import generate_aptamers, evaluate_aptamer

def main():
    aptamers = generate_aptamers(
        length=30,
        count=5,
        gc_content_range=(0.45, 0.55)
    )
    
    print("Generated Aptamers:")
    for i, aptamer in enumerate(aptamers):
        print(f"Aptamer {i+1}:")
        print(f"  Sequence: {aptamer['sequence']}")
        print(f"  Length: {aptamer['length']}")
        print(f"  GC Content: {aptamer['gc_content']:.2f}")
        print(f"  Folding Score: {aptamer['folding_score']:.4f}")
    
    # Fix: Use os.path to create correct file path
    output_file = os.path.join(os.path.dirname(__file__), "example_data.json")
    with open(output_file, "w") as f:
        json.dump(aptamers, f, indent=2)
    
    # Evaluate the first aptamer
    first_aptamer = aptamers[0]["sequence"]
    target_properties = {
        "target_length": 30,
        "target_gc_content": 0.5
    }
    
    score = evaluate_aptamer(first_aptamer, target_properties)
    print(f"\nEvaluation score for first aptamer: {score:.4f}")

if __name__ == "__main__":
    main()