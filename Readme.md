# Aptamer Generator and SELEX Simulator

A Python package for generating and evaluating candidate DNA aptamer sequences, including an in-silico SELEX simulation.

## Installation

You can install directly from GitHub:

```bash
pip install git+https://github.com/staru09/aptamer-generator.git
```

Alternatively, if you have cloned the repository and are working on the source code, you can install in editable mode:

```bash
pip install -e .
```

## Features

1. **Aptamer Generation**
   - Random sequence generation with GC content constraints
   - Seed-based sequence variation
   - Structure and folding potential evaluation

2. **In-silico SELEX**
   - Multi-round selection simulation
   - Target-based binding affinity
   - PCR amplification with mutations
   - Population diversity tracking

## Usage

### Basic Aptamer Generation

```python
from aptamer_generator.generator import generate_aptamers, evaluate_aptamer

# Generate aptamers
aptamers = generate_aptamers(
    length=30,  # Length of each aptamer
    count=10,   # Number of aptamers to generate
    gc_content_range=(0.4, 0.6)  # Desired GC content range
)

# Evaluate an aptamer
target_properties = {
    "target_length": 30,
    "target_gc_content": 0.5
}
score = evaluate_aptamer(aptamers[0]["sequence"], target_properties)
```

### SELEX Simulation

```python
from aptamer_generator.selex import InSilicoSELEX

# Initialize SELEX with target sequence
selex = InSilicoSELEX(target="GGTTGGTGTGGTTGG")  # Example: Thrombin binding aptamer

# Run SELEX simulation
final_library = selex.run(
    num_rounds=8,
    library_size=1000,
    selection_pressure=0.1,
    mutation_rate=0.05,
    gc_range=(0.4, 0.6),
    seq_length=30
)

# Plot evolution progress
selex.plot_progress()
```

## Examples

The package includes several example scripts:

1. Basic Aptamer Generation:
```bash
python examples/examples.py
```

2. SELEX Simulation:
```bash
python examples/examples_selex.py
```

Example outputs are saved as JSON files in the examples directory.

## Running Tests

```bash
# Run all tests
pytest

# Run specific test files
pytest tests/test_generator.py
pytest tests/test_selex.py
```
