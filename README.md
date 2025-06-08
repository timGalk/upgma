# UPGMA Phylogenetic Tree Generator

A Python application for constructing phylogenetic trees using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm. The program provides a user-friendly graphical interface for sequence analysis and tree visualization.

## Features

- **Multiple Input Methods**:
  - DNA/Protein sequences in FASTA format
  - Pre-computed distance matrices
- **Tree Construction**:
  - UPGMA algorithm implementation using Biopython
  - Automatic sequence alignment and distance calculation
  - Interactive tree visualization
- **Export Options**:
  - Save tree in Newick format
  - Export tree as high-quality image (PNG, PDF, SVG)
  - Generate detailed analysis reports
  - Save complete analysis package

## Installation

1. Ensure you have Python 3.7+ installed
2. Clone this repository or download the source code
3. Create and activate a virtual environment (recommended):
   ```bash
   python -m venv .venv
   # On Windows:
   .venv\Scripts\activate
   # On Unix/MacOS:
   source .venv/bin/activate
   ```
4. Install required packages:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Start the application:
   ```bash
   python App.py
   ```

2. **Input Data**:
   - Choose input type (sequences or distance matrix)
   - Load data from file or enter directly
   - For sequences, use FASTA format
   - For distance matrix, use space-separated values

3. **Generate Tree**:
   - Click "Generate Tree" to construct the phylogenetic tree
   - The tree will be displayed in the main window
   - Progress and information will be shown in the information panel

4. **Export Results**:
   - Save tree file (Newick format)
   - Export tree image
   - Generate analysis report
   - Save complete analysis package

## Input Format Examples

### FASTA Format (Sequences)
```
>Sequence1
ATGCATGC
>Sequence2
ATGCATGT
>Sequence3
ATGCATGA
```

### Distance Matrix Format
```
Sequence1 Sequence2 Sequence3
0.0 0.2 0.4
0.2 0.0 0.3
0.4 0.3 0.0
```

## Output Files

1. **Tree File (.nwk)**:
   - Newick format tree representation
   - Compatible with most phylogenetic software

2. **Tree Image (.png/.pdf/.svg)**:
   - High-resolution visualization
   - Customizable format options

3. **Analysis Report (.txt)**:
   - Detailed information about the analysis
   - Input data summary
   - Tree statistics
   - Sequence information (if applicable)

4. **Complete Analysis Package**:
   - All output files in a single directory
   - Timestamped for easy organization
   - Includes tree file, image, and report

## Dependencies

- biopython >= 1.81
- numpy >= 1.21.0
- matplotlib >= 3.4.0
- tk >= 0.1.0

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. 