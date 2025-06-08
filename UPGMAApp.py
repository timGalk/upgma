#!/usr/bin/env python3
"""
UPGMA Phylogenetic Tree Builder
A comprehensive tool for constructing phylogenetic trees using the UPGMA method.
Supports both sequence-based analysis and pre-computed distance matrices.

Author: Bioinformatics Tool
Version: 1.0
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import numpy as np
import pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import UPGMATreeConstructor, DistanceMatrix
from Bio.Phylo import draw_ascii, write
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import networkx as nx
from io import StringIO
import os
import re


class UPGMATreeBuilder:
    """Main class for UPGMA phylogenetic tree construction."""
    
    def __init__(self):
        """Initialize the UPGMA Tree Builder."""
        self.sequences = {}
        self.distance_matrix = None
        self.tree = None
        self.sequence_names = []
        
        # Setup GUI
        self.setup_gui()
    
    def setup_gui(self):
        """Create the main GUI interface."""
        self.root = tk.Tk()
        self.root.title("UPGMA Phylogenetic Tree Builder")
        self.root.geometry("1200x800")
        self.root.configure(bg='#f0f0f0')
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Tab 1: Input Data
        self.input_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.input_frame, text="Input Data")
        self.setup_input_tab()
        
        # Tab 2: Results
        self.results_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.results_frame, text="Results")
        self.setup_results_tab()
        
        # Tab 3: Tree Visualization
        self.viz_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.viz_frame, text="Tree Visualization")
        self.setup_visualization_tab()
        
        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Ready")
        self.status_bar = ttk.Label(self.root, textvariable=self.status_var, 
                                   relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
    
    def setup_input_tab(self):
        """Setup the input data tab."""
        # Input method selection
        method_frame = ttk.LabelFrame(self.input_frame, text="Input Method")
        method_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.input_method = tk.StringVar(value="sequences")
        ttk.Radiobutton(method_frame, text="Sequences (FASTA)", 
                       variable=self.input_method, value="sequences",
                       command=self.toggle_input_method).pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(method_frame, text="Distance Matrix", 
                       variable=self.input_method, value="matrix",
                       command=self.toggle_input_method).pack(side=tk.LEFT, padx=10)
        
        # Sequence input frame
        self.seq_frame = ttk.LabelFrame(self.input_frame, text="Sequence Input")
        self.seq_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # File input
        file_frame = ttk.Frame(self.seq_frame)
        file_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Button(file_frame, text="Load FASTA File", 
                  command=self.load_fasta_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(file_frame, text="Load Sample Data", 
                  command=self.load_sample_sequences).pack(side=tk.LEFT, padx=5)
        
        # Text input
        text_frame = ttk.Frame(self.seq_frame)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        ttk.Label(text_frame, text="Or paste sequences in FASTA format:").pack(anchor=tk.W)
        self.seq_text = scrolledtext.ScrolledText(text_frame, height=15, width=80)
        self.seq_text.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Matrix input frame
        self.matrix_frame = ttk.LabelFrame(self.input_frame, text="Distance Matrix Input")
        self.matrix_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        matrix_file_frame = ttk.Frame(self.matrix_frame)
        matrix_file_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Button(matrix_file_frame, text="Load Matrix File", 
                  command=self.load_matrix_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(matrix_file_frame, text="Load Sample Matrix", 
                  command=self.load_sample_matrix).pack(side=tk.LEFT, padx=5)
        
        ttk.Label(self.matrix_frame, text="Or paste distance matrix (CSV format):").pack(anchor=tk.W, padx=5)
        self.matrix_text = scrolledtext.ScrolledText(self.matrix_frame, height=15, width=80)
        self.matrix_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Control buttons
        button_frame = ttk.Frame(self.input_frame)
        button_frame.pack(fill=tk.X, padx=10, pady=10)
        
        ttk.Button(button_frame, text="Build Tree", 
                  command=self.build_tree).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear All", 
                  command=self.clear_all).pack(side=tk.LEFT, padx=5)
        
        # Initially hide matrix frame
        self.matrix_frame.pack_forget()
    
    def setup_results_tab(self):
        """Setup the results tab."""
        # Distance matrix display
        matrix_display_frame = ttk.LabelFrame(self.results_frame, text="Distance Matrix")
        matrix_display_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.matrix_display = scrolledtext.ScrolledText(matrix_display_frame, height=15)
        self.matrix_display.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Tree display
        tree_display_frame = ttk.LabelFrame(self.results_frame, text="UPGMA Tree (ASCII)")
        tree_display_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.tree_display = scrolledtext.ScrolledText(tree_display_frame, height=15)
        self.tree_display.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Save buttons
        save_frame = ttk.Frame(self.results_frame)
        save_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Button(save_frame, text="Save Tree (Newick)", 
                  command=self.save_tree_newick).pack(side=tk.LEFT, padx=5)
        ttk.Button(save_frame, text="Save Results (Text)", 
                  command=self.save_results_text).pack(side=tk.LEFT, padx=5)
    
    def setup_visualization_tab(self):
        """Setup the tree visualization tab."""
        # Create matplotlib figure
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, self.viz_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Visualization controls
        viz_controls = ttk.Frame(self.viz_frame)
        viz_controls.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Button(viz_controls, text="Refresh Visualization", 
                  command=self.update_visualization).pack(side=tk.LEFT, padx=5)
        ttk.Button(viz_controls, text="Save Plot", 
                  command=self.save_plot).pack(side=tk.LEFT, padx=5)
    
    def toggle_input_method(self):
        """Toggle between sequence and matrix input methods."""
        if self.input_method.get() == "sequences":
            self.seq_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
            self.matrix_frame.pack_forget()
        else:
            self.matrix_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
            self.seq_frame.pack_forget()
    
    def load_fasta_file(self):
        """Load sequences from a FASTA file."""
        file_path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fas *.fa"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                sequences = []
                for record in SeqIO.parse(file_path, "fasta"):
                    sequences.append(f">{record.id}\n{str(record.seq)}")
                
                self.seq_text.delete(1.0, tk.END)
                self.seq_text.insert(1.0, "\n".join(sequences))
                self.status_var.set(f"Loaded {len(sequences)} sequences from {os.path.basename(file_path)}")
                
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load FASTA file: {str(e)}")
    
    def load_sample_sequences(self):
        """Load sample sequences for demonstration."""
        sample_sequences = """
>Human
ATGCGATCGATCGATCG
>Chimp
ATGCGATCGATCGATCG
>Mouse
ATGCGATCGTTCGATCG
>Rat
ATGCGATCGTTCGATCG
>Dog
ATGCGATCAATCGATCG
        """.strip()
        
        self.seq_text.delete(1.0, tk.END)
        self.seq_text.insert(1.0, sample_sequences)
        self.status_var.set("Loaded sample sequences")
    
    def load_matrix_file(self):
        """Load distance matrix from a file."""
        file_path = filedialog.askopenfilename(
            title="Select matrix file",
            filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                with open(file_path, 'r') as f:
                    content = f.read()
                
                self.matrix_text.delete(1.0, tk.END)
                self.matrix_text.insert(1.0, content)
                self.status_var.set(f"Loaded matrix from {os.path.basename(file_path)}")
                
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load matrix file: {str(e)}")
    
    def load_sample_matrix(self):
        """Load sample distance matrix for demonstration."""
        sample_matrix = """
,Human,Chimp,Mouse,Rat,Dog
Human,0,0.1,0.7,0.8,0.6
Chimp,0.1,0,0.8,0.9,0.7
Mouse,0.7,0.8,0,0.2,0.9
Rat,0.8,0.9,0.2,0,1.0
Dog,0.6,0.7,0.9,1.0,0
        """.strip()
        
        self.matrix_text.delete(1.0, tk.END)
        self.matrix_text.insert(1.0, sample_matrix)
        self.status_var.set("Loaded sample distance matrix")
    
    def parse_sequences(self, fasta_text):
        """Parse FASTA sequences from text."""
        sequences = {}
        current_seq = ""
        current_id = ""
        
        for line in fasta_text.strip().split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current_id and current_seq:
                    sequences[current_id] = current_seq
                current_id = line[1:].strip()
                current_seq = ""
            else:
                current_seq += line
        
        if current_id and current_seq:
            sequences[current_id] = current_seq
        
        return sequences
    
    def calculate_pairwise_distances(self, sequences):
        """Calculate pairwise distances between sequences using global alignment."""
        seq_names = list(sequences.keys())
        n = len(seq_names)
        distances = np.zeros((n, n))
        
        # Initialize pairwise aligner
        aligner = PairwiseAligner()
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        
        self.status_var.set("Calculating pairwise distances...")
        
        for i in range(n):
            for j in range(i + 1, n):
                seq1 = Seq(sequences[seq_names[i]])
                seq2 = Seq(sequences[seq_names[j]])
                
                # Perform pairwise alignment
                alignments = aligner.align(seq1, seq2)
                best_alignment = alignments[0]
                
                # Calculate distance as 1 - (similarity score / max possible score)
                max_len = max(len(seq1), len(seq2))
                max_score = max_len * aligner.match_score
                similarity = best_alignment.score / max_score if max_score > 0 else 0
                distance = max(0, 1 - similarity)
                
                distances[i][j] = distance
                distances[j][i] = distance
        
        return distances, seq_names
    
    def parse_distance_matrix(self, matrix_text):
        """Parse distance matrix from CSV text."""
        try:
            # Read CSV data
            lines = matrix_text.strip().split('\n')
            if not lines:
                raise ValueError("Empty matrix")
            
            # Parse header
            header = lines[0].split(',')
            seq_names = [name.strip() for name in header[1:]]  # Skip first empty column
            
            # Parse data
            n = len(seq_names)
            distances = np.zeros((n, n))
            
            for i, line in enumerate(lines[1:]):
                if i >= n:
                    break
                values = line.split(',')
                for j, val in enumerate(values[1:]):  # Skip row name
                    if j >= n:
                        break
                    distances[i][j] = float(val.strip())
            
            return distances, seq_names
            
        except Exception as e:
            raise ValueError(f"Failed to parse distance matrix: {str(e)}")
    
    def validate_distance_matrix(self, distances, seq_names):
        """Validate the distance matrix."""
        n = len(seq_names)
        
        # Check dimensions
        if distances.shape != (n, n):
            raise ValueError("Distance matrix dimensions don't match sequence names")
        
        # Check symmetry
        if not np.allclose(distances, distances.T):
            raise ValueError("Distance matrix is not symmetric")
        
        # Check diagonal zeros
        if not np.allclose(np.diag(distances), 0):
            raise ValueError("Distance matrix diagonal should be zero")
        
        # Check non-negative values
        if np.any(distances < 0):
            raise ValueError("Distance matrix contains negative values")
        
        return True
    
    def build_tree(self):
        """Build UPGMA tree from input data."""
        try:
            if self.input_method.get() == "sequences":
                # Parse sequences
                fasta_text = self.seq_text.get(1.0, tk.END).strip()
                if not fasta_text:
                    messagebox.showerror("Error", "Please provide sequences")
                    return
                
                sequences = self.parse_sequences(fasta_text)
                if len(sequences) < 2:
                    messagebox.showerror("Error", "At least 2 sequences are required")
                    return
                
                # Calculate distances
                distances, seq_names = self.calculate_pairwise_distances(sequences)
                
                self.sequences = sequences
                self.sequence_names = seq_names
                
            else:
                # Parse distance matrix
                matrix_text = self.matrix_text.get(1.0, tk.END).strip()
                if not matrix_text:
                    messagebox.showerror("Error", "Please provide distance matrix")
                    return
                
                distances, seq_names = self.parse_distance_matrix(matrix_text)
                self.sequence_names = seq_names
            
            # Validate matrix
            self.validate_distance_matrix(distances, seq_names)
            self.distance_matrix = distances
            
            # Build UPGMA tree using BioPython
            self.status_var.set("Building UPGMA tree...")
            
            # Convert to BioPython DistanceMatrix format
            matrix_data = []
            for i in range(len(distances)):
                row = []
                for j in range(i + 1):
                    row.append(distances[i][j])
                matrix_data.append(row)
            
            dm = DistanceMatrix(seq_names, matrix_data)
            constructor = UPGMATreeConstructor()
            self.tree = constructor.upgma(dm)
            
            # Display results
            self.display_results()
            self.update_visualization()
            
            self.status_var.set("UPGMA tree built successfully")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to build tree: {str(e)}")
            self.status_var.set("Error building tree")
    
    def display_results(self):
        """Display results in the results tab."""
        # Display distance matrix
        self.matrix_display.delete(1.0, tk.END)
        if self.distance_matrix is not None:
            # Create formatted distance matrix display
            df = pd.DataFrame(self.distance_matrix, 
                            index=self.sequence_names, 
                            columns=self.sequence_names)
            matrix_str = df.round(4).to_string()
            self.matrix_display.insert(1.0, "Distance Matrix:\n\n" + matrix_str)
        
        # Display tree
        self.tree_display.delete(1.0, tk.END)
        if self.tree:
            # ASCII representation
            tree_str = StringIO()
            draw_ascii(self.tree, file=tree_str)
            ascii_tree = tree_str.getvalue()
            
            self.tree_display.insert(1.0, "UPGMA Tree (ASCII representation):\n\n" + ascii_tree)
    
    def update_visualization(self):
        """Update the tree visualization."""
        if not self.tree:
            return
        
        self.ax.clear()
        
        try:
            # Create networkx graph from tree
            G = nx.Graph()
            pos = {}
            labels = {}
            
            # Add nodes and edges
            self._add_tree_to_graph(self.tree.root, G, pos, labels, 0, 0, 1)
            
            # Draw the tree
            nx.draw(G, pos, ax=self.ax, with_labels=True, labels=labels,
                   node_color='lightblue', node_size=500, font_size=8,
                   font_weight='bold', edge_color='gray')
            
            self.ax.set_title("UPGMA Phylogenetic Tree", fontsize=14, fontweight='bold')
            self.ax.set_aspect('equal')
            
        except Exception as e:
            self.ax.text(0.5, 0.5, f"Visualization error: {str(e)}", 
                        transform=self.ax.transAxes, ha='center', va='center')
        
        self.canvas.draw()
    
    def _add_tree_to_graph(self, node, G, pos, labels, x, y, width):
        """Recursively add tree nodes to networkx graph."""
        node_id = id(node)
        G.add_node(node_id)
        pos[node_id] = (x, y)
        
        if node.is_terminal():
            labels[node_id] = node.name
        else:
            labels[node_id] = ""
        
        if hasattr(node, 'clades') and node.clades:
            child_width = width / len(node.clades)
            start_x = x - width/2 + child_width/2
            
            for i, child in enumerate(node.clades):
                child_x = start_x + i * child_width
                child_y = y - 0.2
                
                child_id = self._add_tree_to_graph(child, G, pos, labels, 
                                                 child_x, child_y, child_width)
                G.add_edge(node_id, child_id)
        
        return node_id
    
    def save_tree_newick(self):
        """Save tree in Newick format."""
        if not self.tree:
            messagebox.showerror("Error", "No tree to save")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Save tree",
            defaultextension=".newick",
            filetypes=[("Newick files", "*.newick"), ("Text files", "*.txt"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                write(self.tree, file_path, "newick")
                self.status_var.set(f"Tree saved to {os.path.basename(file_path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save tree: {str(e)}")
    
    def save_results_text(self):
        """Save complete results to text file."""
        if not self.tree:
            messagebox.showerror("Error", "No results to save")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Save results",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    f.write("UPGMA Phylogenetic Tree Analysis Results\n")
                    f.write("=" * 50 + "\n\n")
                    
                    # Input sequences
                    if self.sequences:
                        f.write("Input Sequences:\n")
                        for name, seq in self.sequences.items():
                            f.write(f">{name}\n{seq}\n")
                        f.write("\n")
                    
                    # Distance matrix
                    if self.distance_matrix is not None:
                        f.write("Distance Matrix:\n")
                        df = pd.DataFrame(self.distance_matrix,
                                        index=self.sequence_names,
                                        columns=self.sequence_names)
                        f.write(df.round(4).to_string())
                        f.write("\n\n")
                    
                    # Tree in Newick format
                    f.write("UPGMA Tree (Newick format):\n")
                    tree_str = StringIO()
                    write(self.tree, tree_str, "newick")
                    f.write(tree_str.getvalue())
                    f.write("\n")
                    
                    # ASCII tree
                    f.write("UPGMA Tree (ASCII representation):\n")
                    tree_str = StringIO()
                    draw_ascii(self.tree, file=tree_str)
                    f.write(tree_str.getvalue())
                
                self.status_var.set(f"Results saved to {os.path.basename(file_path)}")
                
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save results: {str(e)}")
    
    def save_plot(self):
        """Save the tree visualization plot."""
        if not self.tree:
            messagebox.showerror("Error", "No tree to save")
            return
        
        file_path = filedialog.asksaveasfilename(
            title="Save plot",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf"), 
                      ("SVG files", "*.svg"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                self.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                self.status_var.set(f"Plot saved to {os.path.basename(file_path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save plot: {str(e)}")
    
    def clear_all(self):
        """Clear all data and results."""
        self.sequences = {}
        self.distance_matrix = None
        self.tree = None
        self.sequence_names = []
        
        self.seq_text.delete(1.0, tk.END)
        self.matrix_text.delete(1.0, tk.END)
        self.matrix_display.delete(1.0, tk.END)
        self.tree_display.delete(1.0, tk.END)
        
        self.ax.clear()
        self.canvas.draw()
        
        self.status_var.set("All data cleared")
    
    def run(self):
        """Start the GUI application."""
        self.root.mainloop()


def main():
    """Main function to run the UPGMA Tree Builder."""
    app = UPGMATreeBuilder()
    app.run()


if __name__ == "__main__":
    main()