import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from Bio import SeqIO, Phylo, pairwise2
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import datetime
from pathlib import Path
import sys
import tempfile

def get_resource_path(relative_path):
    """Get absolute path to resource, works for dev and for PyInstaller"""
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

# Set Biopython data path
if getattr(sys, 'frozen', False):
    # Running as compiled executable
    bio_data_path = get_resource_path('Bio')
    os.environ['BIO_DATA_PATH'] = bio_data_path

class PhyloApp:
    """
    A graphical application for constructing phylogenetic trees using the UPGMA algorithm.
    
    This class provides a user-friendly interface for:
    - Loading DNA/protein sequences or distance matrices
    - Computing phylogenetic trees using UPGMA
    - Visualizing and exporting results
    
    Attributes:
        root (tk.Tk): The main application window
        sequences (list): List of loaded sequences
        seq_ids (list): List of sequence identifiers
        distance_matrix (DistanceMatrix): Loaded distance matrix
        current_tree (Tree): The generated phylogenetic tree
        current_figure (Figure): The matplotlib figure for tree visualization
        option (StringVar): Input type selection (sequences/matrix)
    """
    
    def __init__(self, root):
        """
        Initialize the PhyloApp application.
        
        Args:
            root (tk.Tk): The main application window
        """
        self.root = root
        self.root.title("UPGMA Phylogenetic Tree Generator")
        self.root.geometry("1200x800")
        
        # Initialize variables
        self.sequences = []
        self.seq_ids = []
        self.distance_matrix = None
        self.current_tree = None
        self.current_figure = None
        self.option = tk.StringVar(value='sequences')
        
        # Create main interface
        self.create_widgets()
        
    def create_widgets(self):
        """
        Create and arrange all GUI widgets.
        
        This method sets up the main application layout including:
        - Control panel for input and operations
        - Display panel for tree visualization
        - Information panel for status updates
        """
        # Create main container with paned window
        main_paned = ttk.PanedWindow(self.root, orient='horizontal')
        main_paned.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Left panel for controls
        left_frame = ttk.Frame(main_paned)
        main_paned.add(left_frame, weight=1)
        
        # Right panel for tree display
        right_frame = ttk.Frame(main_paned)
        main_paned.add(right_frame, weight=3)
        
        self.create_control_panel(left_frame)
        self.create_display_panel(right_frame)
        
    def create_control_panel(self, parent):
        """
        Create the control panel with input options and operation buttons.
        
        Args:
            parent (ttk.Frame): The parent frame for the control panel
        """
        # Main control frame
        control_frame = ttk.LabelFrame(parent, text="Controls", padding=10)
        control_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Input options
        input_frame = ttk.LabelFrame(control_frame, text="Input Type", padding=5)
        input_frame.pack(fill='x', pady=(0, 10))
        
        ttk.Radiobutton(input_frame, text="DNA/Protein Sequences (FASTA)", 
                       variable=self.option, value='sequences',
                       command=self.toggle_input_type).pack(anchor='w')
        ttk.Radiobutton(input_frame, text="Distance Matrix", 
                       variable=self.option, value='matrix',
                       command=self.toggle_input_type).pack(anchor='w')
        
        # Text input frame for sequences
        self.sequence_text_frame = ttk.LabelFrame(control_frame, text="Manual Sequence Input", padding=5)
        self.sequence_text_frame.pack(fill='x', pady=(0, 10))
        
        # Add example text for sequences
        sequence_example = ">Sequence1\nATCGATCGATCG\n>Sequence2\nGCTAGCTAGCTA"
        
        self.sequence_text = scrolledtext.ScrolledText(self.sequence_text_frame, height=10, width=50)
        self.sequence_text.pack(fill='both', expand=True, padx=5, pady=5)
        self.sequence_text.insert(tk.END, sequence_example)
        
        # Add confirmation button for sequence input
        sequence_button_frame = ttk.Frame(self.sequence_text_frame)
        sequence_button_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(sequence_button_frame, text="✓ Confirm Sequences", 
                  command=self.confirm_text_input, width=20).pack(side='left', padx=5)
        ttk.Button(sequence_button_frame, text="Clear", 
                  command=self.clear_text_input, width=20).pack(side='left', padx=5)
        
        # Text input frame for distance matrix
        self.matrix_text_frame = ttk.LabelFrame(control_frame, text="Manual Distance Matrix Input (Lower Triangular Format)", padding=5)
        self.matrix_text_frame.pack(fill='x', pady=(0, 10))
        
        # Add example text for distance matrix in lower triangular format
        matrix_example = """Taxon1
Taxon2 0.5
Taxon3 0.7 0.3
Taxon4 0.8 0.4 0.2"""
        
        self.matrix_text = scrolledtext.ScrolledText(self.matrix_text_frame, height=10, width=50)
        self.matrix_text.pack(fill='both', expand=True, padx=5, pady=5)
        self.matrix_text.insert(tk.END, matrix_example)
        
        # Add help text
        help_text = "Format: Each line contains taxon name followed by distances to previous taxa"
        help_label = ttk.Label(self.matrix_text_frame, text=help_text, wraplength=400)
        help_label.pack(padx=5, pady=2)
        
        # Add confirmation button for matrix input
        matrix_button_frame = ttk.Frame(self.matrix_text_frame)
        matrix_button_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(matrix_button_frame, text="✓ Confirm Matrix", 
                  command=self.confirm_matrix_input, width=20).pack(side='left', padx=5)
        ttk.Button(matrix_button_frame, text="Clear", 
                  command=self.clear_matrix_input, width=20).pack(side='left', padx=5)
        
        # File operations
        file_frame = ttk.LabelFrame(control_frame, text="File Operations", padding=5)
        file_frame.pack(fill='x', pady=(0, 10))
        
        ttk.Button(file_frame, text="📁 Load Input File", 
                  command=self.load_file, width=45).pack(pady=2)
        ttk.Button(file_frame, text="🌳 Generate Tree", 
                  command=self.generate_tree, width=45).pack(pady=2)
        
        # Export options
        export_frame = ttk.LabelFrame(control_frame, text="Export Options", padding=5)
        export_frame.pack(fill='x', pady=(0, 10))
        
        ttk.Button(export_frame, text="💾 Save Tree (Newick)", 
                  command=self.save_tree_file, width=45).pack(pady=2)
        ttk.Button(export_frame, text="🖼️ Save Image (PNG)", 
                  command=self.save_image, width=45).pack(pady=2)
        ttk.Button(export_frame, text="📄 Generate Report", 
                  command=self.generate_report, width=45).pack(pady=2)
        ttk.Button(export_frame, text="📊 Save Complete Analysis", 
                  command=self.save_complete_analysis, width=45).pack(pady=2)
        
        # Status and info
        info_frame = ttk.LabelFrame(control_frame, text="Information", padding=5)
        info_frame.pack(fill='both', expand=True)
        
        self.info_text = scrolledtext.ScrolledText(info_frame, height=10, width=30)
        self.info_text.pack(fill='both', expand=True)
        
        # Initial info
        self.update_info("Ready to load sequences or distance matrix.")
        
        # Initialize input type
        self.toggle_input_type()

    def toggle_input_type(self):
        """Toggle visibility of input options based on selected type"""
        if self.option.get() == 'sequences':
            self.sequence_text_frame.pack(fill='x', pady=(0, 10))
            self.matrix_text_frame.pack_forget()
        else:
            self.sequence_text_frame.pack_forget()
            self.matrix_text_frame.pack(fill='x', pady=(0, 10))

    def create_display_panel(self, parent):
        """
        Create the display panel for tree visualization.
        
        Args:
            parent (ttk.Frame): The parent frame for the display panel
        """
        # Tree display area
        display_frame = ttk.LabelFrame(parent, text="Phylogenetic Tree", padding=5)
        display_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Create matplotlib figure
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, display_frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        
        # Initial empty plot
        self.ax.text(0.5, 0.5, 'Load data and generate tree to display here', 
                    ha='center', va='center', transform=self.ax.transAxes, 
                    fontsize=14, alpha=0.5)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.canvas.draw()
        
    def update_info(self, message):
        """
        Update the information panel with a new message.
        
        Args:
            message (str): The message to display
        """
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        self.info_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.info_text.see(tk.END)
        self.root.update()
        
    def confirm_text_input(self):
        """Process and validate the text input sequences"""
        text_input = self.sequence_text.get(1.0, tk.END).strip()
        if not text_input:
            messagebox.showerror("Error", "Please enter sequences in FASTA format")
            return
            
        try:
            # Create a temporary file with the text input
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_file:
                temp_file.write(text_input)
                temp_file_path = temp_file.name
            
            self.sequences = []
            self.seq_ids = []
            
            for record in SeqIO.parse(temp_file_path, "fasta"):
                self.sequences.append(str(record.seq))
                self.seq_ids.append(record.id)
            
            os.unlink(temp_file_path)  # Clean up the temporary file
            
            if not self.sequences:
                messagebox.showerror("Error", "No valid sequences found in the input")
                return
                
            self.update_info(f"Successfully loaded {len(self.sequences)} sequences from text input")
            self.update_info(f"Sequence IDs: {', '.join(self.seq_ids[:5])}{'...' if len(self.seq_ids) > 5 else ''}")
            messagebox.showinfo("Success", f"Successfully loaded {len(self.sequences)} sequences")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to parse sequences: {str(e)}")
            self.update_info(f"Error parsing sequences: {str(e)}")

    def clear_text_input(self):
        """Clear the text input field"""
        self.sequence_text.delete(1.0, tk.END)
        self.sequences = []
        self.seq_ids = []
        self.update_info("Text input cleared")

    def confirm_matrix_input(self):
        """Process and validate the distance matrix input in lower triangular format"""
        text_input = self.matrix_text.get(1.0, tk.END).strip()
        if not text_input:
            messagebox.showerror("Error", "Please enter a distance matrix")
            return
            
        try:
            # Parse the input
            lines = text_input.splitlines()
            if len(lines) < 2:
                raise ValueError("Matrix must have at least two taxa")
                
            # Get taxon names and build matrix
            self.seq_ids = []
            matrix = []
            
            for i, line in enumerate(lines):
                parts = line.strip().split()
                if not parts:
                    continue
                    
                # First part is always the taxon name
                taxon = parts[0]
                self.seq_ids.append(taxon)
                
                # Get distances
                distances = [float(x) for x in parts[1:]]
                
                # Validate number of distances
                if len(distances) != i:
                    raise ValueError(f"Line {i+1}: Expected {i} distance(s) for {taxon}, got {len(distances)}")
                
                # Add distances to matrix
                matrix.append(distances)
            
            if len(self.seq_ids) < 2:
                raise ValueError("Matrix must have at least 2 taxa")
                
            # Create full matrix from lower triangular
            full_matrix = []
            for i in range(len(self.seq_ids)):
                row = [0.0] * (i + 1)  # Initialize with zeros
                for j in range(i):
                    row[j] = matrix[i][j]
                full_matrix.append(row)
            
            # Create distance matrix
            self.distance_matrix = DistanceMatrix(self.seq_ids, full_matrix)
            
            self.update_info(f"Successfully loaded distance matrix for {len(self.seq_ids)} taxa")
            self.update_info(f"Taxa: {', '.join(self.seq_ids[:5])}{'...' if len(self.seq_ids) > 5 else ''}")
            messagebox.showinfo("Success", f"Successfully loaded distance matrix for {len(self.seq_ids)} taxa")
            
        except ValueError as e:
            messagebox.showerror("Error", f"Invalid matrix format: {str(e)}")
            self.update_info(f"Error parsing distance matrix: {str(e)}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to parse distance matrix: {str(e)}")
            self.update_info(f"Error parsing distance matrix: {str(e)}")

    def clear_matrix_input(self):
        """Clear the matrix input field"""
        self.matrix_text.delete(1.0, tk.END)
        self.distance_matrix = None
        self.seq_ids = []
        self.update_info("Matrix input cleared")

    def load_file(self):
        """
        Load input data from a file.
        
        Supports:
        - FASTA files for sequences
        - Text files for distance matrices
        """
        file_path = filedialog.askopenfilename(
            title="Select input file",
            filetypes=[
                ("FASTA files", "*.fasta *.fas *.fa"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )
        
        if not file_path:
            return

        try:
            if self.option.get() == 'sequences':
                self.sequences = []
                self.seq_ids = []
                
                for record in SeqIO.parse(file_path, "fasta"):
                    self.sequences.append(str(record.seq))
                    self.seq_ids.append(record.id)
                
                self.update_info(f"Loaded {len(self.sequences)} sequences from {os.path.basename(file_path)}")
                self.update_info(f"Sequence IDs: {', '.join(self.seq_ids[:5])}{'...' if len(self.seq_ids) > 5 else ''}")
                
            else:  # distance matrix
                with open(file_path, 'r') as f:
                    lines = f.read().splitlines()
                    self.seq_ids = lines[0].split()
                    matrix = [list(map(float, line.split())) for line in lines[1:]]
                    self.distance_matrix = DistanceMatrix(self.seq_ids, matrix)
                
                self.update_info(f"Loaded distance matrix for {len(self.seq_ids)} taxa from {os.path.basename(file_path)}")
                self.update_info(f"Taxa: {', '.join(self.seq_ids[:5])}{'...' if len(self.seq_ids) > 5 else ''}")
                
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {str(e)}")
            self.update_info(f"Error loading file: {str(e)}")

    def compute_distance_matrix(self):
        """
        Compute pairwise distances between sequences.
        
        Returns:
            DistanceMatrix: A matrix of pairwise distances between sequences
        """
        self.update_info("Computing pairwise distances...")
        size = len(self.sequences)
        matrix = [[0] * i for i in range(1, size + 1)]

        total_pairs = sum(range(size))
        current_pair = 0

        for i in range(size):
            for j in range(i):
                current_pair += 1
                if current_pair % 10 == 0:  # Update progress every 10 pairs
                    progress = (current_pair / total_pairs) * 100
                    self.update_info(f"Computing distances... {progress:.1f}% complete")
                
                align = pairwise2.align.globalxx(self.sequences[i], self.sequences[j], one_alignment_only=True)
                score = align[0].score
                max_len = max(len(self.sequences[i]), len(self.sequences[j]))
                distance = 1 - (score / max_len)
                matrix[i][j] = distance

        self.update_info("Distance matrix computation complete.")
        return DistanceMatrix(self.seq_ids, matrix)

    def generate_tree(self):
        """
        Generate a phylogenetic tree using the UPGMA algorithm.
        
        The tree is generated either from sequences or a pre-computed distance matrix,
        depending on the selected input type.
        """
        try:
            if self.option.get() == 'sequences':
                if not self.sequences:
                    messagebox.showerror("Error", "No sequences loaded.")
                    return
                distance_matrix = self.compute_distance_matrix()
            else:
                if not self.distance_matrix:
                    messagebox.showerror("Error", "No distance matrix loaded.")
                    return
                distance_matrix = self.distance_matrix

            self.update_info("Constructing UPGMA tree...")
            constructor = DistanceTreeConstructor()
            self.current_tree = constructor.upgma(distance_matrix)

            # Clear and redraw the plot
            self.ax.clear()
            Phylo.draw(self.current_tree, axes=self.ax)
            self.ax.set_title("UPGMA Phylogenetic Tree", fontsize=16, fontweight='bold')
            
            # Improve plot appearance
            self.fig.tight_layout()
            self.canvas.draw()

            self.update_info("Tree generation complete!")
            self.update_info(f"Tree contains {len(self.seq_ids)} taxa")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate tree: {str(e)}")
            self.update_info(f"Error generating tree: {str(e)}")

    def save_tree_file(self):
        """
        Save the current tree in Newick format.
        """
        if not self.current_tree:
            messagebox.showwarning("Warning", "No tree to save. Generate a tree first.")
            return
            
        file_path = filedialog.asksaveasfilename(
            title="Save tree file",
            defaultextension=".nwk",
            filetypes=[
                ("Newick format", "*.nwk"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )
        
        if file_path:
            try:
                Phylo.write(self.current_tree, file_path, "newick")
                self.update_info(f"Tree saved to {os.path.basename(file_path)}")
                messagebox.showinfo("Success", f"Tree saved successfully to {os.path.basename(file_path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save tree: {str(e)}")

    def save_image(self):
        """
        Save the current tree visualization as an image file.
        """
        if not self.current_tree:
            messagebox.showwarning("Warning", "No tree to save. Generate a tree first.")
            return
            
        file_path = filedialog.asksaveasfilename(
            title="Save tree image",
            defaultextension=".png",
            filetypes=[
                ("PNG files", "*.png"),
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("All files", "*.*")
            ]
        )
        
        if file_path:
            try:
                # Save with high DPI for better quality
                self.fig.savefig(file_path, dpi=300, bbox_inches='tight', 
                               facecolor='white', edgecolor='none')
                self.update_info(f"Image saved to {os.path.basename(file_path)}")
                messagebox.showinfo("Success", f"Image saved successfully to {os.path.basename(file_path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save image: {str(e)}")

    def generate_report(self):
        """
        Generate a detailed analysis report of the current tree.
        """
        if not self.current_tree:
            messagebox.showwarning("Warning", "No tree to report on. Generate a tree first.")
            return
            
        file_path = filedialog.asksaveasfilename(
            title="Save analysis report",
            defaultextension=".txt",
            filetypes=[
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )
        
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    f.write("UPGMA Phylogenetic Tree Analysis Report\n")
                    f.write("=" * 50 + "\n\n")
                    f.write(f"Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"Number of taxa: {len(self.seq_ids)}\n")
                    f.write(f"Input type: {'Sequences' if self.option.get() == 'sequences' else 'Distance Matrix'}\n\n")
                    
                    f.write("Taxa included:\n")
                    for i, taxon in enumerate(self.seq_ids, 1):
                        f.write(f"{i:2d}. {taxon}\n")
                    
                    f.write(f"\nTree in Newick format:\n")
                    f.write(f"{self.current_tree.format('newick')}\n\n")
                    
                    if self.option.get() == 'sequences':
                        f.write("Sequence information:\n")
                        for i, (seq_id, seq) in enumerate(zip(self.seq_ids, self.sequences)):
                            f.write(f"{seq_id}: {len(seq)} bp/aa\n")
                    
                    f.write(f"\nAnalysis method: UPGMA (Unweighted Pair Group Method with Arithmetic Mean)\n")
                    f.write(f"Distance calculation: {'Pairwise sequence alignment' if self.option.get() == 'sequences' else 'User-provided matrix'}\n")
                
                self.update_info(f"Report saved to {os.path.basename(file_path)}")
                messagebox.showinfo("Success", f"Report saved successfully to {os.path.basename(file_path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save report: {str(e)}")

    def save_complete_analysis(self):
        """
        Save a complete analysis package including tree file, image, and report.
        """
        if not self.current_tree:
            messagebox.showwarning("Warning", "No tree to save. Generate a tree first.")
            return
            
        folder_path = filedialog.askdirectory(title="Select folder for complete analysis")
        
        if folder_path:
            try:
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                base_name = f"phylo_analysis_{timestamp}"
                
                # Save tree file
                tree_file = os.path.join(folder_path, f"{base_name}.nwk")
                Phylo.write(self.current_tree, tree_file, "newick")
                
                # Save image
                image_file = os.path.join(folder_path, f"{base_name}.png")
                self.fig.savefig(image_file, dpi=300, bbox_inches='tight', 
                               facecolor='white', edgecolor='none')
                
                # Save report
                report_file = os.path.join(folder_path, f"{base_name}_report.txt")
                with open(report_file, 'w') as f:
                    f.write("UPGMA Phylogenetic Tree Analysis Report\n")
                    f.write("=" * 50 + "\n\n")
                    f.write(f"Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"Number of taxa: {len(self.seq_ids)}\n")
                    f.write(f"Input type: {'Sequences' if self.option.get() == 'sequences' else 'Distance Matrix'}\n\n")
                    
                    f.write("Files generated:\n")
                    f.write(f"- Tree file: {base_name}.nwk\n")
                    f.write(f"- Tree image: {base_name}.png\n")
                    f.write(f"- This report: {base_name}_report.txt\n\n")
                    
                    f.write("Taxa included:\n")
                    for i, taxon in enumerate(self.seq_ids, 1):
                        f.write(f"{i:2d}. {taxon}\n")
                    
                    f.write(f"\nTree in Newick format:\n")
                    f.write(f"{self.current_tree.format('newick')}\n\n")
                    
                    if self.option.get() == 'sequences':
                        f.write("Sequence information:\n")
                        for i, (seq_id, seq) in enumerate(zip(self.seq_ids, self.sequences)):
                            f.write(f"{seq_id}: {len(seq)} bp/aa\n")
                
                self.update_info(f"Complete analysis saved to {folder_path}")
                self.update_info(f"Files: {base_name}.nwk, {base_name}.png, {base_name}_report.txt")
                messagebox.showinfo("Success", f"Complete analysis saved successfully!\n\nFiles saved:\n- {base_name}.nwk\n- {base_name}.png\n- {base_name}_report.txt")
                
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save complete analysis: {str(e)}")

# Run the application
if __name__ == "__main__":
    root = tk.Tk()
    app = PhyloApp(root)
    root.mainloop()