import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from Bio import SeqIO, Phylo, pairwise2
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import datetime
from pathlib import Path

class PhyloApp:
    def __init__(self, root):
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
        # Main control frame
        control_frame = ttk.LabelFrame(parent, text="Controls", padding=10)
        control_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Input options
        input_frame = ttk.LabelFrame(control_frame, text="Input Type", padding=5)
        input_frame.pack(fill='x', pady=(0, 10))
        
        ttk.Radiobutton(input_frame, text="DNA/Protein Sequences (FASTA)", 
                       variable=self.option, value='sequences').pack(anchor='w')
        ttk.Radiobutton(input_frame, text="Distance Matrix", 
                       variable=self.option, value='matrix').pack(anchor='w')
        
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
        
    def create_display_panel(self, parent):
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
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        self.info_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.info_text.see(tk.END)
        self.root.update()
        
    def load_file(self):
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