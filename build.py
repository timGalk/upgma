import os
import sys
import shutil
import subprocess
from pathlib import Path
import site

def get_venv_path():
    """Get the path to the virtual environment."""
    # Check if we're running in a virtual environment
    if hasattr(sys, 'real_prefix') or (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
        return sys.prefix
    return None

def get_biopython_data_path():
    """Get the path to Biopython data files."""
    venv_path = get_venv_path()
    if venv_path:
        # Use virtual environment's site-packages
        if os.name == 'nt':  # Windows
            site_packages = os.path.join(venv_path, 'Lib', 'site-packages')
        else:  # Unix/MacOS
            site_packages = os.path.join(venv_path, 'lib', f'python{sys.version_info.major}.{sys.version_info.minor}', 'site-packages')
    else:
        # Use system site-packages
        site_packages = site.getsitepackages()[0]
    
    # Construct path to Biopython data
    bio_path = os.path.join(site_packages, 'Bio')
    return bio_path

def clean_build_dirs():
    """Clean build and dist directories."""
    dirs_to_clean = ['build', 'dist']
    for dir_name in dirs_to_clean:
        if os.path.exists(dir_name):
            print(f"Cleaning {dir_name} directory...")
            shutil.rmtree(dir_name)

def create_executable():
    """Create executable using PyInstaller."""
    print("Creating executable...")
    
    # Get Biopython data path
    bio_path = get_biopython_data_path()
    
    # Get Python executable path
    python_exe = sys.executable
    
    # PyInstaller command
    cmd = [
        python_exe,  # Use the virtual environment's Python
        '-m', 'PyInstaller',
        '--name=UPGMA_Tree_Builder',
        '--windowed',  # No console window
        '--onefile',   # Single executable
        '--icon=icon.ico' if os.path.exists('icon.ico') else '',
        '--add-data=README.md;.',  # Include README
        f'--add-data={os.path.join(bio_path, "Align", "substitution_matrices", "data")};Bio/Align/substitution_matrices/data',  # Include substitution matrices
        f'--add-data={os.path.join(bio_path, "Align", "substitution_matrices")};Bio/Align/substitution_matrices',  # Include substitution matrices directory
        '--clean',     # Clean PyInstaller cache
        'App.py'
    ]
    
    # Remove empty arguments
    cmd = [arg for arg in cmd if arg]
    
    # Run PyInstaller
    subprocess.run(cmd, check=True)

def create_version_file():
    """Create version file with build information."""
    version_info = f"""# Version information
VSVersionInfo(
  ffi=FixedFileInfo(
    filevers=(1, 0, 0, 0),
    prodvers=(1, 0, 0, 0),
    mask=0x3f,
    flags=0x0,
    OS=0x40004,
    fileType=0x1,
    subtype=0x0,
    date=(0, 0)
  ),
  kids=[
    StringFileInfo([
      StringTable(
        u'040904B0',
        [StringStruct(u'CompanyName', u'UPGMA Tree Builder'),
         StringStruct(u'FileDescription', u'UPGMA Phylogenetic Tree Generator'),
         StringStruct(u'FileVersion', u'1.0.0'),
         StringStruct(u'InternalName', u'UPGMA_Tree_Builder'),
         StringStruct(u'LegalCopyright', u'Copyright (c) 2024'),
         StringStruct(u'OriginalFilename', u'UPGMA_Tree_Builder.exe'),
         StringStruct(u'ProductName', u'UPGMA Tree Builder'),
         StringStruct(u'ProductVersion', u'1.0.0')])
    ]),
    VarFileInfo([VarStruct(u'Translation', [1033, 1200])])
  ]
)
"""
    with open('version_info.txt', 'w') as f:
        f.write(version_info)

def main():
    """Main build process."""
    print("Starting build process...")
    
    # Check if we're in a virtual environment
    venv_path = get_venv_path()
    if venv_path:
        print(f"Using virtual environment: {venv_path}")
    else:
        print("Warning: Not running in a virtual environment!")
        response = input("Do you want to continue anyway? (y/n): ")
        if response.lower() != 'y':
            print("Build cancelled.")
            return
    
    # Clean previous builds
    clean_build_dirs()
    
    # Create version file
    create_version_file()
    
    # Create executable
    create_executable()
    
    print("\nBuild completed successfully!")
    print("Executable can be found in the 'dist' directory.")

if __name__ == "__main__":
    main() 