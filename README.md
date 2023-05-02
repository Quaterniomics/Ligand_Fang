![image](https://user-images.githubusercontent.com/111631655/235791690-829bb78a-086e-493f-b063-7369d0105fe2.png)
Ligand Fang
Ligand Fang is a Python-based tool designed to identify potential protein targets for a given ligand, focusing on proteins related to dementia. The tool performs data preprocessing, feature extraction, and similarity scoring to compare the input ligand with a dataset of known protein structures. Ligand Fang provides a similarity map that visualizes the relationships between the ligand and the proteins, facilitating the identification of promising protein targets for further study.

Features
Fetches ligand and protein data from external sources (ChEMBL and PDB)
Calculates molecular descriptors using RDKit
Preprocesses data and normalizes feature values
Applies dimensionality reduction techniques (PCA and UMAP)
Computes similarity scores using Euclidean and Wasserstein distances
Visualizes results with a color-coded similarity map
Dependencies
Ligand Fang requires the following Python packages:

pandas
numpy
requests
BeautifulSoup4
RDKit
scikit-learn
UMAP
matplotlib
To install these dependencies, run the following command:

<pre>
```
pip install pandas numpy requests beautifulsoup4 rdkit scikit-learn umap-learn matplotlib
```
</pre>

Usage
Clone the Ligand Fang repository or download the ligandfang.py script.
Ensure that all dependencies are installed (see "Dependencies" above).
Run the script with the following command:

<pre>
```
python ligandfang.py
```
</pre>

When prompted, enter a 3-digit ligand code (e.g., "GNT") or enter "ran" to select a random ligand.
The script will fetch the relevant data, perform preprocessing and analysis, and display a similarity map showing the relationships between the input ligand and the protein dataset.
Notes
As of May 2023, Ligand Fang does not consider ionized ligands.

For additional details on the methods and techniques used in Ligand Fang, please refer to the comments and documentation within the ligandfang.py script.

