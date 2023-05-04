import requests
import json
import os
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils import molecular_weight, ProtParam
from Bio import SeqIO
#NOTE: DEPENDENCIES
#You MUST run using python3.10, as 3.11 does not support current 'numba' version
#In the environment, the following dependencies may need to be manually installed:
#umap-learn, pot, matplotlib, scikit-learn, pandas, prody, rdkit, chembl_webresource_client, biopython, numba, scipy, numpy, requests

category_space = ["Alzheimer's", "FTD", "Parkinson's", "Huntington's", "Cancer", "Sirtuins", "HIV", "Other"]
#perfect_category_space = {"Dementia": ["Frontotemporal Dementia", "Alzheimer's Disease", "Lewy Body Dementia", "Vascular Dementia", "Mixed Dementia", "Parkinson's Disease Dementia", "Huntington's Disease", "Creutzfeldt-Jakob Disease", "Normal Pressure Hydrocephalus", "Wernicke-Korsakoff Syndrome", "Traumatic Brain Injury", "HIV-Associated Neurocognitive Disorder", "Chronic Traumatic Encephalopathy", "Down Syndrome", "Mild Cognitive Impairment", "Posterior Cortical Atrophy", "Primary Progressive Aphasia", "Progressive Supranuclear Palsy", "Corticobasal Degeneration", "Multiple Sclerosis", "Prion Disease", "Amyotrophic Lateral Sclerosis", "Frontotemporal Lobar Degeneration", "Dementia with Lewy Bodies", "Parkinson's Disease", "Parkinsonism", "Corticobasal Syndrome", "Progressive Nonfluent Aphasia", "Semantic Dementia", "Posterior Cortical Atrophy", "Primary Progressive Aphasia", "Progressive Supranuclear Palsy", "Corticobasal Degeneration", "Multiple Sclerosis", "Prion Disease", "Amyotrophic Lateral Sclerosis", "Frontotemporal Lobar Degeneration", "Dementia with Lewy Bodies", "Parkinson's Disease", "Parkinsonism", "Corticobasal Syndrome", "Progressive Nonfluent Aphasia", "Semantic Dementia", "Posterior Cortical Atrophy", "Primary Progressive Aphasia", "Progressive Supranuclear Palsy", "Corticobasal Degeneration", "Multiple Sclerosis", "Prion Disease", "Amyotrophic Lateral Sclerosis", "Frontotemporal Lobar Degeneration", "Dementia with Lewy Bodies", "Parkinson's Disease", "Parkinsonism", "Corticobasal Syndrome", "Progressive Nonfluent Aphasia", "Semantic Dementia", "Posterior Cortical Atrophy", "Primary Progressive Aphasia", "Progressive Supranuclear Palsy", "Corticobasal Degeneration", "Multiple Sclerosis", "Prion Disease", "Amyotrophic Lateral Sclerosis", "Frontotemporal Lobar Degeneration", "Dementia with Lewy Bodies", "Parkinson's Disease", "Parkinsonism", "Corticobasal Syndrome","]}

unique_ligands = ["NAG", "GNT", "C6K", "PEG", "C6H", "4I4", "GOL", "TAM", "PLP",
"IOD", "MPD","IHP", "DMS", "ILE", "STL", "IPA", "PIT",
"AR6", "UNX", "L3D", "LFO", "EPE", "7PE", "9A5", "PLP", "142"]

cmd = input("Enter a unique 3-letter ligand code: ").upper()
if cmd in unique_ligands:
    print("Ligand is unique")
else:
    print("Ligand is not in Ligand-Fang/MedStructLib database")

ligand_code = cmd

def fetch_pdb_entry(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"No data was fetched from RCSB PDB for {pdb_id}")
        return None

def fetch_fasta(pdb_id):
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.content.decode('utf-8')
    else:
        print(f"No fasta data was fetched from RCSB PDB for {pdb_id}")
        return None

def extract_relevant_features(data, sequence):
    relevant_features = {}

    entry_info = data.get("rcsb_entry_info", {})
    relevant_features["molecular_weight"] = entry_info.get("molecular_weight")
    relevant_features["deposited_polymer_monomer_count"] = entry_info.get("deposited_polymer_monomer_count")
    relevant_features["deposited_modeled_polymer_monomer_count"] = entry_info.get("deposited_modeled_polymer_monomer_count")
    relevant_features["deposited_unmodeled_polymer_monomer_count"] = entry_info.get("deposited_unmodeled_polymer_monomer_count")
    relevant_features["disulfide_bond_count"] = entry_info.get("disulfide_bond_count")
    relevant_features["inter_mol_covalent_bond_count"] = entry_info.get("inter_mol_covalent_bond_count")
    relevant_features["inter_mol_metalic_bond_count"] = entry_info.get("inter_mol_metalic_bond_count")
    relevant_features["polymer_entity_count_protein"] = entry_info.get("polymer_entity_count_protein")
    seq_length = len(sequence)
    

    # pp = ProtParam.ProteinAnalysis(sequence)
    # charge_distribution =

    relevant_features["sequence_length"] = seq_length
    #relevant_features["charge_distribution"] = charge_distribution

    return relevant_features

pdb_ids = ["7U1O", "7U11", "1QTI", "6EZG", "4BCB", "4BCD", "6C95", "7JQR", "3T4G", "6GHX", "4HDA", "4HD8", "3PKI",
           "6NUJ", "4TSX", "2P5E", "2P5W", "2PYE", "5NUU", "1JS6", "1JS3"]
protein_data = {}

for pdb_id in pdb_ids:
    entry_data = fetch_pdb_entry(pdb_id)
    fasta_seq = fetch_fasta(pdb_id)

    if fasta_seq:
        relevant_features = extract_relevant_features(entry_data, fasta_seq)
        protein_data[pdb_id] = relevant_features

#print(protein_data)
print("Protein data:")
for pdb_id, relevant_features in protein_data.items():
    print(f"{pdb_id}:")
    for feature, value in relevant_features.items():
        print(f"  {feature}: {value}")

#NOTE: Fetch and organize ligand data
import requests
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Crippen
from Bio import SeqIO
from Bio.PDB import PDBParser, DSSP
import prody
import subprocess
from Bio.PDB import PDBParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBList


#NOTE: Med_Struct_Lib (Ligand Fang) is a program that I am working on which will
#search for a receptor or enzyme related to dementia, for example, given a relevant PDB 'unique ligand'.
#The 'Feature Transport' algorithm performs optimal transport on the dimensionality-reduced feature spaces for the ligand and the receptor.
#The program will then return a truncated list of receptors or enzymes that are most compatible with the algorithm-selected proteins.

import numpy as np


def calculate_descriptors(mol):
    descriptors = {
        'Molecular weight': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'Number of H-bond donors': Descriptors.NumHDonors(mol),
        'Number of H-bond acceptors': Descriptors.NumHAcceptors(mol),
        'Number of rotatable bonds': Descriptors.NumRotatableBonds(mol),
        'Number of rings': Descriptors.RingCount(mol),
        'Topological polar surface area (TPSA)': rdMolDescriptors.CalcTPSA(mol),
        'Molar refractivity': Crippen.MolMR(mol),
        #'Number of atoms': mol.GetNumAtoms(),
        #'Electrostactic charge': Chem.GetFormalCharge(mol),
        'Number of heavy atoms': mol.GetNumHeavyAtoms()


        #'Molecular fingerprint': Chem.RDKFingerprint(mol),
        
    }
    
    return descriptors


def fetch_smiles_from_chembl(ligand_code):
    molecule = new_client.molecule
    results = molecule.search(ligand_code)
    
    if results:
        return results[0]['molecule_structures']['canonical_smiles']
    else:
        return None

def fetch_smiles_from_pubchem(ligand_code):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand_code}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']
    else:
        return None



#NOTE: This is a toy dataset being used as I do not have the capacity as a 
# beginner to effectively integrate an entire protein library of SQL into
#  my code. I will just be using
# toy data, as in, procured data, to test my code. The data is genuine; however,
#it cannot be used to test ALL possible ligands. It is a small sample size.

#The data is not quite random. There exists a non-implemented 'category space' which is mutual among both ligand and protein.   
#There are several dementia categories, a sirtuin category, a cancer category, an HIV category, and some outliers.
#The category space is not implemented in this code, but it will be implemented in the future.

#If Ligandfang/MedStructLib's feature extraction goes wrong and the evidence points toward an erroneous protein-ligand match,
#It would be helpful to see that at least the category space is aligned properly. Why? To show similar proteins and ligand
#matches exist in the same category space. This is a good sanity check.


#NOTE: The following code is a toy dataset. It is not random. It is procured data.
#unique_ligands = ["NAG", "GNT", "C6K", "PEG", "C6H", "4I4", "GOL", "TAM", "PLP", "IOD", "MPD","IHP"]
#NOTE: To decrease model complexity, decrease category and feature space by removing proteins from this list; or, uncomment above.
#and comment out the below assignment by dragging the 3 lines below and commenting out with Ctrl + / or Cmd + /.
unique_ligands = ["NAG", "GNT", "C6K", "PEG", "C6H", "4I4", "GOL", "TAM", "PLP",
"IOD", "MPD","IHP", "DMS", "ILE", "STL", "IPA", "PIT",
"AR6", "UNX", "L3D", "LFO", "EPE", "7PE", "9A5", "PLP", "142"] #NOTE: Start by adding more ligand descriptors, then add more ligands, then add more proteins

ligand_data = {}
original_ligand_data = ligand_data.copy()

for ligand_code in unique_ligands:
    smiles = None

    try:
        smiles = fetch_smiles_from_chembl(ligand_code)
    except Exception as e:
        print("Error fetching SMILES from ChEMBL:", e)

    if smiles is None:
        try:
            smiles = fetch_smiles_from_pubchem(ligand_code)
        except Exception as e:
            print("Error fetching SMILES from PubChem:", e)

    if smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        descriptors = calculate_descriptors(mol)
        ligand_data[ligand_code] = descriptors

print("Ligand data:")
print()
print()
for ligand_code, descriptors in ligand_data.items():
    print()
    print(f"{ligand_code}:")
    print()
    for descriptor, value in descriptors.items():
        print(f"{descriptor}: {value}")

print()
print()
print()
print("Building feature landscapes, please wait...")



#NOTE: Preprocessing Feature Spaces
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

# Convert ligand_data dictionary to DataFrame
ligand_df = pd.DataFrame.from_dict(ligand_data, orient='index')

# Normalize ligand data
scaler_ligand = MinMaxScaler()
scaled_ligand_data = scaler_ligand.fit_transform(ligand_df)
scaled_ligand_df = pd.DataFrame(scaled_ligand_data, index=ligand_df.index, columns=ligand_df.columns)

# Convert protein_data dictionary to DataFrame
protein_df = pd.DataFrame.from_dict(protein_data, orient='index')

# Normalize protein data
scaler_protein = MinMaxScaler()
scaled_protein_data = scaler_protein.fit_transform(protein_df)
scaled_protein_df = pd.DataFrame(scaled_protein_data, index=protein_df.index, columns=protein_df.columns)

print("Scaled ligand data:")
print(scaled_ligand_df)

print("\nScaled protein data:")
print(scaled_protein_df)

from sklearn.decomposition import PCA

# Perform PCA on both the scaled ligand and protein data arrays
n_components = 2  # Adjust the number of components based on your desired dimensionality
pca = PCA(n_components=n_components)
pca_ligand = pca.fit_transform(scaled_ligand_data)
pca_protein = pca.fit_transform(scaled_protein_data)

# Convert the PCA results to DataFrames with appropriate index
pca_ligand_df = pd.DataFrame(pca_ligand, index=list(ligand_data.keys()), columns=[f'PC{i+1}' for i in range(n_components)])
pca_protein_df = pd.DataFrame(pca_protein, index=list(protein_data.keys()), columns=[f'PC{i+1}' for i in range(n_components)])

print("PCA results for ligand features:")
print(pca_ligand_df)
print("\nPCA results for protein features:")
print(pca_protein_df)

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

# Plot PCA results for ligands (blue dots)
ax.scatter(pca_ligand_df['PC1'], pca_ligand_df['PC2'], color='blue', label='Ligands')
for i, txt in enumerate(pca_ligand_df.index):
    ax.annotate(txt, (pca_ligand_df['PC1'].iloc[i], pca_ligand_df['PC2'].iloc[i]), fontsize=8)

# Plot PCA results for proteins (red dots)
ax.scatter(pca_protein_df['PC1'], pca_protein_df['PC2'], color='red', label='Proteins')
for i, txt in enumerate(pca_protein_df.index):
    ax.annotate(txt, (pca_protein_df['PC1'].iloc[i], pca_protein_df['PC2'].iloc[i]), fontsize=8)

# Customize plot appearance
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
#n_components = 3, so PC3 is the third component
#ax.set_label('PC3')
ax.legend()
plt.title('PCA results for ligand and protein features')
plt.show()

#NOTE: Add cancer proteins and ligands, and separate into 4-5 (ligand) categories:
# # Dementia, Cancer, Longevity, Other, and Protein

#Feature Mapping
import numpy as np
from scipy.spatial.distance import cdist
from matplotlib.colors import LinearSegmentedColormap

# Calculate pairwise Euclidean distances
distance_matrix = cdist(pca_ligand_df, pca_protein_df, metric='euclidean')

# Normalize distances to range [0, 1]
normalized_distances = (distance_matrix - distance_matrix.min()) / (distance_matrix.max() - distance_matrix.min())

# Create a feature colormap
cmap = LinearSegmentedColormap.from_list("skyblue_yellow", ["skyblue", "yellow"])

# Generate the colored 2D feature map
fig, ax = plt.subplots()
im = ax.imshow(normalized_distances, cmap=cmap)
cbar = fig.colorbar(im, ax=ax)
plt.xticks(np.arange(len(pca_protein_df)), pca_protein_df.index, rotation=90)
plt.yticks(np.arange(len(pca_ligand_df)), pca_ligand_df.index)
plt.xlabel("Proteins")
plt.ylabel("Ligands")
plt.title("Feature Map of Ligand-Protein Similarities")
plt.show()

#Ligand-Protein Matcher
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import numpy as np
from sklearn.decomposition import PCA



# Perform PCA on ligand and protein features
pca = PCA(n_components=2)
pca_ligand = pca.fit_transform(pca_ligand_df)
pca_protein = pca.fit_transform(pca_protein_df)

# Concatenate PCA results into a single variable
pca_combined = np.concatenate((pca_ligand, pca_protein), axis=0)

# Now you should be able to use pca_combined
ligand_points = pca_combined[:len(scaled_ligand_data)]

# Assuming you already have pca_combined (output of PCA)
ligand_points = pca_combined[:len(scaled_ligand_data)]
protein_points = pca_combined[len(scaled_ligand_data):]

# Calculate Euclidean distances
euclidean_distances = cdist(ligand_points, protein_points, metric='euclidean')

# Calculate cosine similarities
cosine_similarities = cosine_similarity(ligand_points, protein_points)

import ot

# Optimal transport distances
ot_distances = np.zeros_like(euclidean_distances)
for i in range(len(ligand_points)):
    for j in range(len(protein_points)):
        ot_distances[i, j] = np.linalg.norm(ligand_points[i] - protein_points[j])

# Normalize scores between 0 and 1
normalized_euclidean = (euclidean_distances - euclidean_distances.min()) / (euclidean_distances.max() - euclidean_distances.min())
normalized_cosine = (cosine_similarities - cosine_similarities.min()) / (cosine_similarities.max() - cosine_similarities.min())
normalized_ot = (ot_distances - ot_distances.min()) / (ot_distances.max() - ot_distances.min())

# Calculate composite score
composite_score = 0.3 * (1 - normalized_euclidean) + 0.4 * normalized_cosine + 0.3 * (1 - normalized_ot)

# Identify potential matches based on the composite score
similarity_threshold = 0.7
potential_matches = np.where(composite_score > similarity_threshold)

# Print out the potential matches
print("Potential matches via PCA-derived analysis:")
for match in zip(*potential_matches):
    ligand_index, protein_index = match
    ligand_name = list(ligand_data.keys())[ligand_index]
    protein_name = list(protein_data.keys())[protein_index]
    print(f"Potential match: {ligand_name} - {protein_name} (Composite score: {composite_score[ligand_index][protein_index]})")
# for match in zip(*potential_matches):
#     ligand_index, protein_index = match
#     ligand_name = list(original_ligand_data.keys())[ligand_index]
#     protein_name = list(scaled_protein_data.keys())[protein_index]
#     print(f"Potential match: {ligand_name} - {protein_name} (Composite score: {composite_score[ligand_index][protein_index]})")

#NOTE: UMAP Visualization & Feature Mapping with UMAP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from umap import UMAP

#We need to manipulate the PCA data to be able to use it in UMAP
#Why? Because the data must be in a df format with similar dimensions
#PCA is sometimes done before UMAP to reduce the dimensions of the data
#Our feature space is around ~20-dimensional and grows with each addition to features,
#, so it makes sense to use PCA followed by UMAP
import numpy as np
import matplotlib.pyplot as plt
from umap import UMAP

# Apply UMAP to pca_combined data
umap = UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean')
umap_embeddings = umap.fit_transform(pca_combined)

# Separate the ligand and protein embeddings
ligand_embeddings = umap_embeddings[:len(pca_ligand)]
protein_embeddings = umap_embeddings[len(pca_ligand):]

# Plot ligand and protein embeddings
plt.figure(figsize=(10, 8))
plt.scatter(ligand_embeddings[:, 0], ligand_embeddings[:, 1], c='blue', label='Ligands', alpha=0.6)
plt.scatter(protein_embeddings[:, 0], protein_embeddings[:, 1], c='red', label='Proteins', alpha=0.6)

# Label each ligand point
for i, label in enumerate(pca_ligand_df.index):  # Replace ligand_labels with your actual list of ligand labels
    plt.annotate(label, (ligand_embeddings[i, 0], ligand_embeddings[i, 1]), fontsize=9, alpha=0.75)

# Label each protein point
for i, label in enumerate(pca_protein_df.index):  # Replace protein_labels with your actual list of protein labels
    plt.annotate(label, (protein_embeddings[i, 0], protein_embeddings[i, 1]), fontsize=9, alpha=0.75)

plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend()
plt.title('UMAP Embeddings of Ligands and Proteins')

plt.show()

#Building UMAP feature map
print("Building UMAP feature map...")
import umap

combined_data_pca = np.concatenate([pca_ligand_df, pca_protein_df], axis=0)
umap_reducer = umap.UMAP()
umap_embedding = umap_reducer.fit_transform(combined_data_pca)

ligand_points = umap_embedding[:len(pca_ligand_df)]
protein_points = umap_embedding[len(pca_ligand_df):]
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

normalized_distances = (distance_matrix - distance_matrix.min()) / (distance_matrix.max() - distance_matrix.min())
cmap = LinearSegmentedColormap.from_list("cyan_yellow", ["cyan", "yellow"])

fig, ax = plt.subplots()
im = ax.imshow(normalized_distances, cmap=cmap)
cbar = fig.colorbar(im, ax=ax)
plt.xticks(np.arange(len(pca_protein_df)), pca_protein_df.index, rotation=90)
plt.yticks(np.arange(len(pca_ligand_df)), pca_ligand_df.index)
plt.xlabel("Proteins")
plt.ylabel("Ligands")
plt.title("Feature Map of Ligand-Protein Similarities (UMAP)")
plt.show()
from scipy.spatial.distance import cdist, cosine
import ot

# Euclidean distances
euclidean_distances = cdist(ligand_points, protein_points, metric='euclidean')

# Cosine similarities
cosine_similarities = 1 - cdist(ligand_points, protein_points, metric='cosine')

# Wasserstein distances
import numpy as np
import ot

# You already have ligand_embeddings and protein_embeddings from the previous code

# Calculate the cost matrix between ligand and protein embeddings
M = ot.dist(ligand_embeddings, protein_embeddings, metric='euclidean')

# Compute the Wasserstein distances between ligand and protein embeddings
import numpy as np
import ot

# You already have ligand_embeddings and protein_embeddings from the previous code

# Compute the Wasserstein distances between ligand and protein embeddings
wasserstein_distances = np.zeros((len(ligand_embeddings), len(protein_embeddings)))

for i, ligand_point in enumerate(ligand_embeddings):
    for j, protein_point in enumerate(protein_embeddings):
        # Calculate the cost matrix between the individual ligand and protein point pair
        M_ij = ot.dist(np.array([ligand_point]), np.array([protein_point]), metric='euclidean')
        wasserstein_distances[i, j] = ot.emd2(np.array([1]), np.array([1]), M_ij)

# Example: Retrieve the Wasserstein distance between the first ligand and the first protein
wasserstein_distance_example = wasserstein_distances[0, 0]
print("Wasserstein distance between the first ligand and the first protein:", wasserstein_distance_example)



#Like last time, write code that iterates through each ligand and protein pair and prints out the composite score
# for each pair. The output should look like this:
#Potential match: IHP - 6EZG (Composite score: 0.9294559015313324)

# Print out the potential matches
print("Potential matches via UMAP-derived analysis:")
for match in zip(*np.where(composite_score < 0.5)):
    ligand_index, protein_index = match
    ligand_name = list(ligand_data.keys())[ligand_index]
    protein_name = list(protein_data.keys())[protein_index]
    print(f"Potential match: {ligand_name} - {protein_name} (Composite score: {composite_score[ligand_index][protein_index]})")

# Print out the matches with the lowest composite scores of 'cmd', the variable for the input ligand
# Find the index of the input ligand in the ligand list
ligand_index = unique_ligands.index(cmd)

# Get the composite scores for the input ligand
input_ligand_scores = composite_score[ligand_index]

# Combine protein names and scores into a list of tuples
protein_score_tuples = list(zip(protein_name, input_ligand_scores))

# Sort the tuples based on the composite scores (ascending order)
sorted_protein_score_tuples = sorted(protein_score_tuples, key=lambda x: x[1])

# Print the matches with the lowest composite scores for the input ligand
print(f"Matches for ligand {cmd} (sorted by composite scores):")
for protein, score in sorted_protein_score_tuples:
    print(f"{protein}: {score}")






#NOTE: There needs to be more features for the ligands and proteins to be able to match them up
#NOTE: Add cancer proteins and ligands, and separate into 4-5 (ligand) categories:
# # Dementia, Cancer, Longevity, Other, and Protein
# #Consider showing the lowest composite scores for both PCA and UMAP (top 4)


