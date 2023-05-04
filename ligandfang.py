import requests
import json
import os
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils import molecular_weight, ProtParam
from Bio import SeqIO
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
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from umap import UMAP
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from matplotlib.colors import LinearSegmentedColormap
import ot
import random
#NOTE: DEPENDENCIES NOTE:
#You MUST run using python3.10, as 3.11 does not support current 'numba' version
#In the environment, the following dependencies may need to be manually installed:
#pubchempy, umap-learn, pot, matplotlib, scikit-learn, pandas, prody, rdkit, chembl_webresource_client, biopython, numba, scipy, numpy, requests

unique_ligands = ["NAG", "GNT", "C6K", "PEG", "C6H", "4I4", "GOL", "TAM", "PLP",
"IOD", "MPD","IHP", "DMS", "ILE", "STL", "IPA", "PIT",
"AR6", "UNX", "L3D", "LFO", "EPE", "7PE", "9A5", "PLP", "142"] 
pdb_ids = ["7U1O", "7U11", "1QTI", "6EZG", "4BCB", "4BCD", "6C95", "7JQR", "3T4G", "6GHX", "4HDA", "4HD8", "3PKI",
      "6NUJ", "4TSX", "2P5E", "2P5W", "2PYE", "5NUU", "1JS6", "1JS3"]
category_space = ["Alzheimer's", "FTD", "Parkinson's", "Huntington's", "Cancer", "Sirtuins", "HIV", "Other"]
#perfect_category_space = {"Dementia": ["Frontotemporal Dementia", "Alzheimer's Disease", "Lewy Body Dementia", "Vascular Dementia", "Mixed Dementia", "Parkinson's Disease Dementia", "Huntington's Disease", "Creutzfeldt-Jakob Disease", "Normal Pressure Hydrocephalus", "Wernicke-Korsakoff Syndrome", "Traumatic Brain Injury", "HIV-Associated Neurocognitive Disorder", "Chronic Traumatic Encephalopathy", "Down Syndrome", "Mild Cognitive Impairment", "Posterior Cortical Atrophy", "Primary Progressive Aphasia", "Progressive Supranuclear Palsy", "Corticobasal Degeneration", "Multiple Sclerosis", "Prion Disease", "Amyotrophic Lateral Sclerosis", "Frontotemporal Lobar Degeneration", "Dementia with Lewy Bodies", "Parkinson's Disease", "Parkinsonism", "Corticobasal Syndrome", "Progressive Nonfluent Aphasia", "Semantic Dementia", "Posterior Cortical Atrophy", "Primary Progressive Aphasia", "Progressive Supranuclear Palsy", "Corticobasal Degeneration", "Multiple Sclerosis", "Prion Disease", "Amyotrophic Lateral Sclerosis", "Frontotemporal Lobar Degeneration", "Dementia with Lewy Bodies", "Parkinson's Disease", "Parkinsonism", "Corticobasal Syndrome", "Progressive Nonfluent Aphasia", "Semantic Dementia", "Posterior Cortical Atrophy", "Primary Progressive Aphasia", "Progressive Supranuclear Palsy", "Corticobasal Degeneration", "Multiple Sclerosis", "Prion Disease", "Amyotrophic Lateral Sclerosis", "Frontotemporal Lobar Degeneration", "Dementia with Lewy Bodies", "Parkinson's Disease", "Parkinsonism", "Corticobasal Syndrome", "Progressive Nonfluent Aphasia", "Semantic Dementia", "Posterior Cortical Atrophy", "Primary Progressive Aphasia", "Progressive Supranuclear Palsy", "Corticobasal Degeneration", "Multiple Sclerosis", "Prion Disease", "Amyotrophic Lateral Sclerosis", "Frontotemporal Lobar Degeneration", "Dementia with Lewy Bodies", "Parkinson's Disease", "Parkinsonism", "Corticobasal Syndrome","]}

#NOTE: Med_Struct_Lib (Ligand Fang) is a program that I am working on which will
#search for a receptor or enzyme related to dementia, for example, given a relevant PDB 'unique ligand'.
#The 'Feature Transport' algorithm performs optimal transport on the dimensionality-reduced feature spaces for the ligand and the receptor.
#The program will then return a truncated list of receptors or enzymes that are most compatible with the algorithm-selected proteins.

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

print("Welcome to Ligand Fang/MedStructLib. This program will search for a receptor or enzyme related to dementia, for example, given a relevant PDB 'unique ligand'.")
print("As of May 2023, Ligand Fang does not consider ionized ligands.")

select = input("Enter 3-digit ligand code. Choose test protein from toy data, such as '9A5', 'NAG' or 'GNT': ").upper()
#select = input("Enter 3-digit ligand code. Enter 'ran' for random ligand: ").upper()


if select == "ran" or select == 'RAN':
    ligand_code = random.choice(unique_ligands)
    mol = ligand_code
    print("Random ligand selected: " + mol)
elif select not in unique_ligands:
    print("Invalid ligand code, please use ligand code from toy data structure 'unique_ligands' in .py file. Exiting program.")
    exit()
else:
    ligand_code = select
    print("Ligand selected: " + select)

"""STEP 1: Fetch PDB entry and FASTA data from RCSB PDB."""


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
    relevant_features["sequence_length"] = seq_length
    #NOTE: LigandFangv2 will consider more descriptors with greater variance,
    #such as: RSA, ASA, ADMET, and other quality features. 
    return relevant_features

protein_data = {}
"""The following loop iterates through our ~30-protein database and utilizes already-established libraries to pull relevant features from the PDB database."""
for pdb_id in pdb_ids:
    entry_data = fetch_pdb_entry(pdb_id)
    fasta_seq = fetch_fasta(pdb_id)

    if fasta_seq:
        relevant_features = extract_relevant_features(entry_data, fasta_seq)
        protein_data[pdb_id] = relevant_features

#NOTE: The following is not necessary for the program to run, but it is useful for debugging purposes.
#Run 'printguts_draft.py' to see the guts of the program.
#print(protein_data)
"""
print("Protein data:")
for pdb_id, relevant_features in protein_data.items():
    print(f"{pdb_id}:")
    for feature, value in relevant_features.items():
        print(f"  {feature}: {value}")
        """

"""NOTE: STEP 2: Feature Extraction. This is the most important step of the program.
The following is the same as step 1, just for the select ligand."""

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

ligand_code = select #if this is not helpful, delete this line
# Fetch SMILES string for the selected ligand
ligand_smiles = fetch_smiles_from_chembl(ligand_code) or fetch_smiles_from_pubchem(ligand_code)

#Step 3: Preprocessing
def preprocess_data(protein_data, ligand_code, ligand_smiles):
    # Step 3A: Convert protein_data dictionary into a Pandas DataFrame
    protein_df = pd.DataFrame.from_dict(protein_data, orient='index')

    # Step 3B: Calculate descriptors for the selected ligand and convert it into a DataFrame
    mol = Chem.MolFromSmiles(ligand_smiles)
    ligand_descriptors = calculate_descriptors(mol)
    ligand_df = pd.DataFrame(ligand_descriptors, index=[ligand_code])

    # Step 3C: Ensure that both DataFrames have the same column names
    common_columns = list(protein_df.columns)
    ligand_df.columns = common_columns

    # Step 3D: Concatenate ligand and protein DataFrames
    combined_df = pd.concat([ligand_df, protein_df])

    # Step 3E: Normalize and scale the combined DataFrame
    scaler = MinMaxScaler()  # Use MinMaxScaler instead of StandardScaler
    scaled_data = scaler.fit_transform(combined_df)

    # Step 3F: Separate the combined DataFrame back into the ligand DataFrame and the protein DataFrame
    scaled_ligand_df = pd.DataFrame(scaled_data[:1], index=ligand_df.index, columns=common_columns)
    scaled_protein_df = pd.DataFrame(scaled_data[1:], index=protein_df.index, columns=common_columns)
    
    return scaled_ligand_df, scaled_protein_df

# Preprocess data
scaled_ligand_df, scaled_protein_df = preprocess_data(protein_data, ligand_code, ligand_smiles)

#Step 4: Euclidean Distance Calculations,
#PCA-derived Feature Similarity Mapping
#and Visualization

def compute_euclidean_distance_and_visualize(ligand_df, protein_df, visualize=True):
    """
    Computes the Euclidean distance between the input ligand and each protein in the protein feature space,
    and prompts the user if they would like to see a Euclidean distance-informed PCA plot with an accompanying similarity map.
    If visualize=True, a PCA plot and similarity map will be created, with the ligand shown in blue and proteins in red.
    Returns a tuple containing the most similar protein match and its similarity score.
    """
    # Compute Euclidean distances
    print("Processing scores (1/4)...")
    distances = np.linalg.norm(protein_df.values - ligand_df.values, axis=1)
    distance_df = pd.DataFrame(distances, index=protein_df.index, columns=['distance'])
    pca = PCA(n_components=2)
    pca.fit(pd.concat([ligand_df, protein_df]))

    # PCA and visualization
    if visualize:
        user_input = input("Would you like to see a Euclidean distance-informed PCA plot with a similarity map? (y/n): ")
        if user_input.lower() == "y":
            # Perform PCA
            combined_df = pd.concat([ligand_df, protein_df])
            pca_result = pca.transform(combined_df)  # Remove the redundant pca.fit_transform() line
            pca_embedding = pd.DataFrame(pca_result, index=combined_df.index)
            # Plot PCA
            fig, ax = plt.subplots()
            ax.scatter(pca_result[0, 0], pca_result[0, 1], c='blue', label=ligand_df.index[0])
            ax.scatter(pca_result[1:, 0], pca_result[1:, 1], c='red', label='Proteins')

            # Label points
            for i, txt in enumerate(combined_df.index):
                ax.annotate(txt, (pca_result[i, 0], pca_result[i, 1]))

            ax.set_title("Euclidean Distance-Informed PCA Plot")
            ax.legend()
            plt.show()

            # Plot similarity map
            fig, ax = plt.subplots()
            cmap = plt.cm.get_cmap("YlOrRd")
            cax = ax.matshow(distance_df.T, cmap=cmap)
            ax.set_xticks(range(len(distance_df.index)))
            ax.set_xticklabels(distance_df.index, rotation=90)
            ax.set_yticks([])
            ax.set_title("Euclidean Distance Similarity Map")
            fig.colorbar(cax)
            plt.show()

        else:
            print("No visualization will be shown.")
    else:
        print("Processing complete.")

    
    # Get most similar protein match and its similarity score
    most_similar_match = distance_df['distance'].idxmin()
    similarity_score = distance_df['distance'].min()

    pca_result = pca.transform(pd.concat([ligand_df, protein_df]))
    pca_embedding = pd.DataFrame(pca_result, index=pd.concat([ligand_df, protein_df]).index)

    return pca_embedding, most_similar_match, similarity_score

# Compute Euclidean distance and visualize
pca_embedding, most_similar_match, similarity_score = compute_euclidean_distance_and_visualize(scaled_ligand_df, scaled_protein_df)
#Step 5: Optimal Transport UMAP Visualization
import umap
from scipy.spatial.distance import cdist
from scipy.stats import wasserstein_distance
import matplotlib.pyplot as plt

def compute_optimal_transport_UMAP_and_visualize(scaled_ligand_df, scaled_protein_df):
    """
    This function computes a Wasserstein-informed UMAP for the selected ligand and proteins.
    It also creates a similarity map using the Wasserstein distances.
    The ligand and protein points are labeled in the UMAP space.
    """

    # Step 5A: Combine the ligand and protein dataframes
    combined_df = pd.concat([scaled_ligand_df, scaled_protein_df])

    # Step 5B: Compute UMAP
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(combined_df)
    umap_embedding = pd.DataFrame(embedding, index=combined_df.index)
    

    # Step 5C: Compute Wasserstein distances
    wasserstein_distances = cdist(scaled_ligand_df, scaled_protein_df, metric=wasserstein_distance)

    # Step 5D: Visualize UMAP
    scatter_proteins = plt.scatter(embedding[1:, 0], embedding[1:, 1], c=wasserstein_distances, cmap='YlOrRd', label='Proteins')
    plt.scatter(embedding[0, 0], embedding[0, 1], c='blue', label='Selected Ligand')

    for i, label in enumerate(combined_df.index):
        plt.annotate(label, (embedding[i, 0], embedding[i, 1]))

    plt.colorbar(scatter_proteins, label='Wasserstein Distance')
    plt.title('Wasserstein-Informed UMAP and Similarity Map')
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    plt.legend()
    plt.show()

    return umap_embedding, wasserstein_distances

# Prompt user if they want to visualize UMAP and similarity map
view_umap = input("Would you like to visualize the Wasserstein-informed UMAP and similarity map? (y/n) ").lower()
if view_umap == 'y':
    print("Processing scores (2/4)...")
    compute_optimal_transport_UMAP_and_visualize(scaled_ligand_df, scaled_protein_df)
else:
    print("No visualization will be shown.")
    print("Processing scores (2/4)...")


#Step 6: Compute Isomap, tSNE, and visualize
from sklearn.manifold import Isomap, TSNE

def compute_isomap_embedding(data, scaled_ligand_df, n_neighbors=5, n_components=2):
    isomap = Isomap(n_neighbors=n_neighbors, n_components=n_components)
    combined_data = pd.concat([scaled_ligand_df, data])
    embedding_array = isomap.fit_transform(combined_data)
    
    # Convert the NumPy array to a DataFrame
    embedding_df = pd.DataFrame(embedding_array, index=combined_data.index)
    
    return embedding_df

def compute_tsne_embedding(data, scaled_ligand_df, n_components=2, perplexity=5):
    tsne = TSNE(n_components=n_components, perplexity=perplexity)
    combined_data = pd.concat([scaled_ligand_df, data])
    embedding_array = tsne.fit_transform(combined_data)
    
    # Convert the NumPy array to a DataFrame
    embedding_df = pd.DataFrame(embedding_array, index=combined_data.index)
    
    return embedding_df


print("Computing Isomap and tSNE embeddings...")

import matplotlib.pyplot as plt

def plot_embedding(embedding, labels, ligand_code, title):
    fig, ax = plt.subplots(figsize=(8, 6))

    for i, label in enumerate(labels):
        if label == ligand_code:
            marker = 'o'
            color = 'blue'
        else:
            marker = 'x'
            color = 'red'

        ax.scatter(embedding[i, 0], embedding[i, 1], c=color, marker=marker, label=label)
        ax.text(embedding[i, 0], embedding[i, 1], label)

    ax.set_title(title)
    plt.show()

def visualize_isomap(data, ligand_code, protein_labels, scaled_ligand_df, n_neighbors=5, n_components=2):
    embedding = compute_isomap_embedding(data, scaled_ligand_df, n_neighbors=n_neighbors, n_components=n_components)
    
    fig, ax = plt.subplots()
    ax.scatter(embedding[0, 0], embedding[0, 1], c='blue', label=ligand_code, marker='o')
    ax.scatter(embedding[1:, 0], embedding[1:, 1], c='red', label='Proteins', marker='o')
    
    for i, txt in enumerate(protein_labels):
        ax.annotate(txt, (embedding[i + 1, 0], embedding[i + 1, 1]))
    
    ax.legend()
    plt.title('Isomap Embedding of Select Ligand and Proteins')
    plt.show()

    combined_df = pd.concat([scaled_ligand_df, scaled_protein_df])
    isomap_embedding = pd.DataFrame(embedding, index=combined_df.index)
    return isomap_embedding, most_similar_match, similarity_score


def visualize_tsne(data, ligand_code, protein_labels, scaled_ligand_df, n_components=2, perplexity=5):
    embedding = compute_tsne_embedding(data, scaled_ligand_df, n_components=n_components, perplexity=perplexity)
    
    fig, ax = plt.subplots()
    ax.scatter(embedding[0, 0], embedding[0, 1], c='blue', label=ligand_code, marker='o')
    ax.scatter(embedding[1:, 0], embedding[1:, 1], c='red', label='Proteins', marker='o')
    
    for i, txt in enumerate(protein_labels):
        ax.annotate(txt, (embedding[i + 1, 0], embedding[i + 1, 1]))
    
    ax.legend()
    plt.title('t-SNE Embedding of Select Ligand and Proteins')
    plt.show()

    combined_df = pd.concat([scaled_ligand_df, scaled_protein_df])
    tsne_embedding = pd.DataFrame(embedding, index=combined_df.index)
    return tsne_embedding, most_similar_match, similarity_score



from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import MinMaxScaler

def create_similarity_map(embedding_tuple):
    embedding_df = embedding_tuple[0]
    receptor = embedding_tuple[1]
    similarity_value = embedding_tuple[2]

    distances = pairwise_distances(embedding_df)
    similarity_map = 1 / (1 + distances)
    return similarity_map


# Plot heatmap function
def plot_1d_heatmap(similarity_map, protein_labels, ligand_code, title):
    fig, ax = plt.subplots(figsize=(10, 2))
    # Use iloc to access the DataFrame values
    im = ax.imshow(similarity_map.iloc[0, 1:].values.reshape(1, -1), cmap="YlOrRd_r", aspect='auto')

    ax.set_xticks(np.arange(len(protein_labels)))
    ax.set_yticks([])

    ax.set_xticklabels(protein_labels)

    plt.xticks(rotation=90)
    plt.title(title)
    plt.colorbar(im, ax=ax)
    plt.show()


# Visualize Isomap embedding and similarity map
isomap_embedding = compute_isomap_embedding(scaled_protein_df, scaled_ligand_df)
tsne_embedding = compute_tsne_embedding(scaled_protein_df, scaled_ligand_df)
plot_1d_heatmap(isomap_embedding, pdb_ids, ligand_code, 'Isomap 1D Similarity Map')
plot_1d_heatmap(tsne_embedding, pdb_ids, ligand_code, 't-SNE 1D Similarity Map')

#Step 7: Ligand-Fang Weighted Meta-Scorer Application
import numpy as np
from sklearn.preprocessing import MinMaxScaler
print("Processing scores (3/4)...")
# Assign weights to each dimensionality reduction technique
weights = {
    "pca": 0.4,
    "umap": 0.3,
    "isomap": 0.2,
    "tsne": 0.1
}
# Compute any finnicky embeddings without visualization
umap_embedding, wasserstein_distances = compute_optimal_transport_UMAP_and_visualize(scaled_ligand_df, scaled_protein_df)

# Compute the weighted composite scores
weighted_composite_scores = (pca_embedding.iloc[1:] * weights['pca'] +
                             umap_embedding.iloc[1:] * weights['umap'] +
                             isomap_embedding.iloc[1:] * weights['isomap'] +
                             tsne_embedding.iloc[1:] * weights['tsne'])


# Set the index to pdb_ids
weighted_composite_scores.index = pdb_ids
# Find the indices of the three most similar protein matches
most_similar_indices = weighted_composite_scores.sum(axis=1).nsmallest(3).index

# Retrieve the most similar protein matches and their similarity scores
top_matches = most_similar_indices.values
top_scores = weighted_composite_scores.sum(axis=1).nsmallest(3).values

# Print in a robotic way
print("Finalizing... (4/4)")
print()

# Print as a loop which says for the first 3 top matches, print the match and score
for i in range(3):
    print("Ligand {} Match #{}: {}, Score: {:.4f}".format(select, i + 1, top_matches[i], top_scores[i]))
    print()
print("Done!")