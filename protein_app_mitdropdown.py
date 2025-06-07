# Import und Setup
from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from periodictable import elements
import streamlit as st
import tempfile
import os

# Definition von OOP-Klassenstruktur
class Atom:
    def __init__(self, element, coord):
        self.element = element.capitalize()
        self.coord = coord

    def get_mass(self):
        try:
            return getattr(elements, self.element).mass
        except AttributeError:
            return 0.0

class Residue:
    def __init__(self, name):
        self.name = name
        self.atoms = []

    def add_atom(self, atom):
        self.atoms.append(atom)

class Chain:
    def __init__(self, chain_id):
        self.chain_id = chain_id
        self.residues = []

    def add_residue(self, residue):
        self.residues.append(residue)

class Protein:
    def __init__(self, name):
        self.name = name
        self.chains = []

    def add_chain(self, chain):
        self.chains.append(chain)

    def get_all_atoms(self):
        return [atom for chain in self.chains for res in chain.residues for atom in res.atoms]

    def calculate_molecular_weight(self):
        return sum(atom.get_mass() for atom in self.get_all_atoms())

# Parser Funktion

def extract_title_from_pdb(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("TITLE"):
                return line[10:].strip().capitalize()
    return "Unbekannter Proteinname"

def parse_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    protein = Protein(os.path.basename(pdb_file))

    for model in structure:
        for chain in model:
            c = Chain(chain.id)
            for residue in chain:
                r = Residue(residue.resname)
                for atom in residue:
                    element = atom.element.strip().capitalize()
                    coord = atom.coord
                    a = Atom(element, coord)
                    r.add_atom(a)
                c.add_residue(r)
            protein.add_chain(c)
    return protein

# Visualisierung der Proteinstruktur

def visualize_protein_3d(protein, selected_chains, selected_elements, color_mode):
    import plotly.graph_objects as go
    import numpy as np
    from collections import defaultdict

    atoms_by_group = defaultdict(list)

    for chain in protein.chains:
        if selected_chains and chain.chain_id not in selected_chains:
            continue
        for residue in chain.residues:
            for atom in residue.atoms:
                if selected_elements and atom.element not in selected_elements:
                    continue
                key = atom.element if color_mode == "Nach Element" else chain.chain_id
                atoms_by_group[(key, chain.chain_id, atom.element)].append(atom.coord)

    element_colors = {
        'H': 'white', 'C': 'gray', 'N': 'blue', 'O': 'red', 'S': 'yellow',
        'P': 'orange', 'Fe': 'darkorange', 'Zn': 'purple', 'Cl': 'green',
        'Na': 'cyan', 'K': 'violet', 'Ca': 'lime', 'Mg': 'teal'
    }

    chain_palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    chain_colors = {}
    fig = go.Figure()

    for i, ((group, chain_id, element), coords) in enumerate(atoms_by_group.items()):
        coords = np.array(coords)
        if color_mode == "Nach Element":
            color = element_colors.get(group, 'black')
        else:
            if group not in chain_colors:
                chain_colors[group] = chain_palette[len(chain_colors) % len(chain_palette)]
            color = chain_colors[group]

        fig.add_trace(go.Scatter3d(
            x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
            mode='markers',
            marker=dict(size=4, color=color),
            name=f'Kette {chain_id} - {element}'
        ))

    fig.update_layout(
        title='3D-Struktur: Ketten & Elemente',
        width=800,
        height=700,
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
        legend=dict(itemsizing='constant')
    )

    st.plotly_chart(fig)

# Streamlit GUI

def run_gui():
    st.title("Proteinstruktur-Analyse")
    st.write("Lade eine PDB-Datei hoch:")

    if st.button("Zurücksetzen"):
        for key in st.session_state.keys():
            del st.session_state[key]
        st.experimental_rerun()

    uploaded_file = st.file_uploader("Wähle eine PDB-Datei", type="pdb")

    if uploaded_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
            tmp.write(uploaded_file.read())
            tmp_path = tmp.name

        protein = parse_pdb(tmp_path)
        title = extract_title_from_pdb(tmp_path)
        st.info(f"**Proteinname laut TITLE-Eintrag:** {title}")

        weight = protein.calculate_molecular_weight()
        st.success(f"Molekulargewicht: {weight:.2f} Da")

        all_chains = sorted(set(chain.chain_id for chain in protein.chains))
        all_elements = sorted(set(atom.element for atom in protein.get_all_atoms()))

        if len(all_chains) > 1:
            if "selected_chains" not in st.session_state:
                st.session_state.selected_chains = all_chains
            selected_chains = st.multiselect("Wähle Ketten aus:", options=all_chains, default=st.session_state.selected_chains, key="selected_chains")
        else:
            selected_chains = all_chains

        selected_elements = st.multiselect("Wähle Elemente aus:", options=all_elements, default=all_elements, key="selected_elements")

        color_mode = st.radio("Farbmodus für Visualisierung:", ["Nach Kette", "Nach Element"], key="color_mode")

        st.write("### 3D-Struktur")
        visualize_protein_3d(protein, selected_chains, selected_elements, color_mode)

if __name__ == '__main__':
    run_gui()
