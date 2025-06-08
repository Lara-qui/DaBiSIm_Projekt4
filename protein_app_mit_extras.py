# 1. Imports und Setup
import os
import tempfile
import numpy as np
import streamlit as st
from Bio.PDB import PDBParser
from periodictable import elements
import plotly.graph_objects as go
from collections import defaultdict

# 2. Datenklassen für Proteinstruktur
class Atom:
    def __init__(self, element, coord):
        self.element = element
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

# 3. PDB Parsen und Titel extrahieren
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
                    a = Atom(element, atom.coord)
                    r.add_atom(a)
                c.add_residue(r)
            protein.add_chain(c)
    return protein

# 4. Protein in 3D visualisieren mit verschiedenen Farbmodi
def visualize_protein_3d(protein, selected_chains, selected_elements, color_mode):
    atoms_by_group = defaultdict(list)

    for chain in protein.chains:
        if selected_chains and chain.chain_id not in selected_chains:
            continue
        for residue in chain.residues:
            for atom in residue.atoms:
                if selected_elements and atom.element not in selected_elements:
                    continue
                if color_mode == "Nach Kette":
                    key = chain.chain_id
                elif color_mode == "Nach Element":
                    key = atom.element
                else:
                    key = (chain.chain_id, atom.element)
                atoms_by_group[key].append(atom.coord)

    element_colors = {
        'H': 'white', 'C': 'gray', 'N': 'blue', 'O': 'red', 'S': 'yellow',
        'P': 'orange', 'Fe': 'darkorange', 'Zn': 'purple', 'Cl': 'green',
        'Na': 'cyan', 'K': 'violet', 'Ca': 'lime', 'Mg': 'teal'
    }
    chain_colors = {
        cid: col for cid, col in zip(
            sorted(set(chain.chain_id for chain in protein.chains)),
            ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        )
    }

    fig = go.Figure()

    for key, coords in atoms_by_group.items():
        coords = np.array(coords)
        if len(coords) == 0:
            continue

        if color_mode == "Nach Kette":
            color = chain_colors.get(key, 'gray')
            name = f"Kette {key}"
            line_color = color
        elif color_mode == "Nach Element":
            color = element_colors.get(key, 'black')
            name = f"Element {key}"
            line_color = color
        else:
            chain_id, element = key
            color = element_colors.get(element, 'black')
            line_color = chain_colors.get(chain_id, 'gray')
            name = f"{chain_id} - {element}"

        fig.add_trace(go.Scatter3d(
            x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
            mode='markers',
            marker=dict(size=4, color=color, line=dict(width=2 if color_mode == "Kombiniert (Kette + Element)" else 0, color=line_color)),
            name=name
        ))

    fig.update_layout(
        title='3D Proteinstruktur', width=800, height=700,
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
        legend=dict(itemsizing='constant')
    )

    st.plotly_chart(fig)

# 5. Streamlit GUI
st.title("Proteinstruktur-Analyse")
uploaded_file = st.file_uploader("Wähle eine PDB-Datei", type="pdb")

if uploaded_file:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
        tmp.write(uploaded_file.read())
        tmp_path = tmp.name

    protein = parse_pdb(tmp_path)
    title = extract_title_from_pdb(tmp_path)
    weight = protein.calculate_molecular_weight()

    st.info(f"**Proteinname:** {title}")
    st.success(f"Molekulargewicht: {weight:.2f} Da")

    all_chains = sorted(set(chain.chain_id for chain in protein.chains))
    all_elements = sorted(set(atom.element for atom in protein.get_all_atoms()))

    if len(all_chains) > 1:
        selected_chains = st.multiselect("Wähle Ketten aus:", options=all_chains, default=all_chains)
    else:
        selected_chains = all_chains

    selected_elements = st.multiselect("Wähle Elemente aus:", options=all_elements, default=all_elements)

    color_mode = st.radio("Farbmodus wählen", ["Nach Kette", "Nach Element", "Kombiniert (Kette + Element)"], horizontal=True)

    st.write("### 3D-Struktur")
    visualize_protein_3d(protein, selected_chains, selected_elements, color_mode)


if __name__ == '__main__':
    run_gui()
