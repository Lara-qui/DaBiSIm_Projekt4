#1. Import und Setup
from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from periodictable import elements
import streamlit as st
import tempfile
import os

#2. Definition von OOP-Klassenstruktur
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

#3. Parser Funktion
def extract_title_from_pdb(file_path):             # Proteinname extrahieren
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("TITLE"):
                return line[10:].strip()
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

# 4. Visualisierung der Proteinstruktur
def visualize_protein_3d(protein):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    xs, ys, zs = [], [], []
    for atom in protein.get_all_atoms():
        x, y, z = atom.coord
        xs.append(x)
        ys.append(y)
        zs.append(z)

    ax.scatter(xs, ys, zs, s=20, c='skyblue', edgecolor='k')
    ax.set_title(f"3D-Struktur von {protein.name}")
    st.pyplot(fig)

#5. Streamlit GUI
def run_gui():
    st.title("Proteinstruktur-Analyse")
    st.write("Lade eine PDB-Datei hoch:")

    uploaded_file = st.file_uploader("Wähle eine PDB-Datei", type="pdb")

    if uploaded_file is not None:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
            tmp.write(uploaded_file.read())
            tmp_path = tmp.name

        protein = parse_pdb(tmp_path)
        title = extract_title_from_pdb(tmp_path)
        st.info(f"**Proteinname:** {title}")
        weight = protein.calculate_molecular_weight()
        st.success(f"Molekulargewicht: {weight:.2f} Da")
        st.write("### 3D-Struktur")
        visualize_protein_3d(protein)

if __name__ == '__main__':
    run_gui()
#plotly einfügen, dropdown menü?
#verschiedene Elemente und Ketten darstellen?
