#1. Import und Setup der nötigen Libraries
from Bio.PDB import PDBParser                 # für den Import einer PDB-Datei
import matplotlib.pyplot as plt               # ermöglicht Plotting
from mpl_toolkits.mplot3d import Axes3D       # für 3D-plotting
from periodictable import elements            # Information über die Atommasse von Elementen
import streamlit as st                        # bildet die GUI Umgebung
import tempfile                               # ermöglicht den Umgang mit temporären Dateien
import os                                     # ermöglicht den Zugriff auf die Datei

#2. Definition von OOP-Klassenstruktur
class Atom:                                   # ordnet dem Begriff "Atom" Koordinaten und das Element zu
    def __init__(self, element, coord):
        self.element = element
        self.coord = coord

    def get_mass(self):                       # ordnet dem Begriff "Atom" die Atommasse abhängig vom Element zu
        try:
            return getattr(elements, self.element).mass
        except AttributeError:
            return 0.0

class Residue:                                # dem Residuum wird ein Name (Aminosäure) und die zugehörigen Atome zugeordnet
    def __init__(self, name):
        self.name = name
        self.atoms = []

    def add_atom(self, atom):
        self.atoms.append(atom)

class Chain:                                  # alpha-/beta-Ketten werden den zugehörigen Residuen zugeordnet
    def __init__(self, chain_id):
        self.chain_id = chain_id
        self.residues = []

    def add_residue(self, residue):
        self.residues.append(residue)

class Protein:                               # dem Protein werden die Ketten zugeordnet, alle Atome einzeln zugeordnet und die Summe der Atommassen aller zugeordneten Atome berechnet
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
def parse_pdb(pdb_file):                     # Daten aus PDB-Datei werden gelesen, QUIET=True unterdrückt Warnungen in Fehlern der Datei, Dateiname wird dem Protein zugeordnet, Protein-Objekt wird erstellt & Ketten, Residuen, Elemente & Atomkoordinaten integriert (+Elementnamen standardisiert)
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
def visualize_protein_3d(protein):            # 3D-Scatterplot des Proteins wird erstellt
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
        weight = protein.calculate_molecular_weight()
        st.success(f"Molekulargewicht: {weight:.2f} Da")
        st.write("### 3D-Struktur")
        visualize_protein_3d(protein)

if __name__ == '__main__':
    run_gui()
#plotly einfügen, dropdown menü?
#verschiedene Elemente und Ketten darstellen?
