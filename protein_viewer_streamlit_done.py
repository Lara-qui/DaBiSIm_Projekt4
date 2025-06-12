# 1. Wichtige Module importieren
import os, tempfile                           # für Dateizugriff und temporäre Dateien
import numpy as np                            # numerische Berechnungen (Koordinaten)
import streamlit as st                        # GUI in Webbrowser
from Bio.PDB import PDBParser                 # Biopython zum Parsen von PDB-Dateien
from periodictable import elements            # enthält Atommasse je Element
import plotly.graph_objects as go             # für 3D-Visualisierung
from collections import defaultdict           # für Gruppierung von Atomen

# 2. Datenklassen definieren
class Atom:
    def __init__(self, element, coord):
        self.element = element.strip().capitalize()  # z. B. "C", "O", "N"
        self.coord = coord                           # 3D-Koordinaten als NumPy-Array

    def get_mass(self):  # Gibt Atommasse aus dem Periodensystem zurück
        try:
            return getattr(elements, self.element).mass
        except AttributeError:
            return 0.0

class Residue:
    def __init__(self, name):
        self.name = name          # z. B. ALA, GLY
        self.atoms = []           # Liste von Atom-Objekten

    def add_atom(self, atom):
        self.atoms.append(atom)

class Chain:
    def __init__(self, chain_id):
        self.chain_id = chain_id  # z. B. "A", "B"
        self.residues = []        # Liste von Residue-Objekten

    def add_residue(self, residue):
        self.residues.append(residue)

class Protein:
    def __init__(self, name):
        self.name = name
        self.chains = []          # Liste aller Ketten

    def add_chain(self, chain):
        self.chains.append(chain)

    def get_all_atoms(self):
        # Gibt alle Atome aus allen Residuen aller Ketten zurück
        return [atom for chain in self.chains for res in chain.residues for atom in res.atoms]

    def calculate_molecular_weight(self):
        # Berechnet das Gesamtgewicht des Proteins aus der Atommasse
        return sum(atom.get_mass() for atom in self.get_all_atoms())

# 3. PDB-Datei einlesen und Proteinstruktur erzeugen
def extract_title_from_pdb(file_path):
    # Liest den TITLE-Eintrag aus der PDB-Datei (falls vorhanden)
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("TITLE"):
                return line[10:].strip().capitalize()
    return "Unbekannter Proteinname"

def parse_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    protein = Protein(os.path.basename(pdb_file))

    for model in structure:  # Modell (z. B. bei NMR mehrere Modelle möglich)
        for chain in model:  # z. B. A, B, ...
            c = Chain(chain.id)
            for residue in chain:
                r = Residue(residue.resname)
                for atom in residue:
                    a = Atom(atom.element, atom.coord)
                    r.add_atom(a)
                c.add_residue(r)
            protein.add_chain(c)
    return protein

# 4. 3D-Visualisierung mit Plotly
def visualize_protein_3d(protein, selected_chains, selected_elements, selected_residues, color_mode, protein_title):
    atoms_by_group = defaultdict(list)  # Gruppierung je nach Farbschema

    # Farbpalette für Residuen (20 Standardfarben)
    color_palette = [
        'red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta', 'lime', 'yellow', 'pink',
        'gold', 'brown', 'gray', 'olive', 'teal', 'navy', 'maroon', 'aqua', 'indigo', 'violet'
    ]
    residue_colors = {}  # Zuordnung Residuum → Farbe
    color_index = 0

    for chain in protein.chains:
        if selected_chains and chain.chain_id not in selected_chains:
            continue
        for residue in chain.residues:
            if selected_residues and residue.name not in selected_residues:
                continue
            # Neue Farbe für bisher nicht vergebene Residuen
            if color_mode == "Nach Residuum" and residue.name not in residue_colors:
                residue_colors[residue.name] = color_palette[color_index % len(color_palette)]
                color_index += 1
            for atom in residue.atoms:
                if selected_elements and atom.element not in selected_elements:
                    continue

                # Gruppierungsschlüssel je nach Modus
                if color_mode == "Nach Kette":
                    key = chain.chain_id
                elif color_mode == "Nach Element":
                    key = atom.element
                elif color_mode == "Nach Residuum":
                    key = residue.name

                atoms_by_group[key].append(atom.coord)

    # Farben für bekannte chemische Elemente
    element_colors = {
        'H': 'white', 'C': 'green', 'N': 'blue', 'O': 'red', 'S': 'yellow',
        'P': 'orange', 'Fe': 'darkgray', 'Zn': 'purple', 'Cl': 'green',
        'Na': 'deepskyblue', 'K': 'violet', 'Ca': 'limegreen', 'Mg': 'teal'
    }

    # Farben je Kette (z. B. A, B, C, ...)
    chain_colors = {
        cid: col for cid, col in zip(
            sorted(set(chain.chain_id for chain in protein.chains)),
            ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        )
    }

    fig = go.Figure()

    # Visualisierung aller Gruppen
    for key, coords in atoms_by_group.items():
        coords = np.array(coords)
        if len(coords) == 0:
            continue

        # Farbauswahl je Modus
        if color_mode == "Nach Kette":
            color = chain_colors.get(key, 'gray')
            name = f"Kette {key}"
        elif color_mode == "Nach Element":
            color = element_colors.get(key, 'black')
            name = f"Element {key}"
        elif color_mode == "Nach Residuum":
            color = residue_colors.get(key, 'gold')
            name = f"Residuum {key}"

        fig.add_trace(go.Scatter3d(
            x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
            mode='markers',
            marker=dict(size=4, color=color),
            name=name
        ))

    # Plot-Layout mit Achsenbeschriftung
    fig.update_layout(
        title=f'3D Struktur: {protein_title}',
        width=800, height=700,
        scene=dict(
            xaxis=dict(title='X', showspikes=False, showaxeslabels=True),
            yaxis=dict(title='Y', showspikes=False, showaxeslabels=True),
            zaxis=dict(title='Z', showspikes=False, showaxeslabels=True)
        ),
        margin=dict(l=0, r=0, b=0, t=40),
        legend=dict(itemsizing='constant'),
    )

    # Reset-Kamera
    if st.button("Ansicht zurücksetzen"):
        fig.update_layout(scene_camera=dict(eye=dict(x=1.5, y=1.5, z=1.5)))

    # HTML-Export
    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_html:
        fig.write_html(tmp_html.name)
        st.download_button(
            label="3D-Visualisierung als HTML speichern",
            data=open(tmp_html.name, 'rb').read(),
            file_name='protein_visualisierung.html',
            mime='text/html'
        )

    # Plot anzeigen
    st.plotly_chart(fig, use_container_width=True)

# 5. Streamlit-Benutzeroberfläche
def run_gui():
    st.title("🧬 Proteinstruktur-Analyse")

    uploaded_file = st.file_uploader("Wähle eine PDB-Datei", type="pdb")
    if uploaded_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
            tmp.write(uploaded_file.read())
            tmp_path = tmp.name

        protein = parse_pdb(tmp_path)
        title = extract_title_from_pdb(tmp_path)
        weight = protein.calculate_molecular_weight()

        # Proteininfos anzeigen
        st.info(f"**Proteinname:** {title}")
        st.success(f"Molekulargewicht: {weight:.2f} Da")

        # Auswahloptionen
        all_chains = sorted(set(chain.chain_id for chain in protein.chains))
        all_elements = sorted(set(atom.element for atom in protein.get_all_atoms()))
        all_residues = sorted(set(res.name for chain in protein.chains for res in chain.residues))

        selected_chains = st.multiselect("Ketten auswählen:", all_chains, default=all_chains)
        selected_elements = st.multiselect("Elemente auswählen:", all_elements, default=all_elements)
        selected_residues = st.multiselect("Residuen auswählen:", all_residues, default=all_residues)

        # Farbschema wählen
        color_mode = st.radio("Farbmodus", ["Nach Kette", "Nach Element", "Nach Residuum"], horizontal=True)

        # Darstellung
        st.write("### 3D-Struktur")
        visualize_protein_3d(protein, selected_chains, selected_elements, selected_residues, color_mode, title)

# 6. App starten
if __name__ == '__main__':
    run_gui()
