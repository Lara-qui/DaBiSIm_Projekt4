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
                key = atom
