'''
# Create a SMARTS mapping (just an example, you can expand this)
smarts_mapping = {
    "H": "[H]",
    "O": "[O]",
    "C": "[C]",
}

# Create a PubChem mapping (just an example, you can expand this)
pubchem_mapping = {
    "H": "1",
    "O": "8",
    "C": "6",
}

# Initialize the database
cat_db = CATDatabase(filepath="chemical_data.xlsx")

# Fetch periodic table data
cat_db.get_periodic_table_data()

# Fetch DeMaLi properties
cat_db.get_deml_properties()

# Fetch reaction SMARTS
cat_db.get_reaction_smarts(smarts_mapping)

# Fetch PubChem notations
cat_db.get_pubchem_notation(pubchem_mapping)

# Show the modified DataFrame
print(cat_db.df)
'''


from matminer.featurizers.composition import ElementProperty
from pymatgen.core.composition import Composition
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

class CATDatabase:
    """
    Class for handling chemical reactions and elements data from an Excel file.
    """
    def __init__(self, filepath):
        """
        Initialize CATDatabase with an Excel file containing chemical data.

        Parameters
        ----------
        filepath : str
            Path to the Excel file containing chemical data.
        """

        self.df = pd.read_excel(filepath)
        self.ptable = Chem.GetPeriodicTable()

    def get_periodic_table_data(self):
        """
        Fetches periodic table data for surface symbols in the dataframe.
        """

        surface_symbols = self.df["Surface"]
        self.df["Surf_Element"], self.df["Surf_Atom_N"], self.df["Surf_Atom_M"], self.df["Surf_Elec_V"], self.df["Surf_VdW_R"], self.df["Surf_Cov_R"], self.df["Surf_Bond_R"], self.df["Surf_Common_Iso"], self.df["Surf_Common_Iso_M"] = zip(
            *surface_symbols.map(self._fetch_element_data))

    def _fetch_element_data(self, symbol):
        """
        Fetches periodic table data for a given element symbol.

        Parameters
        ----------
        symbol : str
            Symbol of the element.

        Returns
        -------
        tuple
            Returns a tuple containing various properties of the element.
        """

        atomic_number = self.ptable.GetAtomicNumber(symbol)
        element_name = self.ptable.GetElementName(atomic_number)
        atomic_mass = self.ptable.GetAtomicWeight(symbol)
        vlc_elec = self.ptable.GetNOuterElecs(symbol)
        vdw_r = self.ptable.GetRvdw(symbol)
        cov_r = self.ptable.GetRcovalent(symbol)
        bond_r = self.ptable.GetRb0(symbol)
        common_isotope = self.ptable.GetMostCommonIsotope(symbol)
        common_isotope_m = self.ptable.GetMostCommonIsotopeMass(symbol)
        return element_name, atomic_number, atomic_mass, vlc_elec, vdw_r, cov_r, bond_r, common_isotope, common_isotope_m

    def get_deml_properties(self):
        """
        Fetches DeMaLi (Density, Melting Point, and Lattice Constants) properties for surface symbols.
        """

        featurizer = ElementProperty.from_preset("deml")
        self.df["Surface"] = self.df["Surface"].apply(Composition)
        self.df = featurizer.featurize_dataframe(self.df, col_id="Surface")

    def get_reaction_smarts(self, smarts_mapping):
        """
        Maps reaction symbols to their respective SMARTS notations.

        Parameters
        ----------
        smarts_mapping : dict
            Dictionary mapping elements to their SMARTS notations.
        """

        reactant_smarts_list, product_smarts_list = zip(
            *self.df["Reaction"].map(lambda x: self._react_from_smarts(x, smarts_mapping)))
        self.df['Reactants_SMARTS'] = reactant_smarts_list
        self.df['Products_SMARTS'] = product_smarts_list

    def _react_from_smarts(self, reaction_str, smarts_mapping):
        """
        Converts a reaction string to its SMARTS notation.

        Parameters
        ----------
        reaction_str : str
            String describing the reaction.
        smarts_mapping : dict
            Dictionary mapping elements to their SMARTS notations.

        Returns
        -------
        tuple
            Returns a tuple of SMARTS notations for reactants and products.
        """

        reactant_str, product_str = reaction_str.split("→")
        reactant_atoms = [atom.strip() for atom in reactant_str.split("+")]
        product_atoms = [atom.strip() for atom in product_str.split("+")]
        reactant_smarts = [smarts_mapping.get(atom, atom) for atom in reactant_atoms]
        product_smarts = [smarts_mapping.get(atom, atom) for atom in product_atoms]
        return reactant_smarts, product_smarts

    def get_pubchem_notation(self, pubchem_mapping):
        """
        Maps reaction symbols to their respective PubChem notations.

        Parameters
        ----------
        pubchem_mapping : dict
            Dictionary mapping elements to their PubChem notations.
        """

        self.df["PubChem_Notation"] = self.df["Reaction"].map(
            lambda x: self._get_pubchem_notation(x, pubchem_mapping))

    def _get_pubchem_notation(self, reaction_str, pubchem_mapping):
        """
        Converts a reaction string to its PubChem notation.

        Parameters
        ----------
        reaction_str : str
            String describing the reaction.
        pubchem_mapping : dict
            Dictionary mapping elements to their PubChem notations.

        Returns
        -------
        tuple
            Returns a tuple of PubChem notations for reactants and products.
        """

        reactant_str, product_str = reaction_str.split("→")
        reactant_atoms = [atom.strip() for atom in reactant_str.split("+")]
        product_atoms = [atom.strip() for atom in product_str.split("+")]
        reactant_notation = [pubchem_mapping.get(atom, atom) for atom in reactant_atoms]
        product_notation = [pubchem_mapping.get(atom, atom) for atom in product_atoms]
        return reactant_notation, product_notation
