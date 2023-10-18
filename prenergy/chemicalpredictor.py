'''
# Import the Mol2VecModel class
# (Assuming the class definition is in a file called `mol2vec_model.py`)
from prenergy import chemicalpredictor.Mol2VecModel

# Initialize the Mol2VecModel with a data file and a model file
data_file = 'data.csv'  # Replace with the actual path to your data CSV file
model_file = 'model.w2v'  # Replace with the actual path to your Word2Vec model file
mol2vec = Mol2VecModel(data_file, model_file)

# Generate molecular vectors
mol2vec.generate_vectors()

# Fit the model and evaluate its performance
trained_model = mol2vec.fit_and_evaluate()

# Now, `trained_model` contains the trained Gradient Boosting Regressor model
# You can use `trained_model` for further predictions or analysis
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from mol2vec.features import mol2alt_sentence, sentences2vec, DfVec
from rdkit import Chem
from gensim.models import word2vec
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, explained_variance_score, max_error, median_absolute_error
import shap


class Mol2VecModel:
    """
    A class to perform molecular vectorization, dimensionality reduction, and energy prediction.
    
    Attributes
    ----------
    df : pandas.DataFrame
        The data frame containing molecular data and energy values.
    model : gensim.models.Word2Vec
        The pretrained Word2Vec model for generating molecular vectors.
    df_array : np.ndarray
        Array of molecules used for vectorization.
    data_final : list
        List to store final average vectors for molecules.
    indices : list
        List to store indices of successfully processed molecules.
    """
    def __init__(self, data_file, model_file):
        """
        Initialize Mol2VecModel with data file and model file.

        Parameters
        ----------
        data_file : str
            The CSV file path containing molecular data.
        model_file : str
            The file path for the pretrained Word2Vec model.
        """

        self.df = pd.read_csv(data_file)
        self.model = word2vec.Word2Vec.load(model_file)
        self._prepare_data()
        
    def _prepare_data(self):
        """
        Prepare the initial data for molecular vectorization.
        Filters and cleans the dataframe.
        """

        # Perform initial data cleaning and filtering
        self.df = self.df[['Energy'] + [col for col in self.df.columns if col != 'Energy']]
        self.df.dropna(inplace=True)
        self.df_array = self.df[['Reactants_SMARTS', 'Products_SMARTS']].values
        
    def generate_vectors(self):
        """
        Generate molecular vectors by converting SMARTS to average vectors.
        """

        # Logic to convert SMARTS to average vectors
        self.data_final = []
        self.indices = []
        
        for indice, val in enumerate(self.df_array):
            try:
                val = [string.replace('[*]', '') for string in val]
                reactant_smarts, product_smarts = self._process_smarts(val)
                avg_vectors = self._convert_to_vectors(reactant_smarts, product_smarts)
                self.data_final.append(avg_vectors)
                self.indices.append(indice)
            except Exception as e:
                print(f"Error at index {indice}: {e}")
                
    def _process_smarts(self, val):
        """
        Process SMARTS data and split them into reactant and product SMARTS.

        Parameters
        ----------
        val : list
            List containing SMARTS for reactants and products.

        Returns
        -------
        tuple
            Tuple of lists containing reactant and product SMARTS.
        """

        reactant_smarts = val[0][2:-2].split(',')
        product_smarts = val[1][2:-2].split(',')
        reactant_smarts = [x.replace("'", "").strip() for x in reactant_smarts if x != '[*]']
        product_smarts = [x.replace("'", "").strip() for x in product_smarts if x != '[*]']
        return reactant_smarts, product_smarts

    def _convert_to_vectors(self, reactant_smarts, product_smarts):
        """
        Convert the reactant and product SMARTS to molecular vectors.

        Parameters
        ----------
        reactant_smarts : list
            List of SMARTS for reactants.
        product_smarts : list
            List of SMARTS for products.

        Returns
        -------
        np.ndarray
            The average vector for reactant and product molecules.
        """

        reactant_sentences = [mol2alt_sentence(Chem.MolFromSmarts(smart), 1) for smart in reactant_smarts]
        product_sentences = [mol2alt_sentence(Chem.MolFromSmarts(smart), 1) for smart in product_smarts]
        reactant_vectors = [DfVec(x).vec for x in sentences2vec(reactant_sentences, self.model, unseen='UNK')]
        product_vectors = [DfVec(x).vec for x in sentences2vec(product_sentences, self.model, unseen='UNK')]
        avg_vectors = ((np.mean(reactant_vectors, axis=0) + np.mean(product_vectors, axis=0)) / 2).flatten()
        return avg_vectors

    def fit_and_evaluate(self):
        """
        Perform dimensionality reduction, train a predictive model, and evaluate its performance.

        Returns
        -------
        sklearn.ensemble.GradientBoostingRegressor
            The trained Gradient Boosting Regressor model.
        """

        # Perform dimensionality reduction and machine learning
        pca_model = PCA(n_components=30)
        tsne_model = TSNE(n_components=2, perplexity=10, n_iter=1000, metric='cosine')
        tsne_pca = tsne_model.fit_transform(pca_model.fit_transform(self.data_final))
        
        # Filter the original DataFrame and add the 2D representation
        self.df = self.df.iloc[self.indices]
        self.df['tsne-2d-one'] = tsne_pca[:, 0]
        self.df['tsne-2d-two'] = tsne_pca[:, 1]
        self.df.drop(columns=['Reactants_SMARTS', 'Products_SMARTS'], inplace=True)
        
        # Data scaling and splitting
        scaler = StandardScaler()
        X = self.df.drop('Energy', axis=1)
        y = self.df['Energy']
        X = scaler.fit_transform(X)
        X_train, X_test, y_train, y_test = train_test_split(X, y)
        
        # Train and evaluate the model
        model = GradientBoostingRegressor()
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        
        # Use different metrics for evaluation
        metrics = {
            'MSE': mean_squared_error,
            'MAE': mean_absolute_error,
            'R^2': r2_score,
            'Explained Variance': explained_variance_score,
            'Max Error': max_error,
            'Median Absolute Error': median_absolute_error
        }

        print("Model: GradientBoosting")
        for metric_name, metric_func in metrics.items():
            metric_value = metric_func(y_test, y_pred)
            print(f"{metric_name}: {metric_value}")

        return model
