### Imports
import os
import requests
import tarfile
import zipfile
import shutil
import pandas as pd
from typing import Dict
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import scipy.sparse as sp
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks, argrelextrema
from scipy.sparse import issparse, csr_matrix
from diptest import diptest
import statsmodels.api as sm
from shapely.geometry import Point, Polygon, MultiPoint
from scipy.spatial import ConvexHull
import re
import gzip

### Get 10x data from URL
def download_and_extract_10x_data(data_path: str, folder_path: str, urls: Dict[str, str]) -> None:
    """
    Downloads files from the given URLs and extracts any .tar.gz archives.

    Args:
        data_path (str): Base directory path where data will be stored.
        folder_path (str): Specific subfolder path where files will be downloaded and extracted.
        urls (Dict[str, str]): A dictionary mapping filenames to their download URLs.

    Example:
        urls = {
            "data.h5": "https://example.com/data.h5",
            "archive1.tar.gz": "https://example.com/archive1.tar.gz",
            "archive2.tar.gz": "https://example.com/archive2.tar.gz"
        }
        download_and_extract_10x_data("./data/", "./data/my_folder/", urls)

    Notes:
        - The function creates `folder_path` if it doesn't exist.
        - Any files ending with `.tar.gz` will be automatically extracted after download.
        - HTTP status codes are checked to ensure successful downloads.
    """
    # Create the folder path if it doesn't exist
    try:
        os.makedirs(folder_path, exist_ok=True)
        print(f"Folder path '{folder_path}' created successfully.")
    except Exception as e:
        print(f"An error occurred while creating the folder path: {e}")

    # Download files from the URLs
    for filename, url in urls.items():
        response = requests.get(url)
        if response.status_code == 200:
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print(f"{filename} successfully downloaded")

            # If the file is a .tar.gz, extract it
            if filename.endswith(".tar.gz"):
                try:
                    with tarfile.open(file_path, "r:gz") as tar:
                        tar.extractall(path=folder_path)
                    print(f"Extracted {filename} to: {folder_path}")
                except Exception as e:
                    print(f"An error occurred while extracting {filename}: {e}")
            
            # If the file is a .zip, extract it
            elif filename.endswith(".zip"):
                try:
                    with zipfile.ZipFile(file_path, 'r') as zip_ref:
                        zip_ref.extractall(folder_path)
                    print(f"Extracted {filename} to: {folder_path}")
                except Exception as e:
                    print(f"An error occurred while extracting {filename}: {e}")

        else:
            print(f"Failed to download {filename}, status code: {response.status_code}")

    # List final contents
    print("Final folder contents:", os.listdir(folder_path))


### Decompress .csv.gz File
def decompress_file(gz_path: str) -> None:
    """
    Decompresses a .gz file to its original CSV form.

    Args:
        gz_path (str): Path to the .gz compressed file.
    """
    csv_path = gz_path.rstrip('.gz')
    with gzip.open(gz_path, 'rb') as f_in:
        with open(csv_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Decompressed {gz_path} to {csv_path}")
    
### Load and Inspect 10x Data (Generic)
def load_and_inspect_10x_csv_data(sample_path: str) -> Dict[str, pd.DataFrame]:
    """
    Loads and inspects all CSV/CSV.GZ files in a given 10x data sample directory.
    Automatically decompresses any .csv.gz files if needed.

    Args:
        sample_path (str): Path to the unzipped 10x data sample folder.

    Returns:
        Dict[str, pd.DataFrame]: Dictionary where keys are file names (without extension)
                                 and values are corresponding pandas DataFrames.
    """
    data_frames = {}

    for file_name in os.listdir(sample_path):
        file_path = os.path.join(sample_path, file_name)

        if file_name.endswith(".csv.gz"):
            csv_path = file_path.rstrip(".gz")

            # Decompress if .csv file doesn't exist
            if not os.path.exists(csv_path):
                decompress_file(file_path)

            # Load DataFrame
            df = pd.read_csv(csv_path)
            key_name = os.path.splitext(os.path.basename(csv_path))[0]
            data_frames[key_name] = df
            print(f"Loaded {key_name} (from .csv.gz) with shape {df.shape}")

        elif file_name.endswith(".csv"):
            df = pd.read_csv(file_path)
            key_name = os.path.splitext(os.path.basename(file_path))[0]
            data_frames[key_name] = df
            print(f"Loaded {key_name} (from .csv) with shape {df.shape}")

    # Quick inspection — show heads
    for name, df in data_frames.items():
        print(f"\n{name} preview:")
        print(df.head())

    return data_frames
def create_adata(sample_path, nucleus_genes_only=False):
    """
    Create an AnnData object from the 10x data files.

    Parameters:
    - sample_path: Path to the sample data folder.
    - nucleus_genes_only: Boolean indicating whether to only consider transcripts within the nucleolus.

    Returns:
    - adata: AnnData object created using the provided sample data.
    """
    # Construct the necessary file paths
    cell_feature_matrix_path = os.path.join(sample_path, "cell_feature_matrix.h5")
    cells_csv_path = os.path.join(sample_path, "cells.csv")
    cells_df = pd.read_csv(cells_csv_path)

    if nucleus_genes_only:
        # Construct the necessary file paths
        transcripts_csv_path = os.path.join(sample_path, "transcripts.csv")
        features_gz_path = os.path.join(sample_path, "cell_feature_matrix","features.tsv.gz")

        # Load additional file for nucleus consideration
        features_df = pd.read_csv(features_gz_path, sep='\t', compression='gzip', header=None)
        transcripts_df = pd.read_csv(transcripts_csv_path)
        
        # Filter for real genes in the assay
        genes = features_df[features_df[2] == "Gene Expression"][1].values

        # Remove transcripts that are not genes and not in the gene list
        transcripts_df_genes = transcripts_df[transcripts_df["feature_name"].isin(genes)]

        # Remove transcripts with qz<20
        transcripts_df_qv = transcripts_df_genes[transcripts_df_genes["qv"]>=20]

        # Remove unassigned transcripts
        transcripts_df_assigned = transcripts_df_qv[transcripts_df_qv["cell_id"] != "UNASSIGNED"]

        # Only keep transcripts that overlaps_nucleus is 1
        new_transcripts_df = transcripts_df_assigned[transcripts_df_assigned["overlaps_nucleus"] == 1]

        # Group by cell_id and count number of transcripts
        transcripts_df_assigned_overlaps_nucleus_grouped = new_transcripts_df.groupby("cell_id").size().reset_index(name="transcripts_count_nucleus")

        # Merge cells_df with transcripts_df_assigned_overlaps_nucleus_grouped
        cells_df_merged = pd.merge(cells_df, transcripts_df_assigned_overlaps_nucleus_grouped, on="cell_id", how="left")

        # Remove transcript_counts column and rename transcripts_count_nucleus as transcript_counts
        new_cells_df = cells_df_merged.drop(columns=["transcript_counts"])
        new_cells_df = new_cells_df.rename(columns={"transcripts_count_nucleus": "transcript_counts"})

        # Create a count matrix based on the new_transcripts_df
        count_matrix = new_transcripts_df.pivot_table(index='cell_id', columns='feature_name', aggfunc='size', fill_value=0)

        # Align cells_df with the count matrix by setting index and reindexing
        new_cells_df = new_cells_df.set_index('cell_id').reindex(count_matrix.index)

        # Create AnnData object
        adata = sc.AnnData(X=count_matrix.values, obs=new_cells_df, var=pd.DataFrame(index=count_matrix.columns))

        adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].to_numpy()
    else:
        adata = sc.read_10x_h5(cell_feature_matrix_path)
        cells_df.set_index(adata.obs_names, inplace=True)
        adata.obs = cells_df.copy()
        adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].to_numpy()
            
    return adata

# Function to convert shapely polygons to a format storable in AnnData
def polygon_to_coords(polygon):
    if polygon.is_empty:
        return None
    else:
        return list(polygon.exterior.coords)
    
def calculate_cell_boundaries(adata, transcripts_df):
    """
    Calculates cell boundaries from transcript data and updates the AnnData object.

    Parameters:
    - adata: AnnData object to be updated.
    - transcripts_df: pandas DataFrame containing transcript coordinates and cell IDs.

    Returns:
    - Updated AnnData object with cell boundaries stored in adata.uns['cell_boundaries'].
    """

    # Group by cell_id and calculate the boundary
    cell_boundaries = {}
    grouped = transcripts_df.groupby('cell_id')
    
    for cell_id, group in grouped:
        points = group[['x_location', 'y_location']].values
        if len(points) < 3:
            continue  # We need at least 3 points to calculate a boundary
        polygon = MultiPoint(points).convex_hull  # Use convex_hull to get the outer boundary
        cell_boundaries[cell_id] = polygon
        
    # Convert polygons to coordinate lists
    cell_boundaries_coords = {k: polygon_to_coords(v) for k, v in cell_boundaries.items()}
    
    # Add cell boundaries to AnnData
    adata.uns['cell_boundaries'] = cell_boundaries_coords
    
    return adata

def calculate_nucleus_boundaries(adata, nucleus_df):
    """
    Calculates nucleus boundaries from nucleus coordinates data and updates the AnnData object.

    Parameters:
    - adata: AnnData object to be updated.
    - nucleus_df: pandas DataFrame containing nucleus coordinates.

    Returns:
    - Updated AnnData object with nucleus boundaries stored in adata.uns['nucleus_boundaries'].
    """
    # Ensure that vertex_x and vertex_y are treated as numeric
    nucleus_df['vertex_x'] = pd.to_numeric(nucleus_df['vertex_x'])
    nucleus_df['vertex_y'] = pd.to_numeric(nucleus_df['vertex_y'])
    
    # Group by cell_id and create a Polygon for each nucleus
    nucleus_polygons = nucleus_df.groupby('cell_id').apply(
        lambda group: Polygon(zip(group['vertex_x'], group['vertex_y']))
    ).to_dict()

    # Convert polygons to coordinate lists
    nucleus_boundaries_coords = {k: polygon_to_coords(v) for k, v in nucleus_polygons.items()}

    # Update adata with nucleus boundaries
    adata.uns['nucleus_boundaries'] = nucleus_boundaries_coords
    
    return adata


# Function to calculate the bandwidth for kernel density estimation based on Silverman's rule of thumb
def calculate_bandwidth_silverman(data):
    return 1.06 * data.std() * len(data) ** (-1 / 5.)

def kde_cv(x, bandwidths):
    scores = []
    for bandwidth in bandwidths:
        # Compute cross validation score for each bandwidth
        score = compute_loo_score(x, bandwidth)
        scores.append(score)
    best_bandwidth = bandwidths[np.argmin(scores)]
    return best_bandwidth, scores

def compute_loo_score(x, bandwidth):
    scores = []
    for i in range(len(x)):
        # Leave one out: use all data points except the ith point
        training_data = np.delete(x, i)
        validation_point = x[i]
        
        # Create KDE with given bandwidth
        kde = gaussian_kde(training_data, bw_method=bandwidth)
        
        # Evaluate the KDE on the left out point
        estimated_density = kde.evaluate([validation_point])[0]
        
        # Actual density if using all data for estimation (this is a bit tricky, usually we don't do this step in LOO)
        kde_all = gaussian_kde(x, bw_method=bandwidth)
        actual_density = kde_all.evaluate([validation_point])[0]
        
        # Compute a score (squared error in this case)
        scores.append((estimated_density - actual_density) ** 2)

    # Return the mean of all squared errors
    return np.mean(scores)

# Function to perform Hartigans' Dip Test for unimodality
def perform_dip_test(data):
    dip_statistic, p_value = diptest(data.values)
    return dip_statistic, p_value

def analyze_gene_expressions(adata, gene, bandwidth=0.3, plot=True, filter_zeros=False):
    x = adata[:, gene].X.toarray().flatten()
    if filter_zeros:
        x = x[x > 0]
    
    if len(x) > 1:
        density = gaussian_kde(x, bw_method=bandwidth)
        xgrid = np.linspace(min(x), max(x), 1000)

        dip, p_value = perform_dip_test(pd.Series(x))
        modality = 'Unimodal' if p_value > 0.05 else 'Multimodal'
        peaks, _ = find_peaks(density(xgrid), distance=20)
        minima = argrelextrema(density(xgrid), np.less)[0]

        highest_peak_xpos = None
        last_min_before_peak = None

        if peaks.size > 0:
            highest_peak = np.argmax(density(xgrid[peaks]))
            highest_peak_xpos = peaks[highest_peak]

        if highest_peak_xpos is not None and minima.size > 0:
            relevant_minima = minima[minima < highest_peak_xpos]
            if relevant_minima.size > 0:
                last_min_before_peak = relevant_minima[-1]

        if plot:
            plt.figure(figsize=(8, 4))
            plt.title(f'Expression Distribution for {gene}')
            plt.plot(xgrid, density(xgrid), label='Density')
            plt.hist(x, bins=50, alpha=0.3, density=True, label='Histogram')
            if highest_peak_xpos is not None:
                plt.axvline(x=xgrid[highest_peak_xpos], color='blue', linestyle=':', label='Highest Peak')
            if last_min_before_peak is not None:
                plt.axvline(x=xgrid[last_min_before_peak], color='red', linestyle='--', label='Minima before Peak')
            plt.xlabel('Expression Level')
            plt.ylabel('Density')
            plt.legend(title=f"Modality: {modality} (p-value={p_value:.2f})")
            plt.show()
            print(f"{gene} has a {modality} modality (p-value={p_value:.2f})")
            if modality == 'Multimodal':
                if last_min_before_peak is not None:
                    print(f"Background threshold for {gene} is {xgrid[last_min_before_peak]:0.3f}")

        if modality == 'Multimodal' and last_min_before_peak is not None:
            return xgrid[last_min_before_peak]

    return None

def plot_spatial_genes(adata, genes, cell_boundaries=False, nucleus_boundaries=False, cmap='light_dark_red'):
    """
    Plots the spatial distribution of transcripts and overlays cell and nucleus boundaries.

    Parameters:
    - adata: AnnData object containing spatial coordinates.
    - genes: List of genes to plot.
    - show_all_transcripts: Whether to show all transcripts or only those associated with cells.
    - cmap: Colormap to use for plotting (default is 'light_dark_blue').
    """
    plt.figure(figsize=(8, 6))
    plt.title('Spatial Distribution of Transcripts')

    if cell_boundaries:
        for cell_id, boundary in adata.uns['cell_boundaries'].items():
            if boundary is not None:
                coords = np.array(boundary)
                plt.plot(coords[:, 0], -coords[:, 1], color='blue', alpha=0.1)

    if nucleus_boundaries:
        for cell_id, boundary in adata.uns['nucleus_boundaries'].items():
            if boundary is not None:
                coords = np.array(boundary)
                plt.plot(coords[:, 0], -coords[:, 1], color='red', alpha=0.1)

    coordinates = adata.obsm['spatial']
        
    for gene in genes:
        if gene in adata.var_names:
            gene_expression = adata[:, gene].X
            if issparse(gene_expression):
                gene_expression = gene_expression.toarray().flatten()  # Ensure it's flattened

            gene_expression = gene_expression.flatten()  # Ensure it's flattened
            
            x = coordinates[:, 0]
            y = -coordinates[:, 1]

            # filter out zero expression values
            mask = gene_expression > 0
            x = x[mask]
            y = y[mask]
            gene_expression = gene_expression[mask]

            # Plotting the spatial distribution of transcripts
            plt.scatter(x, y, c=gene_expression, cmap=cmap, s=1, label=gene)  # Use custom colormap
            plt.colorbar(label='Gene Expression Level')

        else:
            print(f"{gene} gene is not found in the dataset.")

    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.show()

def plot_spatial_transcripts(adata, transcripts_df, genes, cell_boundaries=False, nucleus_boundaries=False):
    """
    Plots the spatial distribution of transcripts and overlays cell and nucleus boundaries.
    
    Parameters:
    - adata: AnnData object containing spatial coordinates.
    - transcripts_df: DataFrame containing transcript information.
    - genes: List of genes to plot.
    - cell_boundaries: Whether to show cell boundaries.
    - nucleus_boundaries: Whether to show nucleus boundaries.
    """
    plt.figure(figsize=(8, 8))

    # Plot cell boundaries
    if cell_boundaries:
        # check if cell boundaries are present in adata.uns
        if 'cell_boundaries' not in adata.uns:
            print('Cell boundaries not found in adata.uns. Please run calculate_cell_boundaries() function first.')
        else:
            for cell_id, boundary in adata.uns['cell_boundaries'].items():
                if boundary is not None:
                    coords = np.array(boundary)
                    plt.plot(coords[:, 0], -coords[:, 1], color='blue', alpha=0.1)

    # Plot nucleus boundaries
    if nucleus_boundaries:
        # check if nucleus boundaries are present in adata.uns
        if 'nucleus_boundaries' not in adata.uns:
            print('Nucleus boundaries not found in adata.uns. Please run calculate_nucleus_boundaries() function first.')
        else:
            for cell_id, boundary in adata.uns['nucleus_boundaries'].items():
                if boundary is not None:
                    coords = np.array(boundary)
                    plt.plot(coords[:, 0], -coords[:, 1], color='gray', alpha=0.1)

    # Subsetting the DataFrame for the selected genes
    transcripts_df_4plot = transcripts_df[transcripts_df['feature_name'].isin(genes) & transcripts_df['cell_id'].isin(adata.obs_names)]

    # Get unique feature names
    feature_names = transcripts_df_4plot['feature_name'].unique()
    colors = plt.cm.jet(np.linspace(0, 1, len(feature_names)))
    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', '+', 'x'] * (len(feature_names) // 11 + 1)
    
    # Map feature names to colors and shapes
    feature_color_shape_map = {feature_name: (color, marker) for feature_name, color, marker in zip(feature_names, colors, markers[:len(feature_names)])}
    
    # Create the scatter plot for transcripts
    for feature_name in feature_names:
        subset = transcripts_df_4plot[transcripts_df_4plot['feature_name'] == feature_name]
        plt.scatter(
            subset['x_location'],
            -subset['y_location'],
            s=0.01,  # Adjust size as needed
            color=feature_color_shape_map[feature_name][0],  # color 
            marker=feature_color_shape_map[feature_name][1],  # shape
            label=feature_name
        )

    # Set plot labels and title
    plt.xlabel('x_location')
    plt.ylabel('y_location')
    plt.title('Spatial Distribution of Transcripts')

    # Create legend with unique entries
    handles, labels = plt.gca().get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    plt.legend(*zip(*unique), markerscale=10)  # Increase markerscale if needed

    # Show the plot
    plt.show()

def process_cells_link(link):
    data = []
    match = re.search(r'&target=([\d._]+)&', link)
    if match:
        coordinates = match.group(1).split('_')
        x_coord, y_coord = coordinates[0], coordinates[1]
        data.append({'X_coordinate': x_coord, 'Y_coordinate': y_coord})
    else:
        print('No X and Y coordinates found in link.')
    cell_info = pd.DataFrame(data)
    return cell_info

def decompress_file(input_path):
    """
    Decompresses a .gz file.

    Parameters:
        input_path (str): The path to the .gz file to decompress.

    Returns:
        str: The path to the decompressed file.
    """
    output_path = os.path.splitext(input_path)[0]

    # Decompress the file
    with gzip.open(input_path, 'rt') as compressed_file:
        with open(output_path, 'w') as decompressed_file:
            decompressed_file.write(compressed_file.read())

    return output_path