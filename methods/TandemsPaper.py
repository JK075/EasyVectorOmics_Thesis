import pandas as pd
import re
import numpy as np



def load_tsv_as_dict(tsv_file):
    """Liest eine TSV-Datei ein und gibt ein Dictionary mit GeneID als Schlüssel und (ProteinID, Spezies) als Werte zurück."""
    gene_dict = {}

    with open(tsv_file, 'r') as file:
        next(file)  # Überspringt die Header-Zeile
        for line in file:
            columns = line.strip().split("\t")
            if len(columns) >= 3:  # Sicherstellen, dass alle Spalten vorhanden sind
                gene_id = str(columns[0]).strip()  # Entferne Leerzeichen und stelle sicher, dass gene_id ein String ist
                protein_id = str(columns[1]).strip()  # Entferne Leerzeichen und stelle sicher, dass protein_id ein String ist
                species = str(columns[2]).strip()  # Entferne Leerzeichen und stelle sicher, dass species ein String ist
                gene_dict[gene_id] = [protein_id, species]

    return gene_dict



def process_gtf(gtf_file):

    """
        Processes a GTF file to filter and extract gene rows.

        Functionality:
        This function reads a GTF file and selects only the rows that correspond to gene features.
        It uses pandas to efficiently filter the data during the reading process, which enhances performance for large files.

        Parameters:
        - gtf_file (str): The path to the GTF file to be processed.

        Returns:
        - DataFrame: A pandas DataFrame containing only the rows where the feature is 'gene'.
        """

    # Filter "gene" rows directly while reading the file for improved speed
    df = pd.read_csv(gtf_file, sep='\t', header=None, comment='#', dtype=str)

    return df[df[2] == 'gene']


def get_harmonic_mean(df):
    # Extrahiere die benötigten Spalten als NumPy-Arrays
    qstart = df[2].to_numpy()
    qend = df[3].to_numpy()
    qlen = df[4].to_numpy()
    sstart = df[5].to_numpy()
    send = df[6].to_numpy()
    slen = df[7].to_numpy()
    pident = df[8].to_numpy()

    # Berechne den Overlap (vektorisiert)
    denominator = qlen + slen

    # Sicherheitscheck: Denominator darf nicht Null sein
    assert np.all(denominator != 0), "Error: Division durch Null im Overlap-Berechnungsprozess!"

    overlap = ((qend - qstart) + (send - sstart)) / denominator * 100

    # Berechne den harmonischen Mittelwert direkt
    harmonic = 2 * (pident * overlap) / (pident + overlap)

    # Füge die berechnete Spalte dem DataFrame hinzu
    df['harmonic_mean'] = harmonic

    return df


def filter_blast_results(blast_file, protein_dict, threshold=0.85):

    """
       Filters BLAST results based on a threshold and protein IDs present in the protein dictionary.

       Functionality:
       This function reads BLAST results from a specified file and filters them based on a given threshold score.
       It checks if both the query and target protein IDs exist in the provided protein dictionary,
       and it retains only those results where the species of both proteins match and their IDs are not equal.

       Parameters:
       - blast_file (str): The path to the BLAST result file to be filtered.
       - protein_dict (dict): A dictionary mapping protein IDs to their parent and species.
       - threshold (float): A minimum percentage threshold for filtering BLAST results (default is 0.85).

       Returns:
       - DataFrame: A pandas DataFrame containing filtered BLAST results, with columns for the query and target IDs only.
       """

    threshold *= 100
    df = pd.read_csv(blast_file, sep='\t', header=None)

    df = get_harmonic_mean(df)

    df = df[df['harmonic_mean'] > threshold]

    # Überprüfe nur die ersten 5 Einträge der Spalte 0

    df[0] = df[0].astype(str)
    df[1] = df[1].astype(str)
    # Filter with vectorized lookup for improved speed
    df = df[df[0].isin(protein_dict) & df[1].isin(protein_dict)]

    df['query_species'] = [protein_dict[str(x)][1] for x in df[0]]
    df['target_species'] = [protein_dict[str(x)][1] for x in df[1]]

    # Keep only rows where species match and IDs are not equal
    df = df[(df['query_species'] == df['target_species']) & (df[0] != df[1])]

    return df.drop(columns=['query_species', 'target_species'])


def process_row(row, processed_df, threshold):
    """
    Processes a single row of filtered BLAST results to determine if two genes are in tandem.

    Parameters:
    - row (tuple): A tuple containing the query and target IDs from the filtered BLAST results.
    - processed_df (DataFrame): A DataFrame that includes processed gene information extracted from a .gtf file, detailing the locations of various genes.
    - threshold (int): The distance threshold (in kilobases) to consider two genes as tandem.

    Returns:
    - tuple: A tuple (query_id, target_id) if the genes are in tandem; otherwise, None.
    """
    query_id, target_id = row[0], row[1]

    # Check if both query and target IDs exist in processed_df
    if query_id in processed_df.index and target_id in processed_df.index:
        query_rows = processed_df.loc[query_id]
        target_rows = processed_df.loc[target_id]

        # Falls nur eine Zeile zurückkommt, sicherstellen, dass wir trotzdem iterieren können
        if isinstance(query_rows, pd.Series):
            query_rows = query_rows.to_frame().T  # In DataFrame umwandeln
        if isinstance(target_rows, pd.Series):
            target_rows = target_rows.to_frame().T  # In DataFrame umwandeln

        # Alle möglichen Kombinationen von query_row und target_row testen
        for _, query_row in query_rows.iterrows():
            for _, target_row in target_rows.iterrows():
                #print(f"Vergleich: target {target_row[0]} mit query {query_row[0]}")

                # Check if the genes are on the same chromosome
                if query_row[0] == target_row[0]:
                    query_start, query_end = int(query_row[3]), int(query_row[4])
                    target_start, target_end = int(target_row[3]), int(target_row[4])

                    # Check if genes are within the threshold distance
                    if ((query_start - threshold * 1000 <= target_end and query_start >= target_start) or
                            (query_end + threshold * 1000 >= target_start and query_end <= target_end)):
                        return query_id, target_id  # Sobald ein Treffer gefunden wird, abbrechen
    return None


def find_tandems(threshold, filtered_blast_df, processed_df):
    """
    Finds tandem gene pairs from the filtered BLAST results.

    Parameters:
    - threshold (int): The distance threshold (in kilobases) to consider two genes as tandem.
    - filtered_blast_df (DataFrame): A DataFrame containing filtered BLAST results.
    - processed_df (DataFrame): A DataFrame containing processed gene information.

    Returns:
    - list: A list of tuples containing pairs of query and target IDs that are in tandem.
    """

    # Set index on processed_df for efficient row lookup
    processed_df['gene_id'] = processed_df[8].str.extract(r"GeneID:(\d+)(?=\D|$)")
    processed_df.set_index('gene_id', inplace=True)

    matching_genid_pairs = []

    # Iterate over the rows of filtered_blast_df
    for row in filtered_blast_df.itertuples(index=False):
        result = process_row(row, processed_df, threshold)
        if result:
            matching_genid_pairs.append(result)

    print(f"Number of matching gene pairs: {len(matching_genid_pairs)}")
    return matching_genid_pairs


def process_genid_pairs(matching_genid_pairs, result_dict):
    """
    Processes tandem gene pairs and groups connected genes into rows based on species.

    Optimized version using dictionaries for faster lookups.

    Parameters:
    - matching_genid_pairs (list): List of tuples with matching query and target IDs.
    - result_dict (dict): Dictionary with protein IDs as keys and lists [parent, species] as values.

    Returns:
    - DataFrame: A pandas DataFrame with columns for Tandems and Species.
    """

    # Dictionary to keep track of which genes are already grouped
    gene_to_row_map = {}
    results = []

    # Iterate over each pair in matching_genid_pairs
    for query_id, target_id in matching_genid_pairs:
        # Get species names for the query from result_dict
        species_name = result_dict.get(query_id, [None, None])[1]

        # Check if either gene is already in the results dictionary
        row_index = gene_to_row_map.get(query_id) or gene_to_row_map.get(target_id)

        if row_index is None:
            # If neither gene is found, create a new row
            new_tandems = {query_id, target_id}
            results.append({
                'Tandems': new_tandems,
                'Species': species_name
            })
            # Store references in the dictionary
            row_index = len(results) - 1
            gene_to_row_map[query_id] = row_index
            gene_to_row_map[target_id] = row_index
        else:
            # If one of the genes is found, update the existing row
            existing_row = results[row_index]
            existing_row['Tandems'].update([query_id, target_id])
            # Update the dictionary for the new genes
            gene_to_row_map[query_id] = row_index
            gene_to_row_map[target_id] = row_index

    # Convert the sets to sorted comma-separated strings
    for row in results:
        row['Tandems'] = ', '.join(sorted(row['Tandems']))

    # Convert the list of results to a DataFrame
    output_df = pd.DataFrame(results, columns=['Tandems', 'Species'])

    return output_df


def calculate_similarity(row):
    """
    Dummy scoring function to calculate similarity score between two sequences.
    Replace this with your own implementation.
    """
    return row[10] #bitscore


def filter_genid_with_highest_mean_score(df_tandems, df_blast):
    """
    Filters the genid with the highest mean score based on similarity scores from the BLAST file.

    Parameters:
    df_tandems (pd.DataFrame): Input DataFrame with genids in the first column (comma-separated) and species in the second column.
    df_blast (pd.DataFrame): BLAST file DataFrame with at least two columns: [target, query].

    Returns:
    pd.DataFrame: Updated DataFrame with the genid having the highest mean score in the third column.
    """
    # Ensure columns are named properly for processing
    df_tandems.columns = ['genids', 'species']

    df_blast[df_blast.columns[0]] = df_blast[df_blast.columns[0]].astype(str)
    df_blast[df_blast.columns[1]] = df_blast[df_blast.columns[1]].astype(str)

    # Rename the first two columns of df_blast to 'target' and 'query', keeping others intact
    df_blast = df_blast.rename(columns={df_blast.columns[0]: 'target', df_blast.columns[1]: 'query'})

    # Initialize a column for storing the highest scoring genid
    df_tandems['highest_mean_genid'] = None

    for index, row in df_tandems.iterrows():
        # Split the comma-separated genids into a list
        genid_list = [genid.strip() for genid in row['genids'].split(',')]
        #print(genid_list)

        highest_mean_score = -1
        best_genid = None

        for genid in genid_list:
            # Filter rows in the BLAST file where the genid is the target
            matching_rows = df_blast[df_blast['target'] == genid]
            print(matching_rows)

            if not matching_rows.empty:
                # Calculate the similarity score for each row
                scores = matching_rows.apply(lambda x: calculate_similarity(x), axis=1)
                print(scores)
                n = len(scores)
                product_of_scores = (1 / n) * scores.iloc[0]

                for score in scores.iloc[1:]:  # Iteriere über alle Scores ab der zweiten Zahl
                    product_of_scores *= score

                # Update the best genid if the current mean score is higher
                if product_of_scores > highest_mean_score:
                    highest_mean_score = product_of_scores
                    best_genid = genid

        # Store the best genid with the highest mean score
        df_tandems.at[index, 'highest_mean_genid'] = best_genid

    return df_tandems



# Main execution



#-------Files------------------------------------------------------------------

#-------Test-----------------------------------------------------------------
"""
fasta_file = r"D:\EasyVectorOmics\Tandems\Test\fasta_test.txt"
gtf_file = r"D:\EasyVectorOmics\Tandems\Test\test_gtf.txt"
blast_file = r"D:\EasyVectorOmics\Tandems\Test\blast_test.txt"
tsv_file = r"D:\EasyVectorOmics\Tandems\Test\tsv_test.txt"
"""
#-----normal------------------------------------------------------


tsv_file = "D:\EasyVectorOmics\TEST\map_gen_protein_species.tsv"
gtf_file = r"D:\EasyVectorOmics\Paper\thesis\gtf\gtf_Canis.gtf"
blast_file = r"D:\EasyVectorOmics\TEST\blast_results.tsv"
output_path = r"D:\EasyVectorOmics\TEST\tandems_out\Tandems_output.tsv"
#tsv_file = r"D:\EasyVectorOmics\Tandems\Orthogroups.tsv"


#----thresholds----------------------------------------------------------------

threshold_blast = 0.89
threshold_kilobase = 200

#-----Methods-------------------------------------------------------------------

print("Start computing Fasta File...")
result_dict = load_tsv_as_dict(tsv_file)
print("done computing Fasta File")

print("start extracting gtf File ...")
processed_df = process_gtf(gtf_file)
print("done extracting gtf File")

print("start filtering Blast File ...")
filtered_blast_df = filter_blast_results(blast_file, result_dict, threshold_blast)
print("done filtering Blast File")

print("start finding tandems ...")
matching_genid_pairs = find_tandems(threshold_kilobase, filtered_blast_df, processed_df)
print("done finding tandems")

print("computing tandems ...")
output_df = process_genid_pairs(matching_genid_pairs, result_dict)

df_blast = pd.read_csv(blast_file, sep='\t', header=None)

final_df = filter_genid_with_highest_mean_score(output_df, df_blast)
output_df.to_csv(output_path, sep='\t', index=False)
print("Done, created outputfile")


#save_dir = 'D:\EasyVectorOmics\Tandems'
#analyze_tandem_clusters(output_df, save_dir)
