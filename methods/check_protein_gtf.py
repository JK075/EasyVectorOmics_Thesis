import pandas as pd


def read_tsv(tsv_file, output_file):
    """
    Liest eine TSV-Datei als Pandas DataFrame ein, filtert sie und speichert die gefilterte Version.
    """
    df = pd.read_csv(tsv_file, sep='\t')

    print("Spalten im DataFrame:", df.columns.tolist())

    # Behalte nur das Protein mit der höchsten Länge pro Gene ID (wenn mehrere existieren)
    df = df.sort_values(by=["Gene ID", "Protein length"], ascending=[True, False])
    df = df.drop_duplicates(subset=["Gene ID"], keep="first")

    # Wähle nur die Spalten Gene ID und Protein accession
    df = df[["Gene ID", "Protein accession"]]

    # Speichere das gefilterte DataFrame als neue TSV-Datei
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Gefilterte TSV-Datei wurde gespeichert als: {output_file}")

    return df


def filter_fasta(fasta_file, filtered_tsv, output_fasta):
    """
    Liest eine FASTA-Datei ein und entfernt Sequenzen, deren Protein accession nicht in der gefilterten TSV-Datei vorkommt.
    """
    # Lade die gefilterten Protein Accessions
    df = pd.read_csv(filtered_tsv, sep='\t')
    valid_accessions = set(df["Protein accession"].tolist())

    # Lese die FASTA-Datei und filtere die Einträge
    with open(fasta_file, 'r') as f_in, open(output_fasta, 'w') as f_out:
        write_entry = False
        sequence_buffer = []
        for line in f_in:
            if line.startswith('>'):
                # Falls vorher eine gültige Sequenz zwischengespeichert wurde, schreibe sie ins Output-File
                if write_entry and sequence_buffer:
                    f_out.writelines(sequence_buffer)

                # Neue Header-Zeile -> Prüfe, ob sie behalten wird
                accession = line.split()[0][1:]  # Entferne das '>' und nehme die ID
                write_entry = accession in valid_accessions

                # Setze den Puffer zurück und speichere den neuen Header
                sequence_buffer = [line]
            else:
                # Füge die Sequenzzeilen dem Puffer hinzu
                sequence_buffer.append(line)

        # Falls die letzte Sequenz im Puffer gültig war, schreibe sie
        if write_entry and sequence_buffer:
            f_out.writelines(sequence_buffer)

    print(f"Gefilterte FASTA-Datei wurde gespeichert als: {output_fasta}")


lupus = r"D:\EasyVectorOmics\Paper\thesis\Canis lupus\ncbi_dataset.tsv"
lupus_protein = r"D:\EasyVectorOmics\Paper\thesis\Canis lupus\protein.faa"
lupus_output = r"D:\EasyVectorOmics\Paper\thesis\Canis lupus\Canis_lupus_protein.faa"

monkey = r"D:\EasyVectorOmics\Paper\thesis\Macaca fascicularis\ncbi_dataset.tsv"
monkey_protein = r"D:\EasyVectorOmics\Paper\thesis\Macaca fascicularis/protein.faa"
monkey_output = r"D:\EasyVectorOmics\Paper\thesis\Macaca fascicularis/Macaca_fascicularis_protein.faa"

mouse = r"D:\EasyVectorOmics\Paper\thesis\Mus musculus\ncbi_dataset.tsv"
mouse_protein = r"D:\EasyVectorOmics\Paper\thesis\Mus musculus\protein.faa"
mouse_output = r"D:\EasyVectorOmics\Paper\thesis\Mus musculus\Mus_musculus_protein.faa"

rat = r"D:\EasyVectorOmics\Paper\thesis\Rattus norvegicus\ncbi_dataset.tsv"
rat_protein = r"D:\EasyVectorOmics\Paper\thesis\Rattus norvegicus\protein.faa"
rat_output = r"D:\EasyVectorOmics\Paper\thesis\Rattus norvegicus\Rattus_norvegicus_protein.faa"

# Beispielaufruf
tsv_file = rat  # Ersetze mit deinem Dateipfad
output_file = r"C:\Users\Jonathan\Downloads\filtered_example.tsv"  # Ersetze mit gewünschtem Ausgabe-Dateipfad
df = read_tsv(tsv_file, output_file)

# FASTA-Filterung
fasta_file = rat_protein # Ersetze mit deinem FASTA-Dateipfad
output_fasta = rat_output  # Ersetze mit gewünschtem Ausgabe-Dateipfad
filter_fasta(fasta_file, output_file, output_fasta)