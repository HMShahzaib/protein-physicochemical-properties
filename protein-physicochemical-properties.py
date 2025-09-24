from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import sys

def compute_protein_properties(fasta_file, output_csv):
    results = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        analysis = ProteinAnalysis(sequence)

        molecular_weight = analysis.molecular_weight()
        theoretical_pI = analysis.isoelectric_point()
        instability_index = analysis.instability_index()
        gravy = analysis.gravy()

        # Compute amino acid composition
        aa_composition = analysis.count_amino_acids()
        A = aa_composition.get('A', 0)  # Alanine
        V = aa_composition.get('V', 0)  # Valine
        I = aa_composition.get('I', 0)  # Isoleucine
        L = aa_composition.get('L', 0)  # Leucine
        total_residues = sum(aa_composition.values())

        # Calculate aliphatic index
        if total_residues > 0:
            aliphatic_index = ((A + 2.9 * V + 3.9 * (I + L)) / total_residues) * 100
        else:
            aliphatic_index = 0

        results.append([
            record.id, molecular_weight, theoretical_pI,
            instability_index, aliphatic_index, gravy
        ])

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Protein_ID", "Molecular_Weight", "Theoretical_pI", "Instability_Index", "Aliphatic_Index", "GRAVY"])
        writer.writerows(results)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.csv")
        sys.exit(1)

    compute_protein_properties(sys.argv[1], sys.argv[2])
