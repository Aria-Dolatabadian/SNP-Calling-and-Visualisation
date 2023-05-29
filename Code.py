
import csv
from Bio import SeqIO

def find_snps(seq1, seq2):
    snps = []
    transitions = ['AG', 'GA', 'CT', 'TC']
    transversions = ['AC', 'CA', 'AT', 'TA', 'GC', 'CG', 'GT', 'TG']

    for i, (base1, base2) in enumerate(zip(seq1, seq2)):
        if base1 != base2:
            position = i + 1
            snp = f"{base1}{position}{base2}"
            if base1 + base2 in transitions:
                snps.append((position, base1, base2, 'Transition'))
            elif base1 + base2 in transversions:
                snps.append((position, base1, base2, 'Transversion'))

    return snps

# Read FASTA files
fasta_file1 = 'file1.fasta'
fasta_file2 = 'file2.fasta'

sequences1 = SeqIO.parse(fasta_file1, 'fasta')
sequences2 = SeqIO.parse(fasta_file2, 'fasta')

# Compare sequences and find SNPs
all_snps = []
for seq1, seq2 in zip(sequences1, sequences2):
    snps = find_snps(seq1.seq, seq2.seq)
    all_snps.extend(snps)

# Write SNPs to CSV file
csv_file = 'snps.csv'
header = ['Position', 'Reference', 'Variant', 'Classification']

with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(header)
    writer.writerows(all_snps)

print(f"SNPs exported to {csv_file}")

import csv
import matplotlib.pyplot as plt

# Read SNPs from CSV file
csv_file = 'snps.csv'

transitions_count = 0
transversions_count = 0

with open(csv_file, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)  # Skip the header row

    for row in reader:
        _, _, _, classification = row
        if classification == 'Transition':
            transitions_count += 1
        elif classification == 'Transversion':
            transversions_count += 1

# Create labels and counts for the pie chart
labels = ['Transitions', 'Transversions']
counts = [transitions_count, transversions_count]

# Plotting the pie chart
plt.figure(figsize=(6, 6))
plt.pie(counts, labels=labels, autopct='%1.1f%%', startangle=90)
plt.title('Transition and Transversion SNPs')

plt.show()


import csv
import matplotlib.pyplot as plt
import numpy as np

# Read SNPs from CSV file
csv_file = 'snps.csv'

snp_types = ['AT', 'TA', 'AC', 'CA', 'CT', 'TC', 'GA', 'AG', 'CG', 'GC', 'TG', 'GT']
snp_counts = {snp_type: 0 for snp_type in snp_types}

with open(csv_file, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)  # Skip the header row

    for row in reader:
        _, reference, variant, _ = row
        snp_type = reference + variant
        if snp_type in snp_types:
            snp_counts[snp_type] += 1

# Create labels and counts for the stacked bar plot
labels = ['SNP Types']
counts = [snp_counts[snp_type] for snp_type in snp_types]

# Set up colors for each SNP type
colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple',
          'cyan', 'magenta', 'gray', 'brown', 'pink', 'lime']

# Plotting the stacked bar plot
plt.figure(figsize=(8, 6))

# Calculate the cumulative counts for each SNP type
cumulative_counts = np.cumsum(counts)

# Plot the stacked bars
bar_width = 0.8
bottom = np.zeros(len(counts))

for i in range(len(counts)):
    plt.bar(labels, counts[i], width=bar_width, bottom=bottom, color=colors[i], edgecolor='black')
    bottom += counts[i]

plt.xlabel('SNP Types')
plt.ylabel('Count')
plt.title('Counts of Different SNP Types')

# Create a legend for the SNP colors
legend_elements = [plt.Rectangle((0, 0), 1, 1, color=color, label=snp_type)
                   for snp_type, color in zip(snp_types, colors)]

# Move the legend out of the plot
plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()
plt.show()

