#!/usr/bin/env python3
"""
This script creates legacy output from plasmidfinder json output.
@author: Kutluhan Incekara
email: kutluhan.incekara@ct.gov
"""

import json
with open('data.json', 'r') as f:
    data = json.load(f)
with open('results_tab.tsv', 'w') as f:
    headers = ['Database', 'Plasmid', 'Identity', 'Query / Template length', 'Contig', 'Position in contig', 'Note', 'Accession number']
    f.write('\t'.join(headers) + '\n')
    for seq_region in data['seq_regions'].values():
        # Extract required information
        database = seq_region['ref_database'][0]
        plasmid = seq_region['name']
        identity = f"{seq_region['identity']:.1f}"
        lengths = f"{seq_region['alignment_length']} / {seq_region['ref_gene_lenght']}"
        contig = seq_region['query_id']
        position = f"{seq_region['query_start_pos']}..{seq_region['query_end_pos']}"
        note = seq_region['note']
        accession = seq_region['ref_acc']        
        # Write the row
        row = [database, plasmid, identity, lengths, contig, position, note, accession]
        f.write('\t'.join(row) + '\n')
# Extract plasmid names from seq_regions
plasmids = [region['name'] for region in data['seq_regions'].values()]
with open('PLASMIDS', 'w') as f:
    f.write(','.join(plasmids))