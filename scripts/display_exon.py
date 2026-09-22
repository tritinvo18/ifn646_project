import requests
import gzip

target_gene = 'SERPINB2'
target_chrom = '18'
chrom_accession = 'NC_000018.10'
gff_file = 'GCF_000001405.26_GRCh38_genomic.gff.gz'

target_exons = []
with gzip.open(gff_file, 'rt') as file:
    for line in file:
        if line.startswith('#'): continue
        if '\texon\t' in line and f'gene={target_gene};' in line and 'NM_001143818.1' in line:
            parts = line.strip().split('\t')
            if parts[0] == chrom_accession:
                target_exons.append({'start': int(parts[3]), 'end': int(parts[4]), 'strand': parts[6]})

first_exon = sorted(target_exons, key=lambda x: x['end'] if x['strand'] == '-' else x['start'], reverse=(target_exons[0]['strand'] == '-'))[0]
strand = first_exon['strand']
start = first_exon['start']
end = first_exon['end']

server = "https://rest.ensembl.org"
endpoint = f"/sequence/region/human/{target_chrom}:{start}..{end}:{1 if strand == '+' else -1}"
response = requests.get(server + endpoint, headers={"Content-Type": "text/plain"})

print(f"--- SERPINB2 First Exon (Canonical Transcript NM_001143818.1) ---")
print(f"Coordinates: chr{target_chrom}:{start}-{end} (Strand {strand})")
print(f"Length: {end - start + 1} bp\n")
print(response.text)
