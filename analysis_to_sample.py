import argparse
import pandas as pd
import requests
import io
import sys
from typing import List, Dict, Optional

def fetch_data(url: str, access_key: str, secret_key: str) -> pd.DataFrame:
    """Helper function to fetch data from API."""
    try:
        response = requests.get(url, auth=(access_key, secret_key))
        response.raise_for_status()
        lines = response.text.split('\n')
        tsv_data = '\n'.join(line for line in lines if not line.startswith('20'))
        df = pd.read_csv(io.StringIO(tsv_data), sep='\t')
        df.columns = df.columns.str.strip()
        return df
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}", file=sys.stderr)
        sys.exit(1)

def fetch_analysis_sets(analysis_set: str, access_key: str, secret_key: str) -> pd.DataFrame:
    """Fetch measurement sets and auxiliary sets for the given analysis set."""
    url = f"https://api.data.igvf.org/multireport.tsv?type=MeasurementSet&input_for=%2Fanalysis-sets%2F{analysis_set}%2F&field=%40id&field=auxiliary_sets&field=input_for"
    df = fetch_data(url, access_key, secret_key)
    
    # Process data
    df['analysis_set_id'] = analysis_set
    df['measurement_sets'] = df['ID'].str.extract(r'/measurement-sets/([^/]+)/')
    
    # Rename and select columns
    return df.rename(columns={'Auxiliary Sets': 'associated_auxiliary_sets'})[
        ['analysis_set_id', 'measurement_sets', 'associated_auxiliary_sets']
    ]

def get_sequence_files(result: pd.DataFrame, access_key: str, secret_key: str) -> pd.DataFrame:
    """Fetch sequence files from measurement sets and auxiliary sets."""
    all_tables = []
    
    for _, row in result.iterrows():
        measurement_id = row['measurement_sets']
        # Fetch measurement set sequence files
        url = f"https://api.data.igvf.org/multireport.tsv?type=SequenceFile&file_set.@id=%2Fmeasurement-sets%2F{measurement_id}%2F&illumina_read_type=*&field=%40id&field=accession&field=illumina_read_type&field=lane&field=md5sum&field=flowcell_id&field=seqspecs"
        df = fetch_data(url, access_key, secret_key)
        df['measurement_sets'] = measurement_id
        df['file_modality'] = 'scRNA'
        all_tables.append(df)
        
        # Fetch auxiliary set sequence files if available
        if pd.notna(row['associated_auxiliary_sets']):
            aux_sets = str(row['associated_auxiliary_sets']).strip().split(',')
            aux_sets = [s.strip() for s in aux_sets]
            
            if len(aux_sets) == 1:
                aux_id = aux_sets[0].split('/')[-2]
                url = f"https://api.data.igvf.org/multireport.tsv?type=SequenceFile&file_set.@id=%2Fauxiliary-sets%2F{aux_id}%2F&illumina_read_type=*&field=%40id&field=accession&field=illumina_read_type&field=lane&field=md5sum&field=seqspecs"
                df = fetch_data(url, access_key, secret_key)
                df['measurement_sets'] = measurement_id
                df['file_modality'] = 'gRNA'
                all_tables.append(df)
            else:
                aux_types = ['gRNA', 'hash']
                for aux_type, aux_set in zip(aux_types, aux_sets):
                    if pd.notna(aux_set):
                        aux_id = aux_set.split('/')[-2]
                        url = f"https://api.data.igvf.org/multireport.tsv?type=SequenceFile&file_set.@id=%2Fauxiliary-sets%2F{aux_id}%2F&illumina_read_type=*&field=%40id&field=accession&field=illumina_read_type&field=lane&field=md5sum&field=seqspecs"
                        df = fetch_data(url, access_key, secret_key)
                        df['measurement_sets'] = measurement_id
                        df['file_modality'] = aux_type
                        all_tables.append(df)
    
    return pd.concat(all_tables, ignore_index=True) if all_tables else pd.DataFrame()

def rearrange_sequence_files(sequence_files: pd.DataFrame) -> pd.DataFrame:
    """Pivot sequence files to wide format based on read type."""
    if sequence_files.empty:
        return pd.DataFrame()
    
    df = sequence_files[sequence_files['Illumina Read Type'].isin(['R1', 'R2'])].copy()
    
    r1_data = df[df['Illumina Read Type'] == 'R1'].copy()
    r2_data = df[df['Illumina Read Type'] == 'R2'].copy()
    
    r1_data = r1_data.rename(columns={'Accession': 'R1_path', 'MD5sum': 'R1_md5sum'})
    r2_data = r2_data.rename(columns={'Accession': 'R2_path', 'MD5sum': 'R2_md5sum'})
    
    r1_cols = ['Lane', 'measurement_sets', 'file_modality', 'R1_path', 'R1_md5sum', 'Flowcell ID', 'Seqspecs']
    r2_cols = ['Lane', 'measurement_sets', 'file_modality', 'R2_path', 'R2_md5sum', 'Flowcell ID', 'Seqspecs']
    
    r1_data = r1_data[r1_cols]
    r2_data = r2_data[r2_cols]
    
    merged_df = pd.merge(
        r1_data, r2_data,
        on=['Lane', 'measurement_sets', 'file_modality', 'Flowcell ID', 'Seqspecs'],
        how='outer'
    )

    # Process Seqspecs
    if 'Seqspecs' in merged_df.columns:
        merged_df['Seqspecs'] = merged_df['Seqspecs'].astype(str).replace('nan', '')
        merged_df.loc[merged_df['Seqspecs'] != '', 'Seqspecs'] = \
            merged_df[merged_df['Seqspecs'] != '']['Seqspecs'].str.extract(r'/configuration-files/([^/]+)/', expand=False)

    merged_df = merged_df.sort_values(['measurement_sets', 'Lane', 'file_modality', 'Flowcell ID'])
    
    final_cols = [
        'R1_path', 'R1_md5sum', 'R2_path', 'R2_md5sum', 
        'measurement_sets', 'Lane', 'file_modality', 
        'Flowcell ID', 'Seqspecs'
    ]
    
    return merged_df[final_cols].reset_index(drop=True)

def main():
    parser = argparse.ArgumentParser(description='Process IGVF analysis sets to sample files')
    
    parser.add_argument('-i', '--analysis-set-id', required=True, help='Analysis set ID')
    parser.add_argument('--access-key', required=True, help='IGVF API access key')
    parser.add_argument('--secret-key', required=True, help='IGVF API secret key')
    parser.add_argument('-per-sample-output', default='per_sample_file.tsv', help='Output file for per-sample data')
    parser.add_argument('-analysis-set-output', default='analysis_sets.tsv', help='Output file for analysis set data')
    
    args = parser.parse_args()

    try:
        # Step 1: Fetch analysis sets
        analysis_sets = fetch_analysis_sets(
            args.analysis_set_id,
            args.access_key,
            args.secret_key
        )
        
        # Save analysis sets data
        analysis_sets.to_csv(args.analysis_set_output, index=False, sep='\t')
        print(f"Successfully saved analysis sets data to {args.analysis_set_output}")
        
        # Step 2: Process sequence data
        sequence_data = get_sequence_files(analysis_sets, args.access_key, args.secret_key)
        
        # Step 3: Rearrange sequence files
        sample_file = rearrange_sequence_files(sequence_data)
        
        # Save per-sample data
        sample_file.to_csv(args.per_sample_output, index=False, sep='\t')
        print(f"Successfully saved per-sample data to {args.per_sample_output}")
        
    except Exception as e:
        print(f"Error processing data: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()