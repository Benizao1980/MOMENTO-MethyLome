REQUIRED_METHYLATION_COLUMNS=["sample_id","platform","motif","modification","called_sites","genomic_sites","methylated_fraction","median_score","source_file","qc_status"]
REQUIRED_METADATA_COLUMNS=["sample_id","species","host","host_group","lineage"]

def missing_columns(columns, required):
    observed=set(columns)
    return [c for c in required if c not in observed]
