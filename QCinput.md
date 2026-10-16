# Input requirements

- One bacterial genome per run; upstream genome QC already completed.
- A local, uncompressed annotated GenBank file (.gb, .gbk or .gbff), or a versioned NCBI assembly accession with annotated GenBank available.
- Supplied CDS annotations are required; this pipeline does not perform gene finding.
- At least one CDS annotation across the genome; individual contigs without CDS are retained.
- Non-empty DNA sequences using valid IUPAC DNA letters.
- Unique sequence identifiers without whitespace.
- CDS locations must be readable, within the sequence bounds and refer to the same sequence record.
- Include all assembly records, including chromosomes and plasmids.
- Raw FASTA alone is not supported.
