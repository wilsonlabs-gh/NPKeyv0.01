# Published BiG-SLiCE output replay

Four BGCs from the BiG-SLiCE 2.0.0 Amycolatopsis benchmark in
[Medema et al.'s supplementary data](https://zenodo.org/records/17969757),
reused under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

`query.db.gz` retains their BGC, CDS, HSP, feature and membership values in the
exact query-mode schema from our pinned upstream commit. The source is clustering
output: we selected the threshold-0.4 run and removed the reference-only
`bgc.dataset_id` column. Scores, ranks, BGC IDs and GCF IDs are unchanged.
This is a query-format fixture, not a downloaded query job or a genome prediction.

`source.json` records the source member, checksums, schema URL and original rows.
Cases cover assigned, unassigned and no-feature BGCs. The fixture uses the published
benchmark's families, not our intended MIBiG reference. Tests require no downloads.
