# Test data

A tiny dataset for `-profile test`: four libraries (ancient and modern, paired- and
single-end) of 3,000 molecules each, a Kraken2/Bracken database of eight genera, and a
host genome. Every step of the workflow runs on it for real in about a minute.

All sequences are synthetic. The genomes and the host are random sequences; taxon IDs
and names are real NCBI ones so the outputs read naturally, but no sequence comes from
those organisms. Each library mixes known shares of host, genera and unclassifiable
noise, so the sample summary can be checked against what went in.

To regenerate everything (needs Docker):

```bash
python3 make_test_data.py
./build_test_db.sh
```
