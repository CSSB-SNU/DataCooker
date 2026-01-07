# Tutorials

Pipeline tutorials live here. Each recipe under `pipelines/recipe/` has its own short guide covering inputs, targets, and execution steps.

## Library overview
- [A3M recipe](tutorials/a3m_recipe_book.md): parse A3M headers/sequences and build an MSA container.
- [AF3 training LMDB](tutorials/build_af3_training.md): filter CIFMol entries, attach metadata, and remove signal peptides.
- [Build sequence hash map](tutorials/build_seq_hash_map.md): create a FASTA hash map from CIFMol-derived sequences.
- [CCD recipe](tutorials/ccd_recipe_book.md): parse chem comp/atom/bond tables into a clean CCD dict.
- [CIF recipe](tutorials/cif_recipe_book.md): full mmCIF ingestion to assembly/contact graphs plus metadata.
- [Extract FASTA](tutorials/extract_fasta.md): derive FASTA strings from CIFMol.
- [Filter A3M LMDB](tutorials/filter_a3m_lmdb.md): trim residue/chain feature containers for AF3.
- [Graph LMDB](tutorials/graph_lmdb.md): build graph bytes per CIFMol with clustering context.
- [Graph LMDB (attached)](tutorials/graph_lmdb_from_attached.md): same as above for attached CIFMol inputs.
- [Load metadata](tutorials/load_metadata.md): load seq/cluster mappings and SignalP outputs.
- [Sequence clustering](tutorials/seq_cluster.md): split FASTA, cluster proteins/antibodies, and write cluster maps.
- [Train/valid graph split](tutorials/train_valid_graph_split.md): build graphs, split train/valid edges, and compute stats.

Use these in order or jump to the pipeline you need. All pages share the same format (purpose, inputs, outputs, steps, usage) for consistency.
