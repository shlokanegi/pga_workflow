
def get_final_targets(wildcards):
    outputs = []
    
    # haplotype sampled graph, indexes, aligned reads with LRS, snarls, and converted to pg.vg
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.gbz", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.dist", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.pg.vg", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.snarls", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.hg2.gaf.gz", sample_id=SAMPLE_IDS, k=K))
    
    # chunking reads and generatingsubgraph
    if config.get("RUN_GBZ_QUERY"):
        # GBZ DB outputs
        outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.chains", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/graph/{sample_id}/gbz_db/{sample_id}-{k}-sampled.hg2.gbz.db", sample_id=SAMPLE_IDS, k=K))
        # GAF DB outputs
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.hg2.sorted.gbwt", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.hg2.sorted.gaf.gz", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.hg2.sorted.gaf.db", sample_id=SAMPLE_IDS, k=K))
        # query outputs
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.gaf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.gfa", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.pg.vg", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.pg.dist", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        
        # anchor generation (with gbz query subgraph or full graph. Depend on config.get("USE_FULL_GRAPH"))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.pkl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/params_run.log", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))              
        if config.get("RUN_DEBUGGING"):
            outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.reads_processed.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    else:
        # vg chunk outputs
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.pg.vg", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.pg.gfa", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.gaf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        
        # Anchor generation with vg chunk index and full graph
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.pg.dist", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.pkl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        if config.get("RUN_DEBUGGING"):
            outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.reads_processed.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/params_run.log", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    # shasta assembly
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/ShastaRun/Assembly.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    # anchor stats
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/extended_anchor_reads_info.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/extended_anchor_stats", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    if config.get("RUN_DEBUGGING"):
        # reliable snarl stats
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/reliable_snarl_stats", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        #--- debugging files
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/debugging/{region_id}_nodes_info.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/debugging/{region_id}_read_traversals.zip", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/debugging/{region_id}_snarls.bandage.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    # assembly alignment
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.paf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    #---- workflow runlog ----#
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pc/pga_run.log", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))

    return outputs
