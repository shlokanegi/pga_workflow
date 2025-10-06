def get_final_targets(wildcards):
    outputs = []
    
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/hg002.chunked.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    #--- hifiasm full assembly
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap1.p_ctg.gfa", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap2.p_ctg.gfa", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap1.p_ctg.fasta", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap2.p_ctg.fasta", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}.hifiasm.fasta", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_hifiasm_to_hg002_minimap_{asm_preset}.bam", sample_id=SAMPLE_IDS, k=K, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_hifiasm_to_hg002_minimap_{asm_preset}.bam.bai", sample_id=SAMPLE_IDS, k=K, asm_preset=ASM_PRESET))
    
    #--- hifiasm r_utg assembly ----#
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.r_utg.fasta", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_r_utg_to_hg002_minimap_{asm_preset}.bam", sample_id=SAMPLE_IDS, k=K, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_r_utg_to_hg002_minimap_{asm_preset}.bam.bai", sample_id=SAMPLE_IDS, k=K, asm_preset=ASM_PRESET))

    #--- kff file for haplotype sampling
    outputs.extend(expand("results_hs/hs-{k}/reads/{sample_id}.ec.fasta", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/hs-{k}/reads/{sample_id}.ec.kff", sample_id=SAMPLE_IDS, k=K))
    
    # haplotype sampled graph, indexes, aligned reads with LRS, snarls, and converted to pg.vg
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.ec.gbz", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.ec.gaf.gz", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.ec.pg.vg", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.ec.snarls", sample_id=SAMPLE_IDS, k=K))
        
    # chunking reads and generatingsubgraph
    if config.get("RUN_GBZ_QUERY"):
        # GBZ DB outputs
        outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.ec.chains", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/graph/{sample_id}/gbz_db/{sample_id}-{k}-sampled.ec.gbz.db", sample_id=SAMPLE_IDS, k=K))
        # GAF DB outputs
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.ec.sorted.gbwt", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.ec.sorted.gaf.gz", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.ec.sorted.gaf.db", sample_id=SAMPLE_IDS, k=K))
        # query outputs
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/query/subgraph.gaf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/query/subgraph.gfa", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/query/subgraph.pg.vg", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/query/subgraph.pg.dist", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        # anchor generation (with gbz query subgraph or full graph. Depend on config.get("USE_FULL_GRAPH"))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/anchors/subgraph.pkl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/anchors/subgraph.anchors.json.extended.jsonl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/anchors/params_run.log", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        if config.get("RUN_DEBUGGING"):
            outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/anchors/subgraph.anchors.json.reads_processed.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    else:
        # vg chunk outputs
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/ec/alignments-combined.processed.sorted.gaf.gz", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/chunk/subgraph.pg.vg", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/chunk/subgraph.pg.gfa", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/chunk/subgraph.gaf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        #--- annotate with refpos ---#
        # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/chunk/subgraph.gam", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/chunk/subgraph.refpos.gam", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/chunk/subgraph.refpos.gaf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        
        # Anchor generation with vg chunk index and full graph
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/chunk/subgraph.pg.dist", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/anchors/subgraph.pkl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/anchors/subgraph.anchors.json.extended.jsonl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        if config.get("RUN_DEBUGGING"):
            outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/anchors/subgraph.anchors.json.reads_processed.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/anchors/params_run.log", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    # shasta assembly
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/shasta/{sample_id}.subregion.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/shasta/ShastaRun/Assembly.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    # anchor stats
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/anchors/extended_anchor_reads_info.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/extended_anchor_stats", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    if config.get("RUN_DEBUGGING"):
        # reliable snarl stats
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/reliable_snarl_stats", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        #--- debugging files
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/debugging/{region_id}_nodes_info.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/debugging/{region_id}_read_traversals.zip", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/debugging/{region_id}_snarls.bandage.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    # #--- displayPaf
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/assembly_alignment/{asm_preset}_plots", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}_displayPaf.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    # assembly alignment
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.paf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    #--- hifiasm
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.coords.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    # #---- hifiasm hap1 and hap2 assemblies ----#
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.subregion.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/hifiasm_assembly/{sample_id}_{region_id}_hifiasm_subregion_to_hg002_minimap_{asm_preset}.paf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/hifiasm_assembly/{sample_id}_{region_id}_hifiasm_subregion_to_hg002_minimap_{asm_preset}.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    #--- hifiasm r_utg assembly ----#
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.r_utg.subregion.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/hifiasm_assembly/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}.paf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/hifiasm_assembly/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET)) 
    #outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/hifiasm_alignment/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}_displayPaf.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    #--- shasta to hifiasm alignment
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/shasta_to_hifiasm_alignment/{sample_id}_{region_id}_shasta_to_hifiasm_minimap_{asm_preset}.paf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/shasta_to_hifiasm_alignment/{sample_id}_{region_id}_shasta_to_hifiasm_{asm_preset}_alignment_plots.pdf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))

    #---- workflow runlog ----#
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/ec/pga_run.log", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    return outputs
