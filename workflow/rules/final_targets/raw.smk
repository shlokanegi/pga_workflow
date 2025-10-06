def get_final_targets(wildcards):
    outputs = []
    # outputs.extend(expand("results/reads/{sample_id}.fasta", sample_id=SAMPLE_IDS))
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
    outputs.extend(expand("results/reads/{sample_id}.kff", sample_id=SAMPLE_IDS))

    # haplotype sampled graph, indexes, aligned reads with LRS, snarls, and converted to pg.vg
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.dist", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.pg.vg", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.snarls", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.longread.withzip.min", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.longread.zipcodes", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/{sample_id}.gaf.zst", sample_id=SAMPLE_IDS, k=K))
    outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.gaf", sample_id=SAMPLE_IDS, k=K))
    
    # chunking reads and generatingsubgraph
    if config.get("RUN_GBZ_QUERY"):
        # GBZ DB outputs
        outputs.extend(expand("results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.chains", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/graph/{sample_id}/gbz_db/{sample_id}-{k}-sampled.gbz.db", sample_id=SAMPLE_IDS, k=K))
        # GAF DB outputs
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.sorted.gbwt", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.sorted.gaf.gz", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.sorted.gaf.db", sample_id=SAMPLE_IDS, k=K))
        # query outputs
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.gaf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.gfa", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.pg.vg", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.pg.dist", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        
        # anchor generation (with gbz query subgraph or full graph. Depend on config.get("USE_FULL_GRAPH"))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.pkl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.extended.jsonl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/anchors/params_run.log", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))              
        if config.get("RUN_DEBUGGING"):
            outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.reads_processed.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    else:
        # vg chunk outputs
        outputs.extend(expand("results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.sorted.gaf.gz", sample_id=SAMPLE_IDS, k=K))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.pg.vg", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.pg.gfa", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.gaf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        #--- annotate with refpos ---#
        # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.gam", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.refpos.gam", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.refpos.gaf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
       
       # Anchor generation with vg chunk index and full graph
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.pg.dist", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.pkl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.extended.jsonl", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        if config.get("RUN_DEBUGGING"):
            outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.reads_processed.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/anchors/params_run.log", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    # shasta assembly
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/shasta/{sample_id}.subregion.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/shasta/ShastaRun/Assembly.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    # anchor stats
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/anchors/extended_anchor_reads_info.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/extended_anchor_stats", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    if config.get("RUN_DEBUGGING"):
        # reliable snarl stats
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/reliable_snarl_stats", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        #--- debugging files
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/debugging/{region_id}_nodes_info.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/debugging/{region_id}_read_traversals.zip", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
        outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/debugging/{region_id}_snarls.bandage.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))
    
    # #--- displayPaf
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}_displayPaf.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    # assembly alignment
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.paf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    #--- hifiasm
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.coords.tsv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    # #---- hifiasm hap1 and hap2 assemblies ----#
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.subregion.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_hifiasm_subregion_to_hg002_minimap_{asm_preset}.paf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_hifiasm_subregion_to_hg002_minimap_{asm_preset}.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    #--- hifiasm r_utg assembly ----#
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.r_utg.subregion.fasta", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}.paf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET)) 
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_alignment/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}_displayPaf.csv", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    
    #--- shasta to hifiasm alignment
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/shasta_to_hifiasm_alignment/{sample_id}_{region_id}_shasta_to_hifiasm_minimap_{asm_preset}.paf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))
    # outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/shasta_to_hifiasm_alignment/{sample_id}_{region_id}_shasta_to_hifiasm_{asm_preset}_alignment_plots.pdf", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID, asm_preset=ASM_PRESET))

    #---- workflow runlog ----#
    outputs.extend(expand("results_hs/hs-{k}/{sample_id}/{region_id}/pga_run.log", sample_id=SAMPLE_IDS, k=K, region_id=REGION_ID))

    return outputs
