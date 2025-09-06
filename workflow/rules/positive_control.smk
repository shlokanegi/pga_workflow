# --------- Rules for running the positive control part of the workflow --------- #

rule prepare_positive_control:
    output:
        sampled_hg2_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.gbz",
        sampled_pg_vg_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.pg.vg",
        snarls_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.snarls",
        alignment_gaf="results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.hg2.gaf",
        alignment_gaf_sorted="results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.hg2.gaf.gz"
    input:
        kff="results/reads/{sample_id}.kff",
        graph_gbz_hg2=config["graph_base_hg2"] + ".gbz",
        graph_hapl_hg2=config["graph_base_hg2"] + ".hapl",
        awk_script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/process_out.awk",
        ont_r10_sequence="/private/groups/migalab/kkyriaki/experiments/data/GIAB_2025/HG002/{sample_id}/{sample_id}.fastq"
    params:
        k=config["HAPLOTYPE_SAMPLING"]["num_haps"],
        diploid_sampling="--diploid_sampling" if config["HAPLOTYPE_SAMPLING"]["diploid_sampling"] else ""
    benchmark: "benchmarks/{sample_id}/hs-{k}/prepare_positive_control.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/prepare_positive_control.log"
    threads: 128
    shell:
        """
        echo "------Haplotype sampling of HG2 graph------"
        vg haplotypes -v 2 -t {threads} --include-reference {params.diploid_sampling} --num-haplotypes {params.k} \
            -i {input.graph_hapl_hg2} -k {input.kff} -g {output.sampled_hg2_gbz} {input.graph_gbz_hg2}

        echo "------Convert the sampled HG2 graph to a PGVG graph------"
        vg convert -t {threads} {output.sampled_hg2_gbz} -p > {output.sampled_pg_vg_hg2}

        echo "------Compute snarls for the sampled HG2 graph------"
        vg snarls -t {threads} -T {output.sampled_hg2_gbz} > {output.snarls_hg2}

        echo "------Align reads to the sampled HG2 graph------"
        vg giraffe -t {threads} -b r10 -f {input.ont_r10_sequence} -Z {output.sampled_hg2_gbz} --output-format gaf 2>> {log} | {input.awk_script} > {output.alignment_gaf}

        echo "------Sort and index the alignment file------"
        vg gamsort -t {threads} -p -G {output.alignment_gaf} | bgzip -c > {output.alignment_gaf_sorted}
        tabix -@ {threads} -p gaf {output.alignment_gaf_sorted}
        """

rule vg_chunk_and_index_positive_control:
	output:
		subgraph_vg_hg2="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.vg",
		chunked_gaf_hg2="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.gaf",
		subgraph_pg_vg_hg2="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.pg.vg",
		subgraph_pg_gfa_hg2="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.pg.gfa",
		subgraph_pg_dist_hg2="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.pg.dist"
	input:
		sorted_gaf_hg2="results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.hg2.gaf.gz",
		sampled_gbz_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.gbz",
		snarls_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.snarls"
	params:
		region="CHM13#0#" + config['region']['chromosome'] + ":" + config['region']['start'] + "-" + config['region']['end'],
		region_underscore="_0_CHM13#0#" + config['region']['chromosome'] + "_" + config['region']['start'] + "_" + config['region']['end'],
		prefix="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph",
		region_id=config["region_id"]
	benchmark:
		"benchmarks/{sample_id}/hs-{k}/{region_id}/vg_chunk_and_index_positive_control.benchmark.txt"
	log:
		"logs/{sample_id}/hs-{k}/{region_id}/vg_chunk_and_index_positive_control.log"
	threads: 128
	shell:
		"""
		vg chunk -a {input.sorted_gaf_hg2} -F -g -x {input.sampled_gbz_hg2} -p {params.region} -S {input.snarls_hg2} --trace -t {threads} -b {params.prefix} > {output.subgraph_vg_hg2} 2> {log}
		mv {params.prefix}{params.region_underscore}.gaf {output.chunked_gaf_hg2}
		vg convert -p {output.subgraph_vg_hg2} > {output.subgraph_pg_vg_hg2}
		vg convert -f {output.subgraph_pg_vg_hg2} > {output.subgraph_pg_gfa_hg2}
		vg index -t {threads} {output.subgraph_pg_vg_hg2} --dist-name {output.subgraph_pg_dist_hg2}
		"""

rule annotate_gam_gaf_with_refpos_for_positive_control:
    output:
        chunked_gaf_uniq="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.uniq.gaf",
        chunked_gam="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.gam",
        annotated_gam="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.refpos.gam",
        annotated_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.refpos.gaf"
    input:
        script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/add_refpos_gam_to_gaf.py",
        chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.gaf",
        sampled_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.gbz"
    threads: 16
    shell:
        """
        sort {input.chunked_gaf} | uniq > {output.chunked_gaf_uniq}
        vg convert --gaf-to-gam {output.chunked_gaf_uniq} {input.sampled_gbz} > {output.chunked_gam}
        vg annotate -x {input.sampled_gbz} -a {output.chunked_gam} -p -P -t {threads} > {output.annotated_gam}
        vg view -a {output.annotated_gam} -j | python {input.script} --gaf {output.chunked_gaf_uniq} --output {output.annotated_gaf}
        """

rule run_anchor_generatation_for_positive_control:
    output:
        anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.pkl",
        anchors_preext="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.jsonl",
        anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl",
        read_processed_tsv="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.reads_processed.tsv",
        pruned_anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.pruned.jsonl",
        params_log="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/params_run.log"
    input:
        sampled_pg_vg_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.pg.vg",
        subgraph_pg_dist_hg2="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.pg.dist",
        chunked_gaf_hg2="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.gaf",
        chunked_fasta_hg2="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta"
    benchmark:
        "benchmarks/{sample_id}/hs-{k}/{region_id}/run_anchor_generatation_for_positive_control.benchmark.txt"
    shell:
        """
        vg_anchor build --graph {input.sampled_pg_vg_hg2} --index {input.subgraph_pg_dist_hg2} \
            --output-prefix results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph

        vg_anchor get-anchors --dictionary {output.anchors_dictionary} --graph {input.sampled_pg_vg_hg2} --alignment {input.chunked_gaf_hg2} --fasta {input.chunked_fasta_hg2}
            --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json
        """


rule chunk_fasta_for_positive_control:
    output:
        chunked_fasta = "results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta"
    input:
        chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.hg2.gaf",
        fasta="results/reads/{sample_id}.fasta"
    container:
        "docker://pegi3s/seqkit:latest"
    shell:
        """
        # Extract READ IDs from chunked GAF and then extract fasta records for these reads
        cut -f1 {input.chunked_gaf} | sort -u -T {config[TMPDIR]}| seqkit grep -f - {input.fasta} > {output.chunked_fasta}
        """

rule run_shasta_assembly_for_positive_control:
    output: "results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/ShastaRun/Assembly.fasta"
    input:
        shasta=config["SHASTA"]["bin"],
        chunked_fasta = "results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta",
        shasta_conf = config["SHASTA"]["conf"],
        anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl",
        # anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.pruned.jsonl"
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/run_shasta_assembly_for_positive_control.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/run_shasta_assembly_for_positive_control.log"
    shell:
        """
        ## make a copy of the anchors.json for Shasta
        cp {input.anchors} results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/shasta/anchors.json
        rm -rf results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/shasta/ShastaRun

        {input.shasta} --input {input.chunked_fasta} --assemblyDirectory results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/shasta/ShastaRun \
        	--config {input.shasta_conf} --anchors results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/shasta/anchors.json --memoryMode filesystem --memoryBacking disk > {log} 2>&1
        """

rule get_extended_anchor_stats_for_positive_control:
    output:
        anchor_reads_info="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/pc/extended_anchor_reads_info.tsv",
        anchor_stats_dir=directory("results_hs/hs-{k}/{sample_id}/{region_id}/pc/extended_anchor_stats")
    input:
        scripts_dir="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts",
        anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl",
        # anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.pruned.jsonl",
        subregion_shasta_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/ShastaRun/Assembly.fasta",
        chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta"
    params:
        region=config['region']['chromosome'] + ":" + config['region']['start'] + "-" + config['region']['end'],
    log: "logs/{sample_id}/hs-{k}/{region_id}/get_extended_anchor_stats_for_positive_control.log"
    shell:
        """
        mkdir -p results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/extended_anchor_stats
        # Since this JSON is confusing to interpret, convert it to a TSV...
        python3 {input.scripts_dir}/preprocess_vganchor_outfiles_extention.py -j {input.anchors}
        # Next, generate anchor stats plots...
        Rscript {input.scripts_dir}/get_extended_anchor_stats.R -d results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc -o ont_r10_{wildcards.region_id} -r {params.region} > {log} 2>&1

        #------- Generate anchor sequence TSV --------#
        master_table_tsv=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/extended_anchor_stats/ont_r10_{wildcards.region_id}_anchors_master_table.tsv
        python {input.scripts_dir}/process_anchor_seqs.py -j {input.anchors} -m $master_table_tsv -f {input.chunked_fasta} -o results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/extended_anchor_stats/ont_r10_{wildcards.region_id}        
        """

rule get_reliable_snarl_stats_for_positive_control:
    output:
        reliable_snarl_stats_dir=directory("results_hs/hs-{k}/{sample_id}/{region_id}/pc/reliable_snarl_stats"),
        snarl_compatibility_fractions="results_hs/hs-{k}/{sample_id}/{region_id}/pc/reliable_snarl_stats/snarl_compatibility_fractions.tsv"
    input:
        anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl",
        scripts_dir="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts",
    log: "logs/{sample_id}/hs-{k}/{region_id}/get_reliable_snarl_stats_for_positive_control.log"
    shell:
        """
        #------- Generate reliable snarl stats --------#
        output_dir=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/reliable_snarl_stats
        reliable_snarls_tsv=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json.reliable_snarls.tsv
        snarl_compatibility_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json.snarl_compatibility.jsonl
        snarl_coverage_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json.snarl_coverage.jsonl
        snarl_allelic_coverage_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json.snarl_allelic_coverage.jsonl
        snarl_allelic_coverage_extended_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json.snarl_allelic_coverage_extended.jsonl
        snarl_variant_type_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json.snarl_variant_type.jsonl
        snarl_positions_tsv=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.sizes.tsv

        Rscript {input.scripts_dir}/reliable_snarl_stats_test.R $reliable_snarls_tsv $snarl_compatibility_json $snarl_coverage_json $snarl_allelic_coverage_json $snarl_allelic_coverage_extended_json $snarl_variant_type_json $snarl_positions_tsv $output_dir >> {log} 2>&1
        """

rule get_debugging_files_for_positive_control:
    output:
        nodes_info_tsv="results_hs/hs-{k}/{sample_id}/{region_id}/pc/debugging/{region_id}_nodes_info.tsv",
        read_traversals_zip="results_hs/hs-{k}/{sample_id}/{region_id}/pc/debugging/{region_id}_read_traversals.zip",
        snarls_bandage_csv="results_hs/hs-{k}/{sample_id}/{region_id}/pc/debugging/{region_id}_snarls.bandage.csv"
    input:
        script_dir="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts",
        read_processed_tsv="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.reads_processed.tsv",
        snarl_compatibility="results_hs/hs-{k}/{sample_id}/{region_id}/pc/reliable_snarl_stats/snarl_compatibility_fractions.tsv"
    shell:
        """
        python {input.script_dir}/process_read_processed_for_bandage.py {input.read_processed_tsv} -o {output.nodes_info_tsv}
        # selected columns for bandage and convert to CSV
        cut -f1,2,5,7 {output.nodes_info_tsv} | tr '\\t' ',' > results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/debugging/{wildcards.region_id}_nodes_info.bandage.csv
        
        # read alignments for bandage
        python {input.script_dir}/extract_read_alignments_for_bandage.py {input.read_processed_tsv} {output.nodes_info_tsv} -o {output.read_traversals_zip}

        # reliable/unreliable snarls for bandage
        snarl_dict="results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.forward_dict.csv"
        python {input.script_dir}/map_snarls_bandage.py $snarl_dict {input.snarl_compatibility} {output.nodes_info_tsv} -o {output.snarls_bandage_csv}
        """


rule align_assembly_to_chunked_reference_and_run_displayPafAlignments_for_positive_control:
    output:
        paf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.paf",
        csv="results_hs/hs-{k}/{sample_id}/{region_id}/pc/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}_displayPaf.csv"    
    input:
        analyzePaf_bin=config["ANALYSEPAF"]["bin"],
        displayPaf_bin=config["DISPLAYPAF"]["bin"],
        hg002_reference_chunked="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/hg002.chunked.fasta",
        subregion_shasta_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/ShastaRun/Assembly.fasta"
    params:
        asm_preset=config["MINIMAP"]["asmPreset"]
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/align_assembly_to_chunked_reference_and_run_displayPafAlignments_for_positive_control_{asm_preset}.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/align_assembly_to_chunked_reference_and_run_displayPafAlignments_for_positive_control_{asm_preset}.log"
    container: "docker://mkolmogo/card_minimap2:2.23"
    threads: 64
    shell:
        """
        minimap2 -t {threads} -I 20G -cx {params.asm_preset} -K 1M --eqx --cs {input.hg002_reference_chunked} {input.subregion_shasta_assembly} > {output.paf}
        {input.analyzePaf_bin} \
            --input {output.paf} \
            --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/assembly_alignment/{wildcards.sample_id}_{wildcards.region_id}_ZOOMED_shasta_to_hg002_minimap_{wildcards.asm_preset} \
            --pixelsPerMb 5000 --minAlignmentQuality 0 --minAlignmentLength 0 > {log} 2>&1

        {input.displayPaf_bin} \
            --paf {output.paf} -r {input.hg002_reference_chunked} -a {input.subregion_shasta_assembly} \
            --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/assembly_alignment/{wildcards.sample_id}_{wildcards.region_id}_ZOOMED_shasta_to_hg002_minimap_{wildcards.asm_preset}_displayPaf \
            --minAlignmentQuality 1 --minAlignmentLength 0 > {log} 2>&1
        """

rule generate_run_summary_positive_control:
    output:
        pga_log="results_hs/hs-{k}/{sample_id}/{region_id}/pc/pga_run.log"
    input:
        params_log="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/params_run.log",
        shasta_conf=config["SHASTA"]["conf"],
        script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/generate_runlog.py"
    params:
        run_mode=config['RUN_MODE'],
        region_id=config['region_id'],
        asm_preset=config['MINIMAP']['asmPreset']
    shell:
        """
        python3 {input.script} --params-log {input.params_log} --shasta-conf {input.shasta_conf} --output-log {output.pga_log} --run-mode {params.run_mode} --region-id {params.region_id} --asm-preset {params.asm_preset}
        """
