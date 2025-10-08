# --------- Rules for running the positive control part of the workflow --------- #

rule prepare_positive_control:
    output:
        sampled_hg2_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.gbz",
        sampled_pg_vg_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.pg.vg",
        snarls_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.snarls",
        distance_index_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.dist",
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
    threads: 32
    resources:
        mem_mb=160000,
        runtime=800*60,
        slurm_partition=choose_partition(800)
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

if config.get("RUN_GBZ_QUERY"):
    ### Use GBZ base for DB construction and querying ###
    
    rule extract_top_level_chains_for_positive_control_graph:
        output:
            chains="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.chains"
        input:
            gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.gbz",
            distance_index="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.dist"
        benchmark: "benchmarks/{sample_id}/hs-{k}/extract_top_level_chains_for_positive_control_graph.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/extract_top_level_chains_for_positive_control_graph.log"
        resources:
            mem_mb=13000,
            runtime=3*60,
            slurm_partition=choose_partition(3)
        shell:
            """
            vg chains -p -o {output.chains} {input.gbz} {input.distance_index} > {log} 2>&1
            """
    
    rule generate_gbz_db_for_positive_control:
        output:
            gbz_db="results_hs/graph/{sample_id}/gbz_db/{sample_id}-{k}-sampled.hg2.gbz.db"
        input:
            gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.gbz",
            chains="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.chains"
        benchmark: "benchmarks/{sample_id}/hs-{k}/generate_gbz_db_for_positive_control.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/generate_gbz_db_for_positive_control.log"
        resources:
            mem_mb=7000,
            runtime=420,
            slurm_partition=choose_partition(7)
        shell:
            """
            gbz2db --chains {input.chains} --output {output.gbz_db} {input.gbz} 2> {log}
            """

    rule generate_gaf_db_for_positive_control:
        output:
            gbwt="results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.hg2.sorted.gbwt",
            gaf="results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.hg2.sorted.gaf.gz",
            gaf_db="results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.hg2.sorted.gaf.db"
        input:
            combined_gaf="results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.hg2.gaf"
        benchmark: "benchmarks/{sample_id}/hs-{k}/generate_gaf_db_for_positive_control.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/generate_gaf_db_for_positive_control.log"
        threads: 2
        resources:
            mem_mb=210000,
            runtime=15000,
            slurm_partition=choose_partition(250)
        shell:
            """
            set -e -o pipefail
            echo "Sorting GAF and creating GBWT index..."
            vg gamsort --progress --threads {threads} --gbwt-output {output.gbwt} --gaf-input {input.combined_gaf} | bgzip --threads 16 > {output.gaf}
            
            echo "Creating GAF database..."
            gaf2db -o {output.gaf_db} -g {output.gbwt} --block-size 10 {output.gaf}
            """
        
    ### QUERY USING GBZ DB
    rule query_gbz_gaf_db_for_positive_control:
        output:
            subgraph_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.gaf",
            subgraph_gfa="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.gfa",
            subgraph_pg_vg="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.pg.vg",
            subgraph_pg_dist="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.pg.dist"
        input:
            gbz_db="results_hs/graph/{sample_id}/gbz_db/{sample_id}-{k}-sampled.hg2.gbz.db",
            gaf_db="results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.hg2.sorted.gaf.db"
        params:
            contig=config['region']['chromosome'],
            interval=config['region']['start'] + ".." + config['region']['end'],
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/query_gbz_gaf_db_for_positive_control.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/{region_id}/query_gbz_gaf_db_for_positive_control.log"
        threads: 16
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=500),
            runtime=20*60,
            slurm_partition=choose_partition(20)
        run:
            import subprocess
            # query the given interval to get the reads and subgraph
            query_cmd = (
                f"query --sample CHM13 --contig {params.contig} --context 0 --snarls --interval {params.interval} "
                f"--gaf-base {input.gaf_db} --gaf-output {output.subgraph_gaf} {input.gbz_db} > {output.subgraph_gfa} 2>> {log}"
            )
            subprocess.run(query_cmd, shell=True, check=True)
            subprocess.run(f"vg convert -g {output.subgraph_gfa} -p > {output.subgraph_pg_vg}", shell=True, check=True)        
            subprocess.run(f"vg index {output.subgraph_pg_vg} --threads {threads} --dist-name {output.subgraph_pg_dist}", shell=True, check=True)

    
    if not config.get("USE_FULL_GRAPH"):

        rule generate_anchors_dictionary_with_subgraph_for_positive_control:
            output:
                anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.pkl"
            input:
                vg_anchors_config=config["VG_ANCHORS"]["conf"],
                subgraph_pg_dist="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.pg.dist",
                subgraph_pg_vg="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.pg.vg"
            benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary_with_subgraph_for_positive_control.benchmark.txt"
            log: "logs/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary_with_subgraph_for_positive_control.log"
            resources:
                mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
                runtime=10*60,
                slurm_partition=choose_partition(10)
            shell:
                """
                echo "Generating anchors dictionary with subgraph from GBZ query..."
                vg-anchors --config {input.vg_anchors_config} build --graph {input.subgraph_pg_vg} --index {input.subgraph_pg_dist} \
                    --output-prefix results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph > {log} 2>&1
                """

        rule get_anchors_from_gaf_with_subgraph_for_positive_control:
            output:
                anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl",
                params_log="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/params_run.log",
                **({} if not config.get("RUN_DEBUGGING") else {"read_processed_tsv": "results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.reads_processed.tsv"})
            input:
                vg_anchors_config=config["VG_ANCHORS"]["conf"],
                anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.pkl",
                chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.gaf",
                subgraph_pg_vg="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.pg.vg",
                chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta"
            params:
                processes=32
            benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf_with_subgraph_for_positive_control.benchmark.txt"
            log: "logs/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf_with_subgraph_for_positive_control.log"
            resources:
                mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
                runtime=30*60,
                slurm_partition=choose_partition(30)
            shell:
                """
                echo "Getting anchors from GAF with subgraph from GBZ query..."
                vg-anchors --config {input.vg_anchors_config} get-anchors \
                    --dictionary {input.anchors_dictionary} \
                    --threads {params.processes} \
                    --graph {input.subgraph_pg_vg} \
                    --alignment {input.chunked_gaf} \
                    --fasta {input.chunked_fasta} \
                    --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json > {log} 2>&1
                """
    
    else:
        rule generate_anchors_dictionary_with_full_graph_and_gbz_query_for_positive_control:
            output:
                anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.pkl"
            input:
                vg_anchors_config=config["VG_ANCHORS"]["conf"],
                subgraph_pg_dist="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.pg.dist",
                sampled_pg_vg_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.pg.vg",
            benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary_with_full_graph_and_gbz_query_for_positive_control.benchmark.txt"
            log: "logs/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary_with_full_graph_and_gbz_query_for_positive_control.log"
            resources:
                mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
                runtime=100*60,
                slurm_partition=choose_partition(100)
            shell:
                """
                echo "Generating anchors dictionary with full graph and GBZ query..."
                vg-anchors --config {input.vg_anchors_config} build --graph {input.sampled_pg_vg_hg2} --index {input.subgraph_pg_dist} \
                    --output-prefix results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph > {log} 2>&1
                """

        rule get_anchors_from_gaf_with_full_graph_and_gbz_query_for_positive_control:
            output:
                anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl",
                params_log="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/params_run.log",
                **({} if not config.get("RUN_DEBUGGING") else {"read_processed_tsv": "results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.reads_processed.tsv"})
            input:
                anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.pkl",
                chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.gaf",
                sampled_pg_vg_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.pg.vg",
                chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta"
            params:
                processes=32
            benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf_with_full_graph_and_gbz_query_for_positive_control.benchmark.txt"
            log: "logs/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf_with_full_graph_and_gbz_query_for_positive_control.log"
            resources:
                mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
                runtime=100*60,
                slurm_partition=choose_partition(100)
            shell:
                """
                echo "Getting anchors from GAF with full graph and GBZ query..."
                vg-anchors --config {input.vg_anchors_config} get-anchors \
                    --dictionary {input.anchors_dictionary} \
                    --threads {params.processes} \
                    --graph {input.sampled_pg_vg_hg2} \
                    --alignment {input.chunked_gaf} \
                    --fasta {input.chunked_fasta} \
                    --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json > {log} 2>&1
                """

else:
    rule vg_chunk_and_index_positive_control:
        output:
            subgraph_vg="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.vg",
            chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.gaf",
            subgraph_pg_vg="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.pg.vg",
            subgraph_pg_gfa="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.pg.gfa",
            subgraph_pg_dist="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.pg.dist"
        input:
            sorted_gaf_hg2="results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.hg2.gaf.gz",
            sampled_gbz_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.gbz",
            snarls_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.snarls"
        params:
            region="CHM13#0#" + config['region']['chromosome'] + ":" + config['region']['start'] + "-" + config['region']['end'],
            region_underscore="_0_CHM13#0#" + config['region']['chromosome'] + "_" + config['region']['start'] + "_" + config['region']['end'],
            prefix="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph",
            region_id=config["region_id"]
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/vg_chunk_and_index_positive_control.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/{region_id}/vg_chunk_and_index_positive_control.log"
        threads: 8
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
            runtime=100*60,
            slurm_partition=choose_partition(100)
        shell:
            """
            vg chunk -a {input.sorted_gaf_hg2} -F -g -x {input.sampled_gbz_hg2} -p {params.region} -S {input.snarls_hg2} --trace -t {threads} -b {params.prefix} > {output.subgraph_vg} 2> {log}
            mv {params.prefix}{params.region_underscore}.gaf {output.chunked_gaf_hg2}
            vg convert -p {output.subgraph_vg} > {output.subgraph_pg_vg}
            vg convert -f {output.subgraph_pg_vg} > {output.subgraph_pg_gfa}
            vg index -t {threads} {output.subgraph_pg_vg} --dist-name {output.subgraph_pg_dist}
            """

    rule annotate_gam_gaf_with_refpos_for_positive_control:
        output:
            chunked_gaf_uniq="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.uniq.gaf",
            chunked_gam="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.gam",
            annotated_gam="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.refpos.gam",
            annotated_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.refpos.gaf"
        input:
            script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/add_refpos_gam_to_gaf.py",
            chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.gaf",
            sampled_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.gbz"
        threads: 16
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=300),
            runtime=20*60,
            slurm_partition=choose_partition(20)
        shell:
            """
            sort {input.chunked_gaf} | uniq > {output.chunked_gaf_uniq}
            vg convert --gaf-to-gam {output.chunked_gaf_uniq} {input.sampled_gbz} > {output.chunked_gam}
            vg annotate -x {input.sampled_gbz} -a {output.chunked_gam} -p -P -t {threads} > {output.annotated_gam}
            vg view -a {output.annotated_gam} -j | python {input.script} --gaf {output.chunked_gaf_uniq} --output {output.annotated_gaf}
            """

    rule run_anchor_generatation_with_full_graph_and_vg_chunk_for_positive_control:
        output:
            anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.pkl",
            anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl",
            params_log="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/params_run.log",
            **({} if not config.get("RUN_DEBUGGING") else {"read_processed_tsv": "results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.reads_processed.tsv"})
        input:
            vg_anchors_config=config["VG_ANCHORS"]["conf"],
            sampled_pg_vg_hg2="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.hg2.pg.vg",
            subgraph_pg_dist="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.pg.dist",
            chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.gaf",
            chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta"
        params:
            processes=32
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/run_anchor_generatation_with_full_graph_and_vg_chunk_for_positive_control.benchmark.txt"
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.5, min_mb=300),
            runtime=120*60,
            slurm_partition=choose_partition(120)
        shell:
            """
            vg-anchors --config {input.vg_anchors_config} build --graph {input.sampled_pg_vg_hg2} --index {input.subgraph_pg_dist} \
                --output-prefix results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph

            vg-anchors --config {input.vg_anchors_config} get-anchors --dictionary {output.anchors_dictionary} --graph {input.sampled_pg_vg_hg2} --alignment {input.chunked_gaf} --fasta {input.chunked_fasta}
                --threads {params.processes} --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/anchors/subgraph.anchors.json
            """


rule chunk_fasta_for_positive_control:
    output:
        chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta"
    input:
        fasta="results/reads/{sample_id}.fasta",
        chunked_gaf=(
            "results_hs/hs-{k}/{sample_id}/{region_id}/pc/query/subgraph.gaf"
            if config.get("RUN_GBZ_QUERY")
            else "results_hs/hs-{k}/{sample_id}/{region_id}/pc/chunk/subgraph.gaf"
        )
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/chunk_fasta_for_positive_control.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/chunk_fasta_for_positive_control.log"
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
		shasta_conf=config["SHASTA"]["conf"],
		chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta",
		anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl"
	benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/run_shasta_assembly_for_positive_control.benchmark.txt"
	log: "logs/{sample_id}/hs-{k}/{region_id}/run_shasta_assembly_for_positive_control.log"
	resources:
		mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.5, min_mb=10000),
		runtime=60*60,
		slurm_partition=choose_partition(60)
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
        anchor_reads_info="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/extended_anchor_reads_info.tsv",
        anchor_stats_dir=directory("results_hs/hs-{k}/{sample_id}/{region_id}/pc/extended_anchor_stats")
    input:
        scripts_dir="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts",
        anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl",
        subregion_shasta_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/ShastaRun/Assembly.fasta",
        chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/{sample_id}.subregion.fasta"
    params:
        region=config['region']['chromosome'] + ":" + config['region']['start'] + "-" + config['region']['end'],
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_extended_anchor_stats_for_positive_control.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/get_extended_anchor_stats_for_positive_control.log"
    resources:
        mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=400),
        runtime=10*60,
        slurm_partition=choose_partition(10)
    shell:
        """
        mkdir -p results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/extended_anchor_stats
        # Since this JSON is confusing to interpret, convert it to a TSV...
        python3 {input.scripts_dir}/preprocess_vganchor_outfiles_extention.py -j {input.anchors}
        # Next, generate anchor stats plots...
        Rscript {input.scripts_dir}/get_extended_anchor_stats.R -d results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc -o ont_r10_{wildcards.region_id} -r {params.region} > {log} 2>&1
        """

# #------- Generate anchor sequence TSV --------# (Can plug into rule above if needed)
# master_table_tsv=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/extended_anchor_stats/ont_r10_{wildcards.region_id}_anchors_master_table.tsv
# python {input.scripts_dir}/process_anchor_seqs.py -j {input.anchors} -m $master_table_tsv -f {input.chunked_fasta} -o results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/extended_anchor_stats/ont_r10_{wildcards.region_id}        


if config.get("RUN_DEBUGGING"):

    rule get_reliable_snarl_stats_for_positive_control:
        output:
            reliable_snarl_stats_dir=directory("results_hs/hs-{k}/{sample_id}/{region_id}/pc/reliable_snarl_stats"),
            snarl_compatibility_fractions="results_hs/hs-{k}/{sample_id}/{region_id}/pc/reliable_snarl_stats/snarl_compatibility_fractions.tsv"
        input:
            anchors="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/subgraph.anchors.json.extended.jsonl",
            scripts_dir="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts",
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_reliable_snarl_stats_for_positive_control.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/{region_id}/get_reliable_snarl_stats_for_positive_control.log"
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
            runtime=10*60,
            slurm_partition=choose_partition(10)
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
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_debugging_files_for_positive_control.benchmark.txt"
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=500),
            runtime=15*60,
            slurm_partition=choose_partition(15)
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

rule align_assembly_to_chunked_reference_for_positive_control:
	output:
		paf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.paf",
		csv="results_hs/hs-{k}/{sample_id}/{region_id}/pc/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.csv"
	input:
		analyzePaf_bin=config["ANALYSEPAF"]["bin"],
		hg002_reference_chunked="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/hg002.chunked.fasta",
		subregion_shasta_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/ShastaRun/Assembly.fasta"
	params:
		asm_preset=config["MINIMAP"]["asmPreset"]
	benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/align_assembly_to_chunked_reference_for_positive_control_{asm_preset}.benchmark.txt"
	log: "logs/{sample_id}/hs-{k}/{region_id}/align_assembly_to_chunked_reference_for_positive_control_{asm_preset}.log"
	container: "docker://mkolmogo/card_minimap2:2.23"
	threads: 8
	resources:
		mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=400),
		runtime=10*60,
		slurm_partition=choose_partition(10)
	shell:
		"""
		minimap2 -t {threads} -I 20G -cx {params.asm_preset} -K 1M --eqx --cs {input.hg002_reference_chunked} {input.subregion_shasta_assembly} > {output.paf}
		{input.analyzePaf_bin} \
			--input {output.paf} \
			--output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/assembly_alignment/{wildcards.sample_id}_{wildcards.region_id}_ZOOMED_shasta_to_hg002_minimap_{wildcards.asm_preset} \
			--pixelsPerMb 600 --minAlignmentQuality 1 --minAlignmentLength 0 > {log} 2>&1
		"""

rule run_displayPafAlignments_for_positive_control:
	output:
		csv="results_hs/hs-{k}/{sample_id}/{region_id}/pc/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}_displayPaf.csv"
	input:
		displayPaf_bin=config["DISPLAYPAF"]["bin"],
		hg002_reference_chunked="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/hg002.chunked.fasta",
		subregion_shasta_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/pc/shasta/ShastaRun/Assembly.fasta",
		paf="results_hs/hs-{k}/{sample_id}/{region_id}/pc/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.paf",
	params:
		asm_preset=config["MINIMAP"]["asmPreset"],
	benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/run_displayPafAlignments_for_positive_control_{asm_preset}.benchmark.txt"
	log: "logs/{sample_id}/hs-{k}/{region_id}/run_displayPafAlignments_for_positive_control_{asm_preset}.log"
	resources:
		mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
		runtime=10*60,
		slurm_partition=choose_partition(10)
	shell:
		"""
		{input.displayPaf_bin} \
			--paf {input.paf} -r {input.hg002_reference_chunked} -a {input.subregion_shasta_assembly} \
			--output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/pc/assembly_alignment/{wildcards.sample_id}_{wildcards.region_id}_ZOOMED_shasta_to_hg002_minimap_{wildcards.asm_preset}_displayPaf \
			--minAlignmentQuality 1 --minAlignmentLength 0 > {log} 2>&1
		"""

rule generate_run_summary_positive_control:
    output:
        pga_log="results_hs/hs-{k}/{sample_id}/{region_id}/pc/pga_run.log"
    input:
        params_log="results_hs/hs-{k}/{sample_id}/{region_id}/pc/anchors/params_run.log",
        shasta_conf=config["SHASTA"]["conf"],
        script="workflow/scripts/generate_runlog.py"
    params:
        run_mode=config['RUN_MODE'],
        region_id=config['region_id'],
        asm_preset=config['MINIMAP']['asmPreset'],
        read_type=config.get('READ_TYPE'),
        run_gbz_query=config.get('RUN_GBZ_QUERY'),
        use_full_graph=config.get('USE_FULL_GRAPH'),
        run_debugging=config.get('RUN_DEBUGGING')
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/generate_run_summary_positive_control.benchmark.txt"
    shell:
        """
        python3 {input.script} --params-log {input.params_log} --shasta-conf {input.shasta_conf} --output-log {output.pga_log} --run-mode {params.run_mode} --region-id {params.region_id} --asm-preset {params.asm_preset} \
            --read-type {params.read_type} --run-gbz-query {params.run_gbz_query} --use-full-graph {params.use_full_graph} --run-debugging {params.run_debugging}
        """
