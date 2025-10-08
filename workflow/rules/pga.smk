# --------- PGA workflow with raw reads --------- #

# rule convert_fq_to_fasta:
#     output:
#         fasta="results/reads/{sample_id}.fasta"
#     input:
#         ont_r10_sequence="/private/groups/migalab/kkyriaki/experiments/data/GIAB_2025/HG002/{sample_id}/{sample_id}.fastq"
#     benchmark:
#         "benchmarks/{sample_id}/convert_fq_to_fasta.benchmark.txt"
#     log:
#         "logs/{sample_id}/convert_fq_to_fasta.log"
#     shell:
#         """
#         # Convert FASTQ to FASTA
#         cat {input.ont_r10_sequence} | awk '{{if(NR%4==1) {{printf(">%s\\n",substr($0,2));}} else if(NR%4==2) print;}}' > {output.fasta}
#         """

rule count_kmers:
    output:
        kff="results/reads/{sample_id}.kff"
    input:
        fasta="results/reads/{sample_id}.fasta"
    benchmark:
        "benchmarks/{sample_id}/convert_fq_to_fasta.benchmark.txt"
    log:
        "logs/{sample_id}/convert_fq_to_fasta.log"
    threads: 32
    resources:
        mem_mb=125000,  # max_rss = 104393.51 MB (≈ 102 GB) + 20% buffer = 125000 MB
        runtime=2100,  # Runtime in seconds
        slurm_partition=choose_partition(40)
    shell:
        """
        # Count kmers from sequencing reads
        mkdir -p results/reads/tmp_{wildcards.sample_id} 
        kmc -k29 -m128 -okff -t{threads} -fa {input.fasta} results/reads/{wildcards.sample_id} results/reads/tmp_{wildcards.sample_id}
        """

rule haplotype_sample_graph:
    output:
        sampled_gbz = "results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz"
    input:
        kff = "results/reads/{sample_id}.kff",
        graph_gbz = config["graph_base"] + ".gbz",
        graph_hapl = config["graph_base"] + ".hapl"
    params:
        k = config["HAPLOTYPE_SAMPLING"]["num_haps"],
        diploid_sampling = "--diploid_sampling" if config["HAPLOTYPE_SAMPLING"]["diploid_sampling"] else ""
    benchmark:
        "benchmarks/{sample_id}/hs-{k}/haplotype_sample_graph.benchmark.txt"
    log:
        "logs/{sample_id}/hs-{k}/haplotype_sample_graph.log"
    threads: 8
    resources:
        mem_mb=90000,
        runtime=2400,
        slurm_partition=choose_partition(40)
    shell:
        """
        vg haplotypes -v 2 -t {threads} --include-reference {params.diploid_sampling} --num-haplotypes {params.k} \
            -i {input.graph_hapl} -k {input.kff} -g {output.sampled_gbz} {input.graph_gbz} > {log} 2>&1
        """

rule vg_convert_graph:
	output:
		sampled_pg_vg="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.pg.vg"
	input:
		sampled_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz"
	benchmark:
		"benchmarks/{sample_id}/hs-{k}/vg_convert_graph.benchmark.txt"
	threads: 2
	resources:
		mem_mb=35000,
		runtime=4800,
		slurm_partition=choose_partition(80)
	shell:
		"""
		vg convert -t {threads} {input.sampled_gbz} -p > {output.sampled_pg_vg}
		"""

rule compute_snarls:
	output:
		snarls="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.snarls"
	input:
		sampled_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz"
	threads: 32
	shell:
		"""
		vg snarls -t {threads} -T {input.sampled_gbz} > {output.snarls}
		"""

rule generate_sampled_graph_indexes:
	output:
		distance_index="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.dist",
		minimizer_index="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.longread.withzip.min",
		zipcodes="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.longread.zipcodes"
	input:
		sampled_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz"
	benchmark:
		"benchmarks/{sample_id}/hs-{k}/generate_sampled_graph_indexes.benchmark.txt"
	log:
		"logs/{sample_id}/hs-{k}/generate_sampled_graph_indexes.log"
	threads: 8
	resources:
		mem_mb=105000,
		runtime=5400,
		slurm_partition=choose_partition(90)
	shell:
		"""
		vg autoindex -G {input.sampled_gbz} -t {threads} -p results_hs/graph/{wildcards.sample_id}/{wildcards.sample_id}-{wildcards.k}-sampled -w lr-giraffe > {log} 2>&1
		"""

rule align_reads_with_vg_giraffe:
    output:
        alignment="results_hs/alignment/hs-{k}/{sample_id}/{sample_id}.gaf.zst"
    input:
        script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/process_out.awk",
        sampled_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz",
        distance_index="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.dist",
        minimizer_index="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.longread.withzip.min",
        zipcodes="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.longread.zipcodes",
        ont_r10_sequence="/private/groups/migalab/kkyriaki/experiments/data/GIAB_2025/HG002/{sample_id}/{sample_id}.fastq"
    benchmark:
        "benchmarks/{sample_id}/hs-{k}/align_reads_with_vg_giraffe.benchmark.txt"
    log:
        "logs/{sample_id}/hs-{k}/align_reads_with_vg_giraffe.log"
    threads: 32
    resources:
        mem_mb=150000,
        runtime=759*60,
        slurm_partition=choose_partition(759)
    shell:
        """
        vg giraffe -t {threads} --parameter-preset r10 -f {input.ont_r10_sequence} -Z {input.sampled_gbz} -d {input.distance_index} -m {input.minimizer_index} -z {input.zipcodes} --output-format gaf 2> {log} | {input.script} | zstd > {output.alignment}
        """

rule parse_gaf:
    output:
        file="results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.gaf"
    input:
        script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/parse_gaf.py",
        gaf="results_hs/alignment/hs-{k}/{sample_id}/{sample_id}.gaf.zst"
    log:
        "logs/{sample_id}/hs-{k}/parse_gaf.log"
    shell:
        """
        python3 {input.script} {input.gaf} > {output.file} 2> {log}
        """


if config.get("RUN_GBZ_QUERY", False):
    ### Use GBZ base for DB construction and querying ###
    
    rule extract_top_level_chains:
        output:
            chains="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.chains"
        input:
            gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz",
            distance_index="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.dist"
        benchmark: "benchmarks/{sample_id}/hs-{k}/extract_top_level_chains.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/extract_top_level_chains.log"
        resources:
            mem_mb=13000,
            runtime=3*60,
            slurm_partition=choose_partition(3)
        shell:
            """
            vg chains -p -o {output.chains} {input.gbz} {input.distance_index} > {log} 2>&1
            """
    
    rule generate_gbz_db:
        output:
            gbz_db="results_hs/graph/{sample_id}/gbz_db/{sample_id}-{k}-sampled.gbz.db"
        input:
            gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz",
            chains="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.chains"
        benchmark: "benchmarks/{sample_id}/hs-{k}/generate_gbz_db.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/generate_gbz_db.log"
        resources:
            mem_mb=7000,
            runtime=420,
            slurm_partition=choose_partition(7)
        shell:
            """
            gbz2db --chains {input.chains} --output {output.gbz_db} {input.gbz} 2> {log}
            """

    rule generate_gaf_db:
        output:
            gbwt="results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.sorted.gbwt",
            gaf="results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.sorted.gaf.gz",
            gaf_db="results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.sorted.gaf.db"
        input:
            combined_gaf="results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.gaf"
        benchmark: "benchmarks/{sample_id}/hs-{k}/generate_gaf_db.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/generate_gaf_db.log"
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
    rule query_gbz_gaf_db:
        output:
            subgraph_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.gaf",
            subgraph_gfa="results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.gfa",
            subgraph_pg_vg="results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.pg.vg",
            subgraph_pg_dist="results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.pg.dist"
        input:
            gbz_db="results_hs/graph/{sample_id}/gbz_db/{sample_id}-{k}-sampled.gbz.db",
            gaf_db="results_hs/alignment/hs-{k}/{sample_id}/gaf_db/alignments-combined.processed.sorted.gaf.db"
        params:
            contig=config['region']['chromosome'],
            interval=config['region']['start'] + ".." + config['region']['end'],
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/query_gbz_gaf_db.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/{region_id}/query_gbz_gaf_db.log"
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

        rule generate_anchors_dictionary_with_subgraph:
            output:
                anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.pkl"
            input:
                vg_anchors_config=config["VG_ANCHORS"]["conf"],
                subgraph_pg_dist="results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.pg.dist",
                subgraph_pg_vg="results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.pg.vg"
            benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary_with_subgraph.benchmark.txt"
            log: "logs/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary_with_subgraph.log"
            resources:
                mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
                runtime=10*60,
                slurm_partition=choose_partition(10)
            shell:
                """
                vg-anchors --config {input.vg_anchors_config} build --graph {input.subgraph_pg_vg} --index {input.subgraph_pg_dist} \
                    --output-prefix results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph > {log} 2>&1
                """

        rule get_anchors_from_gaf_with_subgraph:
            output:
                anchors="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.extended.jsonl",
                params_log="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/params_run.log",
                **({} if not config.get("RUN_DEBUGGING") else {"read_processed_tsv": "results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.reads_processed.tsv"})
            input:
                vg_anchors_config=config["VG_ANCHORS"]["conf"],
                anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.pkl",
                chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.gaf",
                subgraph_pg_vg="results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.pg.vg",
                chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/shasta/{sample_id}.subregion.fasta"
            params:
                processes=32
            benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf_with_subgraph.benchmark.txt"
            log: "logs/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf_with_subgraph.log"
            resources:
                mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
                runtime=30*60,
                slurm_partition=choose_partition(30)
            shell:
                """
                vg-anchors --config {input.vg_anchors_config} get-anchors \
                    --dictionary {input.anchors_dictionary} \
                    --threads {params.processes} \
                    --graph {input.subgraph_pg_vg} \
                    --alignment {input.chunked_gaf} \
                    --fasta {input.chunked_fasta} \
                    --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.anchors.json > {log} 2>&1
                """
    
    else:
        rule generate_anchors_dictionary_with_full_graph_and_gbz_query:
            output:
                anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.pkl"
            input:
                vg_anchors_config=config["VG_ANCHORS"]["conf"],
                subgraph_pg_dist="results_hs/{k}/{sample_id}/{region_id}/query/subgraph.pg.dist",
                sampled_pg_vg="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.pg.vg",
            benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary_with_full_graph_and_gbz_query.benchmark.txt"
            log: "logs/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary_with_full_graph_and_gbz_query.log"
            resources:
                mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
                runtime=100*60,
                slurm_partition=choose_partition(100)
            shell:
                """
                echo "Generating anchors dictionary with full graph and GBZ query..."
                vg-anchors --config {input.vg_anchors_config} build --graph {input.sampled_pg_vg} --index {input.subgraph_pg_dist} \
                    --output-prefix results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph > {log} 2>&1
                """

        rule get_anchors_from_gaf_with_full_graph_and_gbz_query:
            output:
                anchors="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.extended.jsonl",
                params_log="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/params_run.log",
                **({} if not config.get("RUN_DEBUGGING") else {"read_processed_tsv": "results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.reads_processed.tsv"})
            input:
                vg_anchors_config=config["VG_ANCHORS"]["conf"],
                anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.pkl",
                chunked_gaf="results_hs/{k}/{sample_id}/{region_id}/query/subgraph.gaf",
                sampled_pg_vg="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.pg.vg",
                chunked_fasta="results_hs/{k}/{sample_id}/{region_id}/shasta/{sample_id}.subregion.fasta"
            params:
                processes=32
            benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf_with_full_graph_and_gbz_query.benchmark.txt"
            log: "logs/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf_with_full_graph_and_gbz_query.log"
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
                    --graph {input.sampled_pg_vg} \
                    --alignment {input.chunked_gaf} \
                    --fasta {input.chunked_fasta} \
                    --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.anchors.json > {log} 2>&1
                """

else:
    rule sort_index_gaf:
        output:
            sorted_gaf = "results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.sorted.gaf.gz",
            gaf_index = "results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.sorted.gaf.gz.tbi"
        input:
            combined_gaf="results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.gaf"
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/sort_index_gaf.benchmark.txt"
        threads: 64
        shell:
            """
            vg gamsort -t {threads} -p -G {input.combined_gaf} | bgzip -c > {output.sorted_gaf}
            tabix -@ {threads} -p gaf {output.sorted_gaf}
            """

    rule vg_chunk:
        output:
            subgraph_vg="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.vg",
            chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.gaf",
            subgraph_pg_vg="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.pg.vg",
            subgraph_pg_gfa="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.pg.gfa"
        input:
            sorted_gaf="results_hs/alignment/hs-{k}/{sample_id}/alignments-combined.processed.sorted.gaf.gz",
            sampled_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz",
            snarls="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.snarls"
        params:
            region="CHM13#0#" + config['region']['chromosome'] + ":" + config['region']['start'] + "-" + config['region']['end'],
            region_underscore="_0_CHM13#0#" + config['region']['chromosome'] + "_" + config['region']['start'] + "_" + config['region']['end'],
            prefix="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph",
            region_id=config["region_id"]
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/vg_chunk.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/{region_id}/vg_chunk.log"
        threads: 8
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
            runtime=100*60,
            slurm_partition=choose_partition(100)
        shell:
            """
            vg chunk -a {input.sorted_gaf} -F -g -x {input.sampled_gbz} -p {params.region} -S {input.snarls} --trace -t {threads} -b {params.prefix} > {output.subgraph_vg}
            mv {params.prefix}{params.region_underscore}.gaf {output.chunked_gaf}
            vg convert -p {output.subgraph_vg} > {output.subgraph_pg_vg}
            vg convert -f {output.subgraph_pg_vg} > {output.subgraph_pg_gfa}
            """

    rule annotate_gam_gaf_with_refpos:
        output:
            chunked_gaf_uniq="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.uniq.gaf",
            chunked_gam="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.gam",
            annotated_gam="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.refpos.gam",
            annotated_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.refpos.gaf"
        input:
            script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/add_refpos_gam_to_gaf.py",
            sampled_gbz="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.gbz",
            chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.gaf"
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/annotate_gam_gaf_with_refpos.benchmark.txt"
        threads: 16
        shell:
            """
            sort {input.chunked_gaf} | uniq > {output.chunked_gaf_uniq}
            vg convert --gaf-to-gam {output.chunked_gaf_uniq} {input.sampled_gbz} > {output.chunked_gam}
            vg annotate -x {input.sampled_gbz} -a {output.chunked_gam} -p -P -t {threads} > {output.annotated_gam}
            vg view -a {output.annotated_gam} -j | python {input.script} --gaf {output.chunked_gaf_uniq} --output {output.annotated_gaf}
            """

    rule vg_index_dist:
        output:
            subgraph_pg_dist="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.pg.dist"
        input:
            subgraph_pg_vg="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.pg.vg"
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/vg_index_dist.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/{region_id}/vg_index_dist.log"
        threads: 16
        shell:
            """
            vg index {input.subgraph_pg_vg} --threads {threads} --dist-name {output.subgraph_pg_dist}
            """

    rule generate_anchors_dictionary:
        output:
            anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.pkl"
        input:
            vg_anchors_config=config["VG_ANCHORS"]["conf"],
            subgraph_pg_dist="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.pg.dist",
            sampled_pg_vg="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.pg.vg"
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/{region_id}/generate_anchors_dictionary.log"
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
            runtime=100*60,
            slurm_partition=choose_partition(100)
        shell:
            """
            vg-anchors --config {input.vg_anchors_config} build --graph {input.sampled_pg_vg} --index {input.subgraph_pg_dist} \
                --output-prefix results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph > {log} 2>&1
            """

    rule get_anchors_from_gaf:
        output:
            anchors="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.extended.jsonl",
            params_log="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/params_run.log",
            **({} if not config.get("RUN_DEBUGGING") else {"read_processed_tsv": "results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.reads_processed.tsv"})
        input:
            vg_anchors_config=config["VG_ANCHORS"]["conf"],
            anchors_dictionary="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.pkl",
            chunked_gaf="results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.gaf",
            sampled_pg_vg="results_hs/graph/{sample_id}/{sample_id}-{k}-sampled.pg.vg",
            chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/shasta/{sample_id}.subregion.fasta"
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/{region_id}/get_anchors_from_gaf.log"
        params:
            processes=64
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=1000),
            runtime=100*60,
            slurm_partition=choose_partition(100)
        shell:
            """
            vg-anchors --config {input.vg_anchors_config} get-anchors \
                --dictionary {input.anchors_dictionary} \
                --threads {params.processes} \
                --graph {input.sampled_pg_vg} \
                --alignment {input.chunked_gaf} \
                --fasta {input.chunked_fasta} \
                --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.anchors.json > {log} 2>&1
            """

rule chunk_fasta:
    output:
        chunked_fasta = "results_hs/hs-{k}/{sample_id}/{region_id}/shasta/{sample_id}.subregion.fasta"
    input:
        fasta="results/reads/{sample_id}.fasta",
        chunked_gaf=(
            "results_hs/hs-{k}/{sample_id}/{region_id}/query/subgraph.gaf"
            if config.get("RUN_GBZ_QUERY", False)
            else "results_hs/hs-{k}/{sample_id}/{region_id}/chunk/subgraph.gaf"
        )
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/chunk_fasta.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/chunk_fasta.log"
    container:
        "docker://pegi3s/seqkit:latest"
    shell:
        """
        # Extract READ IDs from chunked GAF and then extract fasta records for these reads
        cut -f1 {input.chunked_gaf} | sort -u -T {config[TMPDIR]}| seqkit grep -f - {input.fasta} > {output.chunked_fasta}
        """

rule run_shasta_assembly:
	output: "results_hs/hs-{k}/{sample_id}/{region_id}/shasta/ShastaRun/Assembly.fasta"
	input:
		shasta=config["SHASTA"]["bin"],
		shasta_conf=config["SHASTA"]["conf"],
		chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/shasta/{sample_id}.subregion.fasta",
		anchors="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.extended.jsonl"
	benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/run_shasta_assembly.benchmark.txt"
	log: "logs/{sample_id}/hs-{k}/{region_id}/run_shasta_assembly.log"
	resources:
		mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=10000),
		runtime=60*60,
		slurm_partition=choose_partition(60)
	shell:
		"""
		## make a copy of the anchors.json for Shasta
		cp {input.anchors} results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/shasta/anchors.json
		rm -rf results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/shasta/ShastaRun
		
		{input.shasta} --input {input.chunked_fasta} --assemblyDirectory results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/shasta/ShastaRun \
			--config {input.shasta_conf} --anchors results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/shasta/anchors.json --memoryMode filesystem --memoryBacking disk > {log} 2>&1
		"""

rule get_extended_anchor_stats:
    output:
        anchor_reads_info="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/extended_anchor_reads_info.tsv",
        anchor_stats_dir=directory("results_hs/hs-{k}/{sample_id}/{region_id}/extended_anchor_stats")
    input:
        scripts_dir="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts",
        anchors="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.extended.jsonl",
        subregion_shasta_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/shasta/ShastaRun/Assembly.fasta",
        chunked_fasta="results_hs/hs-{k}/{sample_id}/{region_id}/shasta/{sample_id}.subregion.fasta"
    params:
        region=config['region']['chromosome'] + ":" + config['region']['start'] + "-" + config['region']['end'],
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_extended_anchor_stats.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/get_extended_anchor_stats.log"
    resources:
        mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=400),
        runtime=10*60,
        slurm_partition=choose_partition(10)
    shell:
        """
        mkdir -p results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/extended_anchor_stats
        # Since this JSON is confusing to interpret, convert it to a TSV...
        python3 {input.scripts_dir}/preprocess_vganchor_outfiles_extention.py -j {input.anchors}
        # Next, generate anchor stats plots...
        Rscript {input.scripts_dir}/get_extended_anchor_stats.R -d results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id} -o ont_r10_{wildcards.region_id} -r {params.region} > {log} 2>&1
        """

# #------- Generate anchor sequence TSV --------# (Can plug into rule above if needed)
# master_table_tsv=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/extended_anchor_stats/ont_r10_{wildcards.region_id}_anchors_master_table.tsv
# python {input.scripts_dir}/process_anchor_seqs.py -j {input.anchors} -m $master_table_tsv -f {input.chunked_fasta} -o results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/extended_anchor_stats/ont_r10_{wildcards.region_id}        



if config.get("RUN_DEBUGGING"):

    rule get_reliable_snarl_stats:
        output:
            reliable_snarl_stats_dir=directory("results_hs/hs-{k}/{sample_id}/{region_id}/reliable_snarl_stats"),
            snarl_compatibility_fractions="results_hs/hs-{k}/{sample_id}/{region_id}/reliable_snarl_stats/snarl_compatibility_fractions.tsv"
        input:
            anchors="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.extended.jsonl",
            scripts_dir="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts",
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_reliable_snarl_stats.benchmark.txt"
        log: "logs/{sample_id}/hs-{k}/{region_id}/get_reliable_snarl_stats.log"
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
            runtime=10*60,
            slurm_partition=choose_partition(10)
        shell:
            """
            #------- Generate reliable snarl stats --------#
            output_dir=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/reliable_snarl_stats
            reliable_snarls_tsv=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.anchors.json.reliable_snarls.tsv
            snarl_compatibility_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.anchors.json.snarl_compatibility.jsonl
            snarl_coverage_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.anchors.json.snarl_coverage.jsonl
            snarl_allelic_coverage_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.anchors.json.snarl_allelic_coverage.jsonl
            snarl_allelic_coverage_extended_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.anchors.json.snarl_allelic_coverage_extended.jsonl
            snarl_variant_type_json=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.anchors.json.snarl_variant_type.jsonl
            snarl_positions_tsv=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.sizes.tsv

            Rscript {input.scripts_dir}/reliable_snarl_stats_test.R $reliable_snarls_tsv $snarl_compatibility_json $snarl_coverage_json $snarl_allelic_coverage_json $snarl_allelic_coverage_extended_json $snarl_variant_type_json $snarl_positions_tsv $output_dir >> {log} 2>&1
            """
    
    rule get_debugging_files:
        output:
            nodes_info_tsv="results_hs/hs-{k}/{sample_id}/{region_id}/debugging/{region_id}_nodes_info.tsv",
            read_traversals_zip="results_hs/hs-{k}/{sample_id}/{region_id}/debugging/{region_id}_read_traversals.zip",
            snarls_bandage_csv="results_hs/hs-{k}/{sample_id}/{region_id}/debugging/{region_id}_snarls.bandage.csv"
        input:
            script_dir="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts",
            read_processed_tsv="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/subgraph.anchors.json.reads_processed.tsv",
            snarl_compatibility="results_hs/hs-{k}/{sample_id}/{region_id}/reliable_snarl_stats/snarl_compatibility_fractions.tsv"
        benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/get_debugging_files.benchmark.txt"
        resources:
            mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
            runtime=15*60,
            slurm_partition=choose_partition(15)
        shell:
            """
            python {input.script_dir}/process_read_processed_for_bandage.py {input.read_processed_tsv} -o {output.nodes_info_tsv}
            # selected columns for bandage and convert to CSV
            cut -f1,2,5,7 {output.nodes_info_tsv} | tr '\\t' ',' > results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/debugging/{wildcards.region_id}_nodes_info.bandage.csv
            
            # read alignments for bandage
            python {input.script_dir}/extract_read_alignments_for_bandage.py {input.read_processed_tsv} {output.nodes_info_tsv} -o {output.read_traversals_zip}

            # reliable/unreliable snarls for bandage
            snarl_dict="results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/anchors/subgraph.forward_dict.csv"
            python {input.script_dir}/map_snarls_bandage.py $snarl_dict {input.snarl_compatibility} {output.nodes_info_tsv} -o {output.snarls_bandage_csv}
            """

rule align_assembly_to_chunked_reference:
	output:
		paf="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.paf",
		csv="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.csv"
	input:
		analyzePaf_bin=config["ANALYSEPAF"]["bin"],
		hg002_reference_chunked="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/hg002.chunked.fasta",
		subregion_shasta_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/shasta/ShastaRun/Assembly.fasta"
	params:
		asm_preset=config["MINIMAP"]["asmPreset"]
	benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/align_assembly_to_reference_{asm_preset}.benchmark.txt"
	log: "logs/{sample_id}/hs-{k}/{region_id}/align_assembly_to_reference_{asm_preset}.log"
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
			--output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/assembly_alignment/{wildcards.sample_id}_{wildcards.region_id}_ZOOMED_shasta_to_hg002_minimap_{wildcards.asm_preset} \
			--pixelsPerMb 600 --minAlignmentQuality 1 --minAlignmentLength 0 > {log} 2>&1
		"""

rule run_displayPafAlignments:
	output:
		csv="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}_displayPaf.csv"
	input:
		displayPaf_bin=config["DISPLAYPAF"]["bin"],
		hg002_reference_chunked="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/hg002.chunked.fasta",
		subregion_shasta_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/shasta/ShastaRun/Assembly.fasta",
		paf="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.paf"
	params:
		asm_preset=config["MINIMAP"]["asmPreset"],
	benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/run_displayPafAlignments_{asm_preset}.benchmark.txt"
	log: "logs/{sample_id}/hs-{k}/{region_id}/run_displayPafAlignments_{asm_preset}.log"
	resources:
		mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
		runtime=10*60,
		slurm_partition=choose_partition(10)
	shell:
		"""
		{input.displayPaf_bin} \
			--paf {input.paf} -r {input.hg002_reference_chunked} -a {input.subregion_shasta_assembly} \
			--output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/assembly_alignment/{wildcards.sample_id}_{wildcards.region_id}_ZOOMED_shasta_to_hg002_minimap_{wildcards.asm_preset}_displayPaf \
			--minAlignmentQuality 1 --minAlignmentLength 0 > {log} 2>&1		
		"""

#-------Using shasta alignment coordinates to hg002, chunk hg002 reference for hifiasm-------#
rule chunk_hg002_reference_for_hifiasm_using_shasta_alignment_coordinates:
    output:
        hg002_reference_chunked_for_hifiasm="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.fasta",
        coords_tsv="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.coords.tsv"
    input:
        hg002_reference=config["HG002v101_ref"],
        script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/chunk_hg002_reference_for_hifiasm_using_shasta_alignments.py",
        shasta_csv="results_hs/hs-{k}/{sample_id}/{region_id}/assembly_alignment/{sample_id}_{region_id}_ZOOMED_shasta_to_hg002_minimap_{asm_preset}.csv",
    params:
        asm_preset=config["MINIMAP"]["asmPreset"]
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/chunk_hg002_reference_for_hifiasm_using_shasta_alignment_coordinates_{asm_preset}.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/chunk_hg002_reference_for_hifiasm_using_shasta_alignment_coordinates_{asm_preset}.log"
    resources:
        mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
        runtime=10*60,
        slurm_partition=choose_partition(10)
    shell:
        """
        python3 {input.script} {input.shasta_csv} {input.hg002_reference} \
            results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/hifiasm_assembly \
            hg002.chunked.{wildcards.asm_preset} > {log} 2>&1
        """


rule extract_hifiasm_subregion_assembly:
    output:
        hifiasm_subregion_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.subregion.fasta"
    input:
        script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/extract_subregion_contigs.py",
        bam="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_hifiasm_to_hg002_minimap_{asm_preset}.bam",
        bai="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_hifiasm_to_hg002_minimap_{asm_preset}.bam.bai",
        coords_tsv="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.coords.tsv",
        hifiasm_fasta="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}.hifiasm.fasta"
    params:
        asm_preset=config["MINIMAP"]["asmPreset"]
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/extract_hifiasm_subregion_assembly_{asm_preset}.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/extract_hifiasm_subregion_assembly_{asm_preset}.log"
    resources:
        mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
        runtime=10*60,
        slurm_partition=choose_partition(10)
    shell:
        """
        python3 {input.script} {input.bam} {input.coords_tsv} {input.hifiasm_fasta} {output.hifiasm_subregion_assembly} > {log} 2>&1
        """


rule align_hifiasm_subregion_assembly_to_chunked_hg002_reference:
    output:
        paf="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_hifiasm_subregion_to_hg002_minimap_{asm_preset}.paf",
        csv="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_hifiasm_subregion_to_hg002_minimap_{asm_preset}.csv"
    input:
        analyzePaf_bin=config["ANALYSEPAF"]["bin"],
        hg002_reference_chunked_for_hifiasm="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.fasta",
        hifiasm_subregion_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.subregion.fasta"
    params:
        asm_preset=config["MINIMAP"]["asmPreset"]
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/align_hifiasm_subregion_assembly_to_chunked_hg002_reference_{asm_preset}.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/align_hifiasm_subregion_assembly_to_chunked_hg002_reference_{asm_preset}.log"
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
        runtime=60*60,
        slurm_partition=choose_partition(60)
    shell:
        """
        minimap2 -t {threads} -I 20G -cx {params.asm_preset} -K 1M --eqx --cs {input.hg002_reference_chunked_for_hifiasm} {input.hifiasm_subregion_assembly} > {output.paf}
        {input.analyzePaf_bin} \
            --input {output.paf} \
            --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/hifiasm_assembly/{wildcards.sample_id}_{wildcards.region_id}_hifiasm_subregion_to_hg002_minimap_{wildcards.asm_preset} \
            --pixelsPerMb 1000 --minAlignmentQuality 0 > {log} 2>&1
        """


#------ Align r_utg hifiasm assembly to chunked hg002 reference for hifiasm assembly -------#
rule extract_r_utg_hifiasm_subregion_assembly:
    output:
        r_utg_hifiasm_subregion_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.r_utg.subregion.fasta"
    input:
        script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/extract_subregion_contigs.py",
        bam="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_r_utg_to_hg002_minimap_{asm_preset}.bam",
        bai="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_r_utg_to_hg002_minimap_{asm_preset}.bam.bai",
        coords_tsv="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.coords.tsv",
        hifiasm_r_utg_fasta="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.r_utg.fasta"
    params:
        asm_preset=config["MINIMAP"]["asmPreset"]
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/extract_r_utg_hifiasm_subregion_assembly_{asm_preset}.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/extract_r_utg_hifiasm_subregion_assembly_{asm_preset}.log"
    resources:
        mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
        runtime=10*60,
        slurm_partition=choose_partition(10)
    shell:
        """
        python3 {input.script} {input.bam} {input.coords_tsv} {input.hifiasm_r_utg_fasta} {output.r_utg_hifiasm_subregion_assembly} > {log} 2>&1
        """

rule align_r_utg_hifiasm_subregion_assembly_to_chunked_hg002_reference:
    output:
        paf="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}.paf",
        csv="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}.csv"
    input:
        analyzePaf_bin=config["ANALYSEPAF"]["bin"],
        hg002_reference_chunked_for_hifiasm="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.fasta",
        r_utg_hifiasm_subregion_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.r_utg.subregion.fasta"
    params:
        asm_preset=config["MINIMAP"]["asmPreset"]
    benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/align_r_utg_hifiasm_subregion_assembly_to_chunked_hg002_reference_{asm_preset}.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/{region_id}/align_r_utg_hifiasm_subregion_assembly_to_chunked_hg002_reference_{asm_preset}.log"
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
        runtime=60*60,
        slurm_partition=choose_partition(60)
    shell:
        """
        minimap2 -t {threads} -I 20G -cx {params.asm_preset} -K 1M --eqx --cs {input.hg002_reference_chunked_for_hifiasm} {input.r_utg_hifiasm_subregion_assembly} > {output.paf}
        {input.analyzePaf_bin} \
            --input {output.paf} \
            --output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/hifiasm_assembly/{wildcards.sample_id}_{wildcards.region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{wildcards.asm_preset} \
            --pixelsPerMb 100 --minAlignmentQuality 1 > {log} 2>&1
        """

rule run_displayPafAlignments_for_r_utg_hifiasm_subregion_assembly:
	output:
		filtered_paf="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}.filtered.paf",
		csv="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_alignment/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}_displayPaf.csv",
	input:
		displayPaf_bin=config["DISPLAYPAF"]["bin"],
		hg002_reference_chunked_for_hifiasm="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/hg002.chunked.{asm_preset}.for_hifiasm.fasta",
		r_utg_hifiasm_subregion_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.r_utg.subregion.fasta",
		paf="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}.paf",
		csv="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}_{region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{asm_preset}.csv",
		py_script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/filter_paf_to_csv_alignments.py"
	params:
		asm_preset=config["MINIMAP"]["asmPreset"],
	benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/run_displayPafAlignments_for_r_utg_hifiasm_subregion_assembly_{asm_preset}.benchmark.txt"
	log: "logs/{sample_id}/hs-{k}/{region_id}/run_displayPafAlignments_for_r_utg_hifiasm_subregion_assembly_{asm_preset}.log"
	resources:
		mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
		runtime=10*60,
		slurm_partition=choose_partition(10)
	shell:
		"""
		# FILTER paf to only include alignments in the CSV
		python {input.py_script} --csv {input.csv} --paf {input.paf} --output {output.filtered_paf}
		
		{input.displayPaf_bin} \
			--paf {output.filtered_paf} -r {input.hg002_reference_chunked_for_hifiasm} -a {input.r_utg_hifiasm_subregion_assembly} \
			--output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/hifiasm_alignment/{wildcards.sample_id}_{wildcards.region_id}_r_utg_hifiasm_subregion_to_hg002_minimap_{wildcards.asm_preset}_displayPaf \
			--minAlignmentQuality 1 --minAlignmentLength 0 > {log} 2>&1		
		"""

#------ Align shasta assembly to hifiasm assembly --------#
rule shasta_to_hifiasm_alignment:
	output:
		shasta_to_hifiasm_alignment_paf="results_hs/hs-{k}/{sample_id}/{region_id}/shasta_to_hifiasm_alignment/{sample_id}_{region_id}_shasta_to_hifiasm_minimap_{asm_preset}.paf",
		shasta_to_hifiasm_alignment_bam="results_hs/hs-{k}/{sample_id}/{region_id}/shasta_to_hifiasm_alignment/{sample_id}_{region_id}_shasta_to_hifiasm_minimap_{asm_preset}.bam",
		shasta_to_hifiasm_alignment_bai="results_hs/hs-{k}/{sample_id}/{region_id}/shasta_to_hifiasm_alignment/{sample_id}_{region_id}_shasta_to_hifiasm_minimap_{asm_preset}.bam.bai",
	input:
		analyzePaf_bin=config["ANALYSEPAF"]["bin"],
		# the subregion assembly from the final hifiasm assembly
		# hifiasm_subregion_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.subregion.fasta",
		# the subregion assembly from the r_utg hifiasm assembly
		hifiasm_subregion_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/hifiasm_assembly/{sample_id}.{asm_preset}.hifiasm.r_utg.subregion.fasta",
		shasta_assembly="results_hs/hs-{k}/{sample_id}/{region_id}/shasta/ShastaRun/Assembly.fasta",
	params:
		asm_preset=config["MINIMAP"]["asmPreset"]
	benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/shasta_to_hifiasm_alignment_{asm_preset}.benchmark.txt"
	log: "logs/{sample_id}/hs-{k}/{region_id}/shasta_to_hifiasm_alignment_{asm_preset}.log"
	threads: 16
	resources:
		mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
		runtime=60*60,
		slurm_partition=choose_partition(60)
	shell:
		"""
		# Generate index for hifiasm assembly
		samtools faidx {input.hifiasm_subregion_assembly}

		minimap2 -t {threads} -I 20G -cx {params.asm_preset} -K 1M --eqx --cs {input.hifiasm_subregion_assembly} {input.shasta_assembly} > {output.shasta_to_hifiasm_alignment_paf}

		minimap2 -t {threads} -I 20G -ax {params.asm_preset} -K 1M --eqx --cs {input.hifiasm_subregion_assembly} {input.shasta_assembly} | \
			samtools view -bS -F 256 - | \
			samtools sort -@ {threads} -o {output.shasta_to_hifiasm_alignment_bam} -
		
		# Index the BAM file
		samtools index {output.shasta_to_hifiasm_alignment_bam} {output.shasta_to_hifiasm_alignment_bai}

		# Generate analysePaf output
		{input.analyzePaf_bin} \
			--input {output.shasta_to_hifiasm_alignment_paf} \
			--output results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/hifiasm_to_shasta_alignment/{wildcards.sample_id}_{wildcards.region_id}_hifiasm_to_shasta_minimap_{wildcards.asm_preset} \
			--pixelsPerMb 5000 --minAlignmentQuality 0 > {log} 2>&1
		"""

rule generate_alignment_plot_for_shasta_to_hifiasm_alignment:
	output:
		alignment_plots_pdf="results_hs/hs-{k}/{sample_id}/{region_id}/shasta_to_hifiasm_alignment/{sample_id}_{region_id}_shasta_to_hifiasm_{asm_preset}_alignment_plots.pdf"
	input:
		r_script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/generate_alignment_diagonal_plot.R",
		shasta_to_hifiasm_alignment_paf="results_hs/hs-{k}/{sample_id}/{region_id}/shasta_to_hifiasm_alignment/{sample_id}_{region_id}_shasta_to_hifiasm_minimap_{asm_preset}.paf"
	params:
		asm_preset=config["MINIMAP"]["asmPreset"]
	benchmark: "benchmarks/{sample_id}/hs-{k}/{region_id}/generate_alignment_plot_for_shasta_to_hifiasm_alignment_{asm_preset}.benchmark.txt"
	log: "logs/{sample_id}/hs-{k}/{region_id}/generate_alignment_plot_for_shasta_to_hifiasm_alignment_{asm_preset}.log"
	resources:
		mem_mb=lambda wildcards, input, attempt: get_mem_mb(input, buffer_factor=1.2, min_mb=100),
		runtime=60*60,
		slurm_partition=choose_partition(60)
	shell:
		"""
		# Create output directory
		mkdir -p results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/shasta_to_hifiasm_alignment
		
		# Generate alignment diagonal plots
		Rscript {input.r_script} \
			{input.shasta_to_hifiasm_alignment_paf} \
			results_hs/hs-{wildcards.k}/{wildcards.sample_id}/{wildcards.region_id}/shasta_to_hifiasm_alignment \
			{wildcards.sample_id}_{wildcards.region_id}_shasta_to_hifiasm_{wildcards.asm_preset} > {log} 2>&1
		"""

rule generate_run_summary:
    output:
        pga_log="results_hs/hs-{k}/{sample_id}/{region_id}/pga_run.log"
    input:
        params_log="results_hs/hs-{k}/{sample_id}/{region_id}/anchors/params_run.log",
        shasta_conf=config["SHASTA"]["conf"],
        script="/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/workflow/scripts/generate_runlog.py"
    params:
        run_mode=config['RUN_MODE'],
        region_id=config['region_id'],
        asm_preset=config['MINIMAP']['asmPreset'],
        read_type=config.get('READ_TYPE'),
        run_gbz_query=config.get('RUN_GBZ_QUERY'),
        use_full_graph=config.get('USE_FULL_GRAPH'),
        run_debugging=config.get('RUN_DEBUGGING')
    shell:
        """
        python3 {input.script} --params-log {input.params_log} --shasta-conf {input.shasta_conf} --output-log {output.pga_log} --run-mode {params.run_mode} --region-id {params.region_id} --asm-preset {params.asm_preset} \
            --read-type {params.read_type} --run-gbz-query {params.run_gbz_query} --use-full-graph {params.use_full_graph} --run-debugging {params.run_debugging}
        """
