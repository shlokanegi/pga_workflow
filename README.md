# Pangenome Guided Genome Assembly Workflow


## Directory Structure
* `workflow/`: Contains the main Snakefile (Other snakefiles are inside `other_snakefiles/`)
* `config/`: Contains snakemake configuration files for different datasets
* `benchmark/`: Runtime and resource usage statistics for individual Snakemake rules
* `logs/`: Standard output and error logs for each rule execution, useful for debugging
* `results*/`: Output files generated by the pipeline.


## Run the workflow
From the pga_workflow/ directory, run:
```sh
snakemake --use-singularity --singularity-args "-B /private/groups/" --configfile config/config_ontR10.yaml --cores 128 --printshellcmds
```
To resume the workflow from where it left off (e.g., after interruption), use:
```sh
snakemake --use-singularity --singularity-args "-B /private/groups/" --configfile config/config_ontR10.yaml --cores 128 --printshellcmds --rerun-incomplete 
## Sometimes, it asks to unlock and then rerun
snakemake --use-singularity --singularity-args "-B /private/groups/" --configfile config/config_ontR10.yaml --cores 128 --printshellcmds --unlock
```
This will rerun any jobs that were incomplete or failed during the previous run.

## Workflow Control with RUN_MODE
Different execution modes are controlled by the `RUN_MODE` parameter in the `config/config_ontR10.yaml` file.

### RUN_MODE Options:
1. **`"no_positive_control"`**: Runs only the sample workflow
2. **`"positive_control_only"`**: Runs only the positive control workflow. Uses the HG002-included graph for analysis
3. **`"all"`** (default): Runs both the sample workflow and the positive control workflow

Edit the `RUN_MODE` parameter in your config file (e.g., `config/config_ontR10.yaml`) and run the workflow as usual.

## Analyzing Anchor Coverage

The workflow generates anchor coverage statistics during the assembly process. These statistics are stored in a JSON file with the extension `.coverage.json` in the anchors output directory. To analyze and visualize these statistics:

**Generate Coverage Statistics**:
```sh
python vg_assembly/analyze_anchor_coverage.py \
    --coverage-file results_hs/hs-16/LGA80510/rccx_test/anchors/subgraph.anchors.json.jsonl.coverage.json \
    --output-plot results_hs/hs-16/LGA80510/rccx_test/anchors/coverage_histogram.png
```
