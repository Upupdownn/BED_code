#!/usr/bin/env python3
'''
End-Motif Decomposition and SVM Workflow

Inputs:
    - input_dir: Directory containing input BAM or TSV files (rows=fragments, columns=chr/start/end/mapq/strand)
    - output_dir: Output directory for all results
    - --tb_file: Reference genome 2bit file for sequence extraction
    - --train_info_file: TSV file with sample info (at least 'id' and 'label' columns for training)
    - --val_info_file: Optional TSV file with sample info for validation
    - --id_col: Column name for sample ID [default: id]
    - --label_col: Column name for binary label (0/1) [default: label]
    - -p/--processes: Number of parallel processes [default: 30]

Description:
    This script orchestrates a workflow for processing BAM or TSV fragment files to extract end-motif features,
    perform Bayesian End-motif Decomposition (BED), and train/validate SVM models for classification.
    It handles conversion from BAM to TSV if needed, feature extraction, merging, reference generation,
    pilot testing for learning rates, full decomposition, and SVM training/validation.

Workflow Steps:
    1. If input is BAM, convert to fragment TSV files.
    2. Extract end-motif frequencies for each sample.
    3. Merge end-motif features and compute MDS features.
    4. Generate reference background and pilot dataset from training info.
    5. Run BED in pilot mode to test learning rates.
    6. Run BED decomposition on all samples.
    7. Train and validate SVM models on EDM and ACW features.
    8. Plot

Outputs:
    - fragment_tsvs/: TSV files from BAM conversion (if applicable)
    - sample_edm_frequencies/: Per-sample end-motif frequency TSVs
    - features/: Merged EDM.tsv, MDS.tsv, and BED outputs (AW.tsv, BW.tsv, etc.)
    - pilot/: Pilot EDM.tsv, BED outputs, and pilot_lr_search.tsv
    - svm_edm/: SVM models and scores for EDM features
    - svm_acw/: SVM models and scores for ACW features

Dependencies:
    Requires external scripts in PATH: bam_to_tsv.py, extract_edm_features.py, merge_edm_features.py,
    prepare_ref_pilot.py, bayesian_decompo.py, svm_train_val.py

Example Usage:
    ./BED_workflow.py /path/to/input /path/to/output --tb_file ref.2bit --train_info_file train.tsv
'''

import argparse
import logging
from pathlib import Path
import shutil
import subprocess
import os


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

script_bam2tsv = '01_bam_to_tsv.py'
script_extract_edm = '02_extract_edm_features.py'
script_merge_edm = '03_merge_edm_features.py'
script_ref_pilot = '04_generate_ref_pilot.py'
script_bed = '05_bayesian_decompo.py'
script_svm_train_val = '06_svm_train_val.py'
script_list = [script_bam2tsv, script_extract_edm, script_merge_edm, script_ref_pilot, script_bed, script_svm_train_val]
scripts = None      # update by check_scripts_available

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("input_dir", help="Input directory containing BAM or TSV files")
    parser.add_argument("output_dir", help="Output directory for results")
    parser.add_argument("--tb_file", required=True, help="Reference genome 2bit file")
    parser.add_argument("--train_info_file", required=True, help="TSV file with at least 'id' and 'label' columns for training")
    parser.add_argument("--val_info_file", default=None, help="TSV file with at least 'id' and 'label' columns for validation")
    parser.add_argument("--id_col", default="id", help="Column name for sample ID [default: id]")
    parser.add_argument("--label_col", default="label", help="Column name for label (0/1) [default: label]")
    parser.add_argument("-p", "--processes", type=int, default=10, help="Number of processes [default: 10]")
    return parser.parse_args()

def check_scripts_available(script_list: list):
    """
    Check if required scripts are available.
    Priority: 1. Current script's directory 2. System PATH
    """
    bin_dir: Path = Path(__file__).parent.resolve() / 'bin'

    found_scripts = {}
    for script in script_list:
        local_path: Path = bin_dir / script
        
        if local_path.exists() and os.access(local_path, os.X_OK):
            found_scripts[script] = str(local_path)
        else:
            system_path = shutil.which(script)
            if system_path:
                found_scripts[script] = system_path
            else:
                raise FileNotFoundError(f"Required script '{script}' not found in '{bin_dir}' or system PATH.")

    return found_scripts
        
def detect_file_type(dir_path):
    """Detect if directory contains BAM or TSV files."""
    for path in Path(dir_path).iterdir():
        if path.is_file():
            filename = path.name.lower()
            if filename.endswith(".bam"):
                return "bam"
            if filename.endswith(".tsv") or ".tsv." in filename:
                return "tsv"
    raise RuntimeError("No tsv or bam file found in directory")

def bam2tsv(bam_dir, tsv_dir, processes=20):
    """Convert BAM files to fragment TSV files."""
    os.makedirs(tsv_dir, exist_ok=True)
    for bam_file in Path(bam_dir).rglob("*.bam"):
        sample_id = bam_file.stem
        tsv_file = os.path.join(tsv_dir, f"{sample_id}.tsv")

        cmd_args = [scripts[script_bam2tsv], bam_file, tsv_file, '-p', str(processes)]
        logger.info(f"Converting {bam_file}")
        subprocess.run(cmd_args, check=True)

def extract_EDM(frag_dir, out_dir, tb_file, processes=20):
    """Extract end-motif features from fragment TSV files."""
    os.makedirs(out_dir, exist_ok=True)
    for frag_file in Path(frag_dir).rglob("*.tsv*"):
        sample_id = frag_file.name.split('.tsv')[0]
        out_file = os.path.join(out_dir, f"{sample_id}.tsv")
        cmd_args = [scripts[script_extract_edm], frag_file, tb_file, out_file, '-p', str(processes)]
        subprocess.run(cmd_args, check=True)

def merge_EDM(sample_edm_dir, edm_file, mds_file):
    """Merge sample end-motif features and compute MDS features."""
    cmd_args = [scripts[script_merge_edm], sample_edm_dir, '--merged_output', edm_file, '--mds_output', mds_file]
    subprocess.run(cmd_args, check=True)

def generate_ref_and_pilot(edm_file, info_file, id_col, label_col, ref_out, pilot_out):
    """Generate reference background and pilot dataset."""
    cmd_args = [scripts[script_ref_pilot], edm_file, info_file, '--id_col', id_col, '--label_col', label_col, 
                '--ref_out', ref_out, '--pilot_out', pilot_out]
    subprocess.run(cmd_args, check=True)


def BED_func(edm_file, ref_file, out_dir, pilot=False, processes=20):
    """Run Bayesian End-motif Decomposition (BED)."""
    cmd_args = [scripts[script_bed], edm_file, ref_file, out_dir, '-p', str(processes)]
    if pilot: cmd_args.append('--pilot')
    subprocess.run(cmd_args, check=True)
    
def SVM_train_val(feature_file, train_info_file, val_info_file, out_dir):
    """Train and validate SVM models."""
    train_score_file = os.path.join(out_dir, f"score_train.tsv")
    val_score_file = os.path.join(out_dir, f"score_val.tsv")
    model_dir = os.path.join(out_dir, "models")

    # Train mode
    cmd_args = [scripts[script_svm_train_val], feature_file, train_info_file, train_score_file, 
                '--model_dir', model_dir, '--mode', 'train']
    subprocess.run(cmd_args, check=True)

    # Validate mode (if validation file provided)
    if val_info_file is not None:
        cmd_args = [scripts[script_svm_train_val], feature_file, val_info_file, val_score_file, 
                    '--model_dir', model_dir, '--mode', 'validate']
        subprocess.run(cmd_args, check=True)


def main():
    args = parse_args()
    input_dir = args.input_dir
    out_dir = args.output_dir
    processes = args.processes
    tb_file = args.tb_file
    train_info_file = args.train_info_file
    val_info_file = args.val_info_file
    id_col, label_col = args.id_col, args.label_col

    # Check if required scripts are available
    global scripts
    scripts = check_scripts_available(script_list)

    # Convert BAM to TSV if necessary
    if detect_file_type(input_dir) == 'bam':
        logger.info(f"Converting BAM to fragment TSV file.")
        frag_dir = os.path.join(out_dir, 'fragment_tsvs')
        bam2tsv(input_dir, frag_dir, processes)
        logger.info(f"BAM to TSV conversion completed.")
    else:
        logger.info(f"Starting workflow with existing fragment TSV files.")
        frag_dir = input_dir

    # Extract end-motif features
    sample_edm_dir = os.path.join(out_dir, "sample_edm_frequencies")
    extract_EDM(frag_dir, sample_edm_dir, tb_file, processes)

    # Merge features and compute MDS
    feature_dir = os.path.join(out_dir, 'features')
    edm_file = os.path.join(feature_dir, 'EDM.tsv')
    mds_file = os.path.join(feature_dir, 'MDS.tsv')
    merge_EDM(sample_edm_dir, edm_file, mds_file)

    # Prepare reference and pilot
    ref_file = os.path.join(out_dir, 'reference_background.npy')
    pilot_out_dir = os.path.join(out_dir, 'pilot')
    pilot_edm_file = os.path.join(pilot_out_dir, 'EDM.tsv')
    generate_ref_and_pilot(edm_file, train_info_file, id_col, label_col, ref_out=ref_file, pilot_out=pilot_edm_file)

    # Run BED in pilot mode
    BED_func(pilot_edm_file, ref_file, pilot_out_dir, pilot=True, processes=processes)

    # Run full BED decomposition
    BED_func(edm_file, ref_file, feature_dir, pilot=False, processes=processes)

    # SVM on EDM features
    svm_edm_dir = os.path.join(out_dir, 'svm_edm')
    SVM_train_val(edm_file, train_info_file, val_info_file, svm_edm_dir)

    # SVM on ACW features
    svm_acw_dir = os.path.join(out_dir, 'svm_acw')
    acw_file = os.path.join(feature_dir, 'ACW.tsv')
    SVM_train_val(acw_file, train_info_file, val_info_file, svm_acw_dir)


if __name__ == "__main__":
    main()

