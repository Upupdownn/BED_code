#!/usr/bin/env python3
"""
Bayesian End-motif Decomposition (BED) for multiple samples.

Inputs:
    - frequency_file:   TSV file with motif frequencies
                        Rows: samples, Columns: 256 4-mer motifs, Values: frequency
    - reference_file:   Reference background motif frequencies (.npy file)
                        Typically from healthy controls, shape: (256,)

Outputs in normal mode (default):
    The following TSV files will be generated in the specified output_dir:
    - AW.tsv          : Aberrant weight (non-background weight) per sample
    - BW.tsv          : Background weight per sample
    - BCW.tsv         : Background conditional weight (alpha)
    - ACW.tsv         : Aberrant conditional weight (1-alpha)
    - BCP.tsv         : Background conditional probability
    - ACP.tsv         : Aberrant conditional probability
    - StopInfo.tsv    : Convergence information (early stopping, final epoch, final loss)

Pilot mode (--pilot):
    Purpose: Tests multiple learning rates on the input dataset to help select an optimal initial learning rate.
    Behavior: Runs the full decomposition process for the same data using different learning rates,
              then computes average convergence epochs and early-stopping rate per lr for comparison.

    Outputs in pilot mode (saved to output_dir):
    - StopInfo_lr_<lr>.tsv : Full StopInfo table for each tested learning rate (e.g., StopInfo_lr_1.0.tsv)
    - pilot_lr_search.tsv  : Summary table with columns:
        - lr                  : Tested learning rate
        - epoch_mean          : Average number of epochs across all samples (lower is better)
        - early_stop_rate_%   : Percentage of samples that reached early stopping (higher is better)

Example usage:
    - Normal mode
    bayesian_decompo.py freq.tsv ref.npy results/ --lr 1.2 --processes 20

    - Pilot mode - learning rate exploration
    bayesian_decompo.py pilot_freq.tsv ref.npy pilot_results/ --pilot --processes 10
    
Supports multi-processing with 'spawn' start method for CUDA compatibility.
"""

import logging
import argparse
import numpy as np
import pandas as pd
import os
import torch
from torch import Tensor
from torch.optim.lr_scheduler import ReduceLROnPlateau
import torch.multiprocessing as tmp


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("frequency_file", help="Input TSV file: rows=samples, columns=256 motifs, value=motif frequency")
    parser.add_argument("reference_file", help="Reference background frequencies (.npy)")
    parser.add_argument("output_dir", help="Output directory for result TSV files")
    parser.add_argument("--alpha_init", type=float, default=4.596, help="Initial alpha value [default: 4.596]")
    parser.add_argument("--max_epochs", type=int, default=10000, help="Maximum epochs per sample [default: 10000]")
    parser.add_argument("--patience", type=int, default=50, help="Early stopping patience [default: 50]")
    parser.add_argument("--lr", type=float, default=1.0, help="Initial learning rate [default: 1.0]")
    parser.add_argument('-p', "--processes", type=int, default=20, help="Number of processes [default: 20]")
    parser.add_argument("--pilot", action="store_true", default=False, 
                        help="Pilot mode: use frequency_file as exploration data to test different lr values")

    return parser.parse_args()

class Utils():
    @staticmethod
    def kmer_list(k=4) -> list:
        """Generate all possible k-mers recursively."""
        bases = "ACGT"
        if k == 1:
            return list(bases)
        else:
            kmers = []
            for k_minus_mer in Utils.kmer_list(k - 1):
                for base in bases:
                    kmers.append(k_minus_mer + base)
            return kmers
        
    @staticmethod
    def BED_res_tidy(BED_res: dict, index):
        """Convert BED results to DataFrames."""
        df_dict = dict()
        for key, value in BED_res.items():
            columns = [key] if key in ['AW', 'BW'] else Utils.kmer_list(4)
            df = pd.DataFrame(value, index=index, columns=columns)
            df.index.name = 'id'
            df_dict[key] = df
        
        return df_dict

def compute_loss(freqs: Tensor, refs: Tensor, alpha_raw: Tensor):
    """Compute MSE loss between predicted and reference frequencies."""
    alpha = torch.sigmoid(alpha_raw)                        # Background conditional weight (BCW)
    BG_weight = (freqs @ alpha.T).diagonal(offset=0)        # Background weight (BW)
    scale_factor = 1.0 / (refs.mean().detach() + 1e-8)
    pred = freqs * alpha / BG_weight.reshape(-1, 1)
    loss = torch.mean(((pred * scale_factor) - (refs * scale_factor)) ** 2)
    return loss

def bayesian_decomposition(
        freqs: np.ndarray,        # Observed motif frequencies (samples x 256)
        refs: np.ndarray,         # Reference background frequencies (1 x 256)
        alpha_init=4.596,
        max_epochs=10000,
        abs_epsilon=1e-7,
        rel_epsilon=1e-10,
        N=50,
        lr=1,
        device='cuda', 
        ) -> tuple[dict, dict]:
    """Perform Bayesian End-motif Decomposition for one or more samples."""
    refs = np.tile(refs, (freqs.shape[0], 1))
    freqs = torch.tensor(freqs, dtype=torch.float64).to(device)
    refs = torch.tensor(refs, dtype=torch.float64).to(device)
    alpha_raw = torch.full(freqs.shape, alpha_init, dtype=torch.float64, requires_grad=True, device=device)

    optimizer = torch.optim.Adam([alpha_raw], lr=lr)
    scheduler = ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5,
        patience=max(10, int(N * 0.6)),
        threshold=1e-8, threshold_mode='abs',
        cooldown=10, min_lr=1e-8
    )

    pre_loss = float('inf')
    counter = 0
    ep = 0
    is_early_stop, epoch_at_stop, loss_at_stop = False, 0, 0
    while ep < max_epochs:
        loss = compute_loss(freqs, refs, alpha_raw)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        curr_loss = loss.item()
        delta = curr_loss - pre_loss
        scheduler.step(curr_loss)

        # Check absolute convergence
        if curr_loss < abs_epsilon:
            is_early_stop, epoch_at_stop, loss_at_stop = True, ep, curr_loss
            break
        
        # Check relative change
        if abs(delta) < rel_epsilon:
            counter += 1
            if counter >= N:
                is_early_stop, epoch_at_stop, loss_at_stop = True, ep, curr_loss
                break
        else:
            counter = 0
        
        pre_loss = curr_loss
        ep += 1
    
    if ep >= max_epochs:
        is_early_stop, epoch_at_stop, loss_at_stop = False, max_epochs, curr_loss
    stop_info = {'early': is_early_stop, 'epoch': epoch_at_stop, 'loss': loss_at_stop}

    alpha = torch.sigmoid(alpha_raw)
    bg_cond_weight = alpha                                  # BCW
    ab_cond_weight = 1 - alpha                              # ACW
    Bg_weight = (freqs @ alpha.T).diagonal(offset=0)        # BW
    Ab_weight = 1 - Bg_weight                               # AW
    bg_cond_prob = freqs * alpha / Bg_weight.reshape(-1, 1)                 # BCP
    ab_cond_prob = freqs * (1 - alpha) / (1 - Bg_weight.reshape(-1, 1))     # ACP

    feature_dict = {
        'BW':Bg_weight.detach().cpu().numpy(),
        'AW':Ab_weight.detach().cpu().numpy(),
        'BCW':bg_cond_weight.detach().cpu().numpy(),
        'ACW':ab_cond_weight.detach().cpu().numpy(),
        'BCP':bg_cond_prob.detach().cpu().numpy(),
        'ACP':ab_cond_prob.detach().cpu().numpy()
    }

    return stop_info, feature_dict

def process_single_sample(args):
    """Entry point for each process."""
    sample_id, freq, ref_bg, params = args
    freq = np.reshape(freq, (1, -1))  # Ensure 2D input for bayesian_decomposition
    stop_info, feature_dict = bayesian_decomposition(
        freq, ref_bg,
        alpha_init=params['alpha_init'],
        max_epochs=params['max_epochs'],
        N=params['N'],
        lr=params['lr'],
        device='cuda' if torch.cuda.is_available() else 'cpu'
    )
    return sample_id, stop_info, Utils.BED_res_tidy(feature_dict, index=[sample_id])


def run_multi_process_BED(freq_df: pd.DataFrame, ref_bg: np.ndarray,
                n_processes=10, alpha_init=4.596, max_epochs=10000, N=50, lr=1.0):
    """Multi-process BED decomposition for all samples in freq_df."""

    n_samples = len(freq_df)
    params = {
        'alpha_init': alpha_init,
        'max_epochs': max_epochs,
        'N': N,
        'lr': lr
    }

    tasks = [
        (freq_df.index[i], freq_df.iloc[i].to_numpy(), ref_bg, params)
        for i in range(n_samples)
    ]

    logger.info(f"Starting multi-processing BED decomposition ({n_processes} processes)...")
    tmp.set_start_method('spawn', force=True)
    
    with tmp.Pool(processes=n_processes) as pool:
        results = pool.map(process_single_sample, tasks)

    # Collect results
    stop_info_list = []
    df_dict = {}
    for sample_id, stop_info, feature_dict in results:
        stop_info['id'] = sample_id
        stop_info_list.append(stop_info)

        for tag, tag_df in feature_dict.items():
            if tag not in df_dict:
                df_dict[tag] = tag_df
            else:
                df_dict[tag] = pd.concat([df_dict[tag], tag_df], axis=0)

    df_dict['StopInfo'] = pd.DataFrame(stop_info_list).set_index('id')

    # Verify result size
    for df in df_dict.values():
        assert len(df) == n_samples, "Result length mismatch"

    return df_dict

def run_BED_on_dataset(freq_file, ref_bg_file, out_dir, alpha_init=4.596, max_epochs=10000, N=50, lr=1, processes=20):
    """Perform BED decomposition on the entire dataset and save results."""
    logger.info(f"Loading frequency matrix: {freq_file}")
    freq_df = pd.read_table(freq_file, header=0, index_col=0)

    logger.info(f"Loading reference background: {ref_bg_file}")
    ref_bg = np.load(ref_bg_file)

    logger.info("Running multi-process BED decomposition...")
    df_dict = run_multi_process_BED(freq_df, ref_bg, processes, alpha_init, max_epochs, N, lr)

    logger.info(f"Saving results to directory: {out_dir}")
    os.makedirs(out_dir, exist_ok=True)
    for key, df in df_dict.items():
        df.to_csv(os.path.join(out_dir, key+'.tsv'), sep='\t', header=True, index=True)
    logger.info("BED decomposition completed.")

def run_BED_on_pilot(freq_file, ref_bg_file, out_dir, alpha_init=4.596, max_epochs=10000, N=50, processes=20):
    "Pilot mode: Tests different learning rates and outputs the average convergence epochs."
    os.makedirs(out_dir, exist_ok=True)

    # The test range can be adjusted by the user.
    # lr_list = [0.001, 0.005, 0.01, 0.1] + [float(round(x, 1)) for x in np.arange(1, 2.1, 0.1)]
    lr_list = [0.001, 0.01, 0.1, 1, 1.5, 2]

    logger.info(f"PILOT MODE STARTED")
    logger.info(f"Loading frequency matrix: {freq_file}")
    freq_df = pd.read_table(freq_file, header=0, index_col=0)

    logger.info(f"Loading reference background: {ref_bg_file}")
    ref_bg = np.load(ref_bg_file)

    logger.info(f"Testing {len(lr_list)} learning rates on {freq_df.shape[0]} samples...")
    logger.info(f"lr candidates: {lr_list}")
    lr_data = []
    for i, lr in enumerate(lr_list, 1):
        logger.info(f"[{i}/{len(lr_list)}] Testing lr = {lr:g}")

        df_dict = run_multi_process_BED(freq_df, ref_bg, processes, alpha_init, max_epochs, N, lr)
        stop_df: pd.DataFrame = df_dict['StopInfo']
        epoch_mean = stop_df['epoch'].mean()
        early_stop_rate = stop_df['early'].mean() * 100

        lr_data.append([lr, epoch_mean, early_stop_rate])
        logger.info(f"mean epochs: {epoch_mean:.1f}, early-stop rate: {early_stop_rate:.1f}%")

        out_path = os.path.join(out_dir, f"StopInfo_lr_{lr:g}.tsv")
        stop_df.to_csv(out_path, sep='\t', header=True, index=True)

    df = pd.DataFrame(lr_data, columns=['lr', 'epoch_mean', 'early_stop_rate_%'])
    out_path = os.path.join(out_dir, "Pilot_lr_search.tsv")
    df.to_csv(out_path, sep='\t', header=True, index=False)

    logger.info(f"Pilot result saved to: {out_path}")
    logger.info("Pilot mode completed.")

def main():
    args = parse_args()

    if args.pilot:
        run_BED_on_pilot(
            freq_file=args.frequency_file,
            ref_bg_file=args.reference_file,
            out_dir=args.output_dir,
            alpha_init=args.alpha_init,
            max_epochs=args.max_epochs,
            N=args.patience,
            processes=args.processes
        )
    else:
        run_BED_on_dataset(
            freq_file=args.frequency_file,
            ref_bg_file=args.reference_file,
            out_dir=args.output_dir,
            alpha_init=args.alpha_init,
            max_epochs=args.max_epochs,
            N=args.patience,
            lr=args.lr,
            processes=args.processes
        )

if __name__ == "__main__":
    main()
