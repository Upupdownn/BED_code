#!/usr/bin/env python3
"""
Train and validate SVM model on features.
Supports train (cross-validation) and validate (load model) modes.
"""

import argparse
import logging
import pandas as pd
import numpy as np
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
from sklearn.model_selection import train_test_split, StratifiedKFold
import joblib
from pathlib import Path
import os


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description="Train and validate SVM model on features.")
    parser.add_argument("feature_file", help="Input TSV file with sample features (samples x features)")
    parser.add_argument("label_file", help="Label TSV file: sample_id, label (0=control, 1=cancer)")
    parser.add_argument("output_file", help="Output file for SVM prediction scores")
    parser.add_argument("--mode", choices=['train', 'validate'], default='train',
                        help="Mode: train (cross-validation), validate (load model)")
    parser.add_argument("--model_dir", default=None, help="Path to model dir (save for train, load for validate)")
    parser.add_argument("--n_fold", type=int, default=5, help="Number of CV folds [default: 5]")
    parser.add_argument("--n_repeat", type=int, default=10, help="Number of CV repeats [default: 10]")
    parser.add_argument("--random_state", type=int, default=42, help="Random seed")
    return parser.parse_args()

def load_data(feature_file, label_file):
    """Load features and labels, align by sample ID."""
    features = pd.read_table(feature_file, header=0, index_col=0)
    labels = pd.read_table(label_file, header=0, index_col=0)
    common_ids = features.index.intersection(labels.index)
    if len(common_ids) == 0:
        raise ValueError("No common sample IDs between features and labels")
    X = features.loc[common_ids].values
    y = labels.loc[common_ids].iloc[:, 0].to_numpy()
    sample_ids = common_ids
    logger.info(f"Loaded {len(common_ids)} samples with {X.shape[1]} features")
    return X, y, sample_ids

def cross_validation(X, y, n_fold=10, n_repeat=10, random_state=42):
    n_prob = np.zeros((n_repeat, len(y)), dtype=np.float64)
    model_list = []
    for i in range(n_repeat):
        kfold = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=i)
        for j, (train_idx, test_idx) in enumerate(kfold.split(X, y)):
            x_train, y_train, x_test = X[train_idx], y[train_idx], X[test_idx]

            model = SVC(kernel='linear', C=1, probability=True, random_state=random_state)
            model.fit(x_train, y_train)
        
            n_prob[i][test_idx] = model.predict_proba(x_test)[:, 1]
            model_list.append(model)
    y_prob = np.mean(n_prob, axis=0)
    auc = roc_auc_score(y, y_prob)
    logger.info(f"CV AUC: {auc:.4f}")

    return auc, y_prob, model_list

def validate_model(model_list, X, y):
    """Evaluate model performance."""
    model_probs = []    # 一行为一个model的结果
    for model in model_list:
        model_probs.append(model.predict_proba(X)[:, 1])
    
    y_prob = np.mean(model_probs, axis=0)
    auc = roc_auc_score(y, y_prob)
    logger.info(f"validate AUC: {auc:.4f}")
    return auc, y_prob

def main():
    args = parse_args()

    if args.model_dir is not None:
        Path(args.model_dir).mkdir(parents=True, exist_ok=True)

    logger.info("Loading data...")
    X, y, sample_ids = load_data(args.feature_file, args.label_file)

    if args.mode == 'train':
        logger.info(f"train with Cross-validation: {args.n_fold}-fold, repeat {args.n_repeat} times")
        auc, scores, model_list = cross_validation(X, y, args.n_fold, args.n_repeat, args.random_state)

        if args.model_dir is not None:
            for i, model in enumerate(model_list):
                joblib.dump(model, os.path.join(args.model_dir, f'svm_model_{i}.pkl'))
            logger.info(f"Model saved in {args.model_dir}.")
        df = pd.DataFrame({'id': sample_ids, 'score': scores})
        df.to_csv(args.output_file, header=True, index=False, sep='\t')
        logger.info(f"SVM prediction scores saved in {args.output_file}.")
    elif args.mode == 'validate':
        if not args.model_dir:
            raise ValueError("--model_dir required for validate mode")
        logger.info(f"Validate mode: loading model from {args.model_dir}")
        model_list = []
        for fname in os.listdir(args.model_dir):
            model = joblib.load(os.path.join(args.model_dir, fname))
            model_list.append(model)
        logger.info(f"Validate mode: load {len(model_list)} models")

        auc, scores = validate_model(model_list, X, y)
        df = pd.DataFrame({'id': sample_ids, 'score': scores})
        df.to_csv(args.output_file, header=True, index=False, sep='\t')
        logger.info(f"SVM prediction scores saved in {args.output_file}.")


if __name__ == "__main__":
    main()

