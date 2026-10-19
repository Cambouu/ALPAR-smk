"""Train XGBoost with fixed hyperparameters on existing pivoted feature tables.

Checks whether the supervisor's best parameters (eta 0.01, 500 trees, depth 3)
change our results. Each (feature table, split, parameter set) gets one fit on the
train strains and is scored on the test strains. No grid search, no selection.

Splits:
  datasail  the pipeline's train.txt/test.txt
  holdout   stratified random 80/20 split (seed 42) over every strain labelled 0/1
            for the antibiotic, like ALPAR's original resampling "holdout"

Strains in a split but missing from a table get zero features, as in the baseline.
Labels: 1 = resistant, 0 = susceptible.
"""
import argparse
import csv
import time
from pathlib import Path

import numpy as np
import sklearn.metrics as skm
import xgboost as xgb
from sklearn.model_selection import train_test_split

PARAM_SETS = {
    # From the supervisor's best_result.txt files.
    "sup_ampicillin": dict(max_depth=3, min_child_weight=5, subsample=0.8, colsample_bytree=0.8, eta=0.01, n_estimators=500),
    "sup_ciprofloxacin": dict(max_depth=3, min_child_weight=1, subsample=0.8, colsample_bytree=0.8, eta=0.01, n_estimators=500),
}


def load_labels(phenotype_table, antibiotic):
    with open(phenotype_table) as fh:
        reader = csv.reader(fh, delimiter="\t")
        column = next(reader).index(antibiotic)
        return {row[0]: int(row[column]) for row in reader if row[column] in ("0", "1")}


def load_table(path):
    """Wide table: header 'hash<TAB>feature...', one row per strain.

    Lines are split on LF only: old GPA cells end in a stray CR ("1 + CR"), which
    text mode would treat as a line break.
    """
    with open(path, newline="\n") as fh:
        n_features = len(fh.readline().rstrip("\r\n").split("\t")) - 1
        n_rows = sum(1 for _ in fh)
    matrix = np.zeros((n_rows, n_features), dtype=np.float32)
    strains = {}
    with open(path, newline="\n") as fh:
        fh.readline()
        for i, line in enumerate(fh):
            cells = line.rstrip("\r\n").split("\t")
            strains[cells[0]] = i
            # float() ignores the stray CR in old GPA cells ("1\r"); empty means 0.
            matrix[i] = np.fromiter(
                (float(c) if c.strip() else 0.0 for c in cells[1:]),
                dtype=np.float32, count=n_features,
            )
    return strains, matrix


def rows_for(strains, index, matrix):
    out = np.zeros((len(strains), matrix.shape[1]), dtype=np.float32)
    for i, s in enumerate(strains):
        if s in index:
            out[i] = matrix[index[s]]
    return out


def write_result(path, y_test, y_hat, params):
    with open(path, "w") as fh:
        fh.write(f"Parameters: {params}\n")
        fh.write(f"Accuracy score: {skm.accuracy_score(y_test, y_hat)}\n")
        fh.write(f"Balanced Accuracy score: {skm.balanced_accuracy_score(y_test, y_hat)}\n")
        fh.write(f"Brier score loss: {skm.brier_score_loss(y_test, y_hat)}\n")
        for avg in ("macro", "micro", "weighted", "binary"):
            fh.write(f"F1 score {avg}: {skm.f1_score(y_test, y_hat, average=avg)}\n")
        fh.write(f"Precision score: {skm.precision_score(y_test, y_hat)}\n")
        fh.write(f"Recall score: {skm.recall_score(y_test, y_hat)}\n")
        fh.write(f"Confusion matrix: {skm.confusion_matrix(y_test, y_hat).tolist()}\n")
        fh.write(f"ROC AUC Score: {skm.roc_auc_score(y_test, y_hat)}\n")
        fh.write(f"Jaccard score: {skm.jaccard_score(y_test, y_hat)}\n")
        fh.write(f"Matthews correlation coefficient: {skm.matthews_corrcoef(y_test, y_hat)}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--antibiotic", required=True)
    ap.add_argument("--phenotype-table", required=True)
    ap.add_argument("--datasail-dir", required=True, help="contains <abx>/train.txt and test.txt")
    ap.add_argument("--table", action="append", required=True, help="NAME=PATH, repeatable")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    abx = args.antibiotic
    outdir = Path(args.outdir) / abx
    outdir.mkdir(parents=True, exist_ok=True)
    labels = load_labels(args.phenotype_table, abx)

    ds = Path(args.datasail_dir) / abx
    splits = {"datasail": (ds.joinpath("train.txt").read_text().split(), ds.joinpath("test.txt").read_text().split())}
    labelled = list(labels)
    tr, te = train_test_split(labelled, test_size=0.2, random_state=args.seed,
                              stratify=[labels[s] for s in labelled])
    splits["holdout"] = (tr, te)
    for name, (tr, te) in splits.items():
        missing = [s for s in tr + te if s not in labels]
        if missing:
            raise SystemExit(f"{name}: {len(missing)} split strains have no 0/1 label")
        print(f"split {name}: train {len(tr)} (R={sum(labels[s] for s in tr)}), "
              f"test {len(te)} (R={sum(labels[s] for s in te)})", flush=True)

    summary = outdir / "summary.tsv"
    if not summary.exists():
        summary.write_text("table\tsplit\tparams\tn_features\tMCC\tBA\tTN\tFP\tFN\tTP\tfit_s\n")

    for spec in args.table:
        table_name, table_path = spec.split("=", 1)
        t0 = time.time()
        index, matrix = load_table(table_path)
        print(f"{table_name}: {matrix.shape} loaded in {time.time() - t0:.0f}s", flush=True)
        for split_name, (tr, te) in splits.items():
            absent = sum(s not in index for s in tr + te)
            X_tr, X_te = rows_for(tr, index, matrix), rows_for(te, index, matrix)
            y_tr = np.array([labels[s] for s in tr])
            y_te = np.array([labels[s] for s in te])
            for param_name, p in PARAM_SETS.items():
                t1 = time.time()
                model = xgb.XGBClassifier(
                    objective="binary:logistic", eval_metric="logloss", seed=args.seed, n_jobs=args.threads,
                    max_depth=p["max_depth"], min_child_weight=p["min_child_weight"], subsample=p["subsample"],
                    colsample_bytree=p["colsample_bytree"], learning_rate=p["eta"], n_estimators=p["n_estimators"],
                ).fit(X_tr, y_tr)
                y_hat = model.predict(X_te)
                fit_s = time.time() - t1
                write_result(outdir / f"{table_name}__{split_name}__{param_name}_Result.txt", y_te, y_hat, p)
                tn, fp, fn, tp = skm.confusion_matrix(y_te, y_hat).ravel()
                mcc = skm.matthews_corrcoef(y_te, y_hat)
                ba = skm.balanced_accuracy_score(y_te, y_hat)
                with open(summary, "a") as fh:
                    fh.write(f"{table_name}\t{split_name}\t{param_name}\t{matrix.shape[1]}\t{mcc:.4f}\t{ba:.4f}\t{tn}\t{fp}\t{fn}\t{tp}\t{fit_s:.0f}\n")
                print(f"{table_name} {split_name} {param_name}: MCC {mcc:.4f} BA {ba:.4f} "
                      f"[[{tn} {fp}] [{fn} {tp}]] zero-filled={absent} fit {fit_s:.0f}s", flush=True)
            del X_tr, X_te
        del matrix, index


if __name__ == "__main__":
    main()
