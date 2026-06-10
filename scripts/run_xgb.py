#!/usr/bin/env python3
import argparse
from pathlib import Path


def require_deps():
    missing = []
    for name in ("numpy", "pandas", "sklearn", "xgboost"):
        try:
            __import__(name)
        except ImportError:
            missing.append(name)
    if missing:
        raise RuntimeError(
            "Missing Python package(s): "
            + ", ".join(missing)
            + ". Install them in the python_bin environment used by FAST_treatment_ML."
        )


def auc(y, pred):
    from sklearn.metrics import roc_auc_score

    if len(set(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, pred))


def build_repeated_folds(fold_df, train_subjects, np):
    train_subjects = [str(subject_id) for subject_id in train_subjects]
    subject_index = {subject_id: index for index, subject_id in enumerate(train_subjects)}
    repeated_folds = []

    fold_df = fold_df.copy()
    fold_df["SUBJECT_ID"] = fold_df["SUBJECT_ID"].astype(str)
    for repeat_id in sorted(fold_df["REPEAT"].unique()):
        repeat_df = fold_df[fold_df["REPEAT"] == repeat_id]
        if len(repeat_df) != len(train_subjects):
            raise RuntimeError(f"XGB repeat {repeat_id} does not contain every training subject.")
        if set(repeat_df["SUBJECT_ID"]) != set(train_subjects):
            raise RuntimeError(f"XGB repeat {repeat_id} subject IDs do not match subjects.csv.")

        foldid = np.empty(len(train_subjects), dtype=int)
        for row in repeat_df.itertuples(index=False):
            foldid[subject_index[row.SUBJECT_ID]] = int(row.FOLD_ID)

        folds = []
        for fold in sorted(np.unique(foldid)):
            valid_idx = np.where(foldid == fold)[0]
            train_idx = np.where(foldid != fold)[0]
            folds.append((train_idx, valid_idx))
        repeated_folds.append((int(repeat_id), folds))

    if not repeated_folds:
        raise RuntimeError("xgb_folds.csv does not contain any CV repeats.")
    return repeated_folds


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--nthread", type=int, default=1)
    parser.add_argument("--n-trials", type=int, default=0)
    args = parser.parse_args()

    require_deps()

    import numpy as np
    import pandas as pd
    import xgboost as xgb

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    x_train_df = pd.read_csv(data_dir / "xgb_train.csv.gz")
    x_test_df = pd.read_csv(data_dir / "xgb_test.csv.gz")
    subjects = pd.read_csv(data_dir / "subjects.csv")
    xgb_folds_df = pd.read_csv(data_dir / "xgb_folds.csv")
    y_train_df = subjects[subjects["SET"] == "train"].reset_index(drop=True)
    y_test_df = subjects[subjects["SET"] == "test"].reset_index(drop=True)

    x_train = x_train_df.to_numpy(dtype=float)
    x_test = x_test_df.to_numpy(dtype=float)
    y_train = y_train_df["TREATMENT_GROUP"].to_numpy(dtype=int)
    y_test = y_test_df["TREATMENT_GROUP"].to_numpy(dtype=int)

    dtrain = xgb.DMatrix(x_train, label=y_train, feature_names=list(x_train_df.columns))
    dtest = xgb.DMatrix(x_test, label=y_test, feature_names=list(x_train_df.columns))
    repeated_folds = build_repeated_folds(
        xgb_folds_df,
        y_train_df["SUBJECT_ID"].tolist(),
        np,
    )

    params = {
        "objective": "binary:logistic",
        "eval_metric": "auc",
        "max_depth": 2,
        "eta": 0.03,
        "min_child_weight": 5,
        "subsample": 0.8,
        "colsample_bytree": 0.6,
        "lambda": 10.0,
        "alpha": 1.0,
        "nthread": max(1, args.nthread),
        "seed": args.seed,
    }

    def evaluate_repeated_cv(trial_params):
        repeat_aucs = []
        repeat_iterations = []
        for repeat_id, folds in repeated_folds:
            cv = xgb.cv(
                trial_params,
                dtrain,
                num_boost_round=500,
                folds=folds,
                early_stopping_rounds=25,
                verbose_eval=False,
                seed=args.seed + repeat_id,
            )
            best_index = int(cv["test-auc-mean"].idxmax())
            repeat_aucs.append(float(cv["test-auc-mean"].iloc[best_index]))
            repeat_iterations.append(best_index + 1)

        return {
            "mean_auc": float(np.mean(repeat_aucs)),
            "sd_auc": float(np.std(repeat_aucs, ddof=1)) if len(repeat_aucs) > 1 else 0.0,
            "best_iteration": int(round(float(np.median(repeat_iterations)))),
            "repeat_aucs": repeat_aucs,
            "repeat_iterations": repeat_iterations,
        }

    tuning_path = out_dir / "tuning.csv"
    if args.n_trials > 0:
        try:
            import optuna
        except ImportError as exc:
            raise RuntimeError("n_trials > 0 requires the Python package 'optuna'.") from exc

        def objective(trial):
            trial_params = dict(params)
            trial_params.update(
                {
                    "max_depth": trial.suggest_int("max_depth", 1, 4),
                    "eta": trial.suggest_float("eta", 0.005, 0.08, log=True),
                    "min_child_weight": trial.suggest_float("min_child_weight", 2, 25, log=True),
                    "subsample": trial.suggest_float("subsample", 0.55, 0.95),
                    "colsample_bytree": trial.suggest_float("colsample_bytree", 0.25, 0.85),
                    "lambda": trial.suggest_float("lambda", 1, 100, log=True),
                    "alpha": trial.suggest_float("alpha", 0.01, 20, log=True),
                }
            )
            result = evaluate_repeated_cv(trial_params)
            trial.set_user_attr("cv_auc_sd", result["sd_auc"])
            trial.set_user_attr("best_iteration", result["best_iteration"])
            for index, value in enumerate(result["repeat_aucs"], start=1):
                trial.set_user_attr(f"repeat_{index}_auc", value)
            for index, value in enumerate(result["repeat_iterations"], start=1):
                trial.set_user_attr(f"repeat_{index}_best_iteration", value)
            return result["mean_auc"]

        sampler = optuna.samplers.TPESampler(seed=args.seed)
        study = optuna.create_study(direction="maximize", sampler=sampler)
        study.optimize(objective, n_trials=args.n_trials)
        params.update(study.best_params)
        best_iteration = int(study.best_trial.user_attrs["best_iteration"])
        cv_auc = float(study.best_value)
        cv_auc_sd = float(study.best_trial.user_attrs["cv_auc_sd"])

        tuning_rows = []
        for trial in study.trials:
            row = {
                "TRIAL": trial.number,
                "STATE": trial.state.name,
                "CV_AUC_MEAN": trial.value,
                "CV_AUC_SD": trial.user_attrs.get("cv_auc_sd"),
                "BEST_ITERATION_MEDIAN": trial.user_attrs.get("best_iteration"),
            }
            for key, value in sorted(trial.user_attrs.items()):
                if key.startswith("repeat_"):
                    row[key.upper()] = value
            row.update({f"PARAM_{key.upper()}": value for key, value in trial.params.items()})
            tuning_rows.append(row)
        pd.DataFrame(tuning_rows).to_csv(tuning_path, index=False)
    else:
        result = evaluate_repeated_cv(params)
        best_iteration = result["best_iteration"]
        cv_auc = result["mean_auc"]
        cv_auc_sd = result["sd_auc"]
        tuning_row = {
            "TRIAL": 0,
            "STATE": "FIXED_PARAMETERS",
            "CV_AUC_MEAN": cv_auc,
            "CV_AUC_SD": cv_auc_sd,
            "BEST_ITERATION_MEDIAN": best_iteration,
        }
        for index, value in enumerate(result["repeat_aucs"], start=1):
            tuning_row[f"REPEAT_{index}_AUC"] = value
        for index, value in enumerate(result["repeat_iterations"], start=1):
            tuning_row[f"REPEAT_{index}_BEST_ITERATION"] = value
        tuning_row.update({f"PARAM_{key.upper()}": value for key, value in params.items()})
        pd.DataFrame([tuning_row]).to_csv(tuning_path, index=False)

    model = xgb.train(params, dtrain, num_boost_round=best_iteration, verbose_eval=False)
    train_pred = model.predict(dtrain)
    test_pred = model.predict(dtest)

    metrics_row = {
        "CV_AUC": cv_auc,
        "CV_AUC_SD": cv_auc_sd,
        "CV_REPEATS": len(repeated_folds),
        "TEST_AUC": auc(y_test, test_pred),
        "INSAMPLE_AUC": auc(y_train, train_pred),
        "BEST_ITERATION": best_iteration,
        "N_FEATURES": x_train.shape[1],
    }
    metrics_row.update({f"PARAM_{key.upper()}": value for key, value in params.items()})
    metrics = pd.DataFrame([metrics_row])
    metrics.to_csv(out_dir / "metrics.csv", index=False)

    predictions = pd.concat(
        [
            pd.DataFrame(
                {
                    "SET": "train",
                    "SUBJECT_ID": y_train_df["SUBJECT_ID"],
                    "FU": y_train_df["FU"],
                    "TREATMENT_GROUP": y_train,
                    "PREDICTED_PROB": train_pred,
                }
            ),
            pd.DataFrame(
                {
                    "SET": "test",
                    "SUBJECT_ID": y_test_df["SUBJECT_ID"],
                    "FU": y_test_df["FU"],
                    "TREATMENT_GROUP": y_test,
                    "PREDICTED_PROB": test_pred,
                }
            ),
        ],
        ignore_index=True,
    )
    predictions.to_csv(out_dir / "predictions.csv", index=False)

    importance = model.get_score(importance_type="gain")
    importance_df = pd.DataFrame(
        [
            {
                "FEATURE_NAME": k,
                "GAIN": v,
                "FEATURE_TYPE": "covariate" if k.startswith("covariate::") else "omics",
            }
            for k, v in importance.items()
        ],
        columns=["FEATURE_NAME", "GAIN", "FEATURE_TYPE"],
    ).sort_values("GAIN", ascending=False)
    importance_df.to_csv(out_dir / "importance.csv", index=False)

    model.save_model(str(out_dir / "model.json"))


if __name__ == "__main__":
    main()
