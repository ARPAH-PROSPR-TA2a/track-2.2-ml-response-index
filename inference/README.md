# FAST Treatment ML Inference

This directory will contain the user-facing inference pipeline.

The inference surface will consume artifacts emitted by training:

- `manifest.json`
- `reports/preprocessing.csv`
- `FU*/xgb/model.json`
- `FU*/enet/weights.csv`

The first implementation target is an equivalence test: score the training
fixture through the inference path and reproduce the predictions written by the
training pipeline. After that test anchors the artifact contract, this README
will document the public inference API.
