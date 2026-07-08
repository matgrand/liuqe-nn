# LIUQE-NN: a 10kHz LIUQE surrogate model

This is the repository for the `liuqe-nn` (or `liuqe-ml`) code. 
`liuqe-nn` is a matlab super fast (20~70 microseconds cpu inference) neural network surrogate for the equilibrium reconstruction code `liuqe`. It is meant for fast shape control, as an alternative to the `LIUQE-RT` block.

## Test the pretrained models
- First of all build the mex functions: to build the surrogate mex function (from `liuqenn_inference.m` and `liuqe_Iy_nn_inference`) do the following steps (after cloning the repo (and switching to branch `paper`), in `matlab` (use `matlab960` if inside `lac`)):
  - convert the pretrained weights from `.json` to `.mat`:
    ```matlab
    run('pretrained/json2mat.m'); 
    ```
  - compile the mex, this will build the main code and another function for plasma current density `Iy`:
    ```matlab
    build_liuqenn_inference_mex
    ```
  - you should now see `liuqenn_inference_mex.mexa64` and `liuqe_Iy_nn_inference.mexa64` mex functions
  - 
- To run an inference speed test:
    ```matlab
    speed_test
    ```
    > note: if you did not run `prepare_dataset.ipynb` before, dont trust the accuracy, only speed. 

- For a quick comparison with `lih` run:
    ```matlab
    lih_comparison
    ```

- To see how to run an inference on a standard `meq` `LX` structure see `liuqe_nn_LX.m`

- For a full comparison of the surrogate performance on real-time generated data run:
    ```matlab
    realtime_analysis
    ```
    This will compare the surrogate against `liuqe`, `liuqe-rt` and `lih` (using `liuqe` as ground-truth)

## Dataset generation and Training
> Note: Training and dataset preparations requires a lot of tuning parameters (read the code to understand how to change them).


- to generate the raw dataset (in a `lac` server with `liuqe` and `mdsplus`) run the script `create_ds.m`. This script goes through the "saved" `liuqe` snapshots of many shots and performs basic filtering using time vector and plasma current to keep only "good" equilibria. It uses the script `ds_prep_save_shot.m` to filter out `nan`s and `inf`s and saves the dataset in `.mat` format. 
    > "saved" means that we do not re-run `liuqe`, but we use the after-shot calculated one saved in `mdsplus`, this drastically reduces dataset generation time.

- to convert the raw dataset in a format suitable for training run the jupyter notebook `prepare_dataset.ipynb`. This will also create a demo dataset that can be tested using `speed_test`.
- to train a model: `train.ipynb`
> Note: the train/dataset code is meant to run fast (reduced epochs, dataset dimension, ...) locally (see `LOCAL` var) and slow  (full/everything maxxed) on a `slurm` equipped server. 













