# LIUQE-ML a 10kHz LIUQE surrogate model

This is the repository for the `liuqe-ml` (or `liuqe-nn` in some places) code. 
`liuqe-ml` is a matlab super fast (30~70 microseconds) neural network surrogate for the equilibrium reconstruction code `liuqe`. It is meant for fast shape control, as an alternative to the `LIUQE-RT` block.

Here there are some

## Test the pretrained models
To build the surrogate mex function (from `liuqenn_inference.m`) do the following steps (after cloning the repo, in matlab):
- convert the pretrained weights from `.json` to `.mat`
    ```matlab
    run('pretrained/json2mat.m'); 
    ```
- compile the mex, this will build the main code and another model for plasma current density `Iy`:
  ```matlab
    build_liuqenn_inference_mex
    ```
- you should now see `liuqenn_inference_mex.mexa64` and `liuqe_Iy_nn_inference.m` mex functions

To run a speed test :
```matlab
    speed_test
```
> note: if you did not run `prepare_dataset.ipynb` before, dont trust the accuracy, only speed. 

To run a comparison with `lih`:
```matlab
    lih_comparison
```

## Training
Dataset generation and preparation steps:
> Note: Training and dataset preparations requires a lot of tuning parameters (read the code to understand how to change them).

> Note: the train/dataset code is meant to run fast locally (see `LOCAL` var) and slowly on a `slurm` equipped server. 

- to generate the raw dataset (in a `lac` server with `liuqe` and `mdsplus`) run:
    ```matlab
    create_ds
    ```
- to convert the raw dataset in a format suitable for training run the jupyter notebook `prepare_dataset.ipynb`. This will also create a demo dataset that can be tested using `speed_test`.
- to train a model: `train.ipynb`