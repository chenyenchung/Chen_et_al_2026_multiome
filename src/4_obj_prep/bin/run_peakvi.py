import argparse
import os
from pathlib import Path

# Set commandline arguments
parser = argparse.ArgumentParser(
    prog='scverse Runner',
    description='Run PeakVI on ATAC counts'
)

parser.add_argument('-f', '--h5mu', required=True)
parser.add_argument('-o', '--outdir', default='peakVI')
parser.add_argument('-k', '--nfactor', type=int, default=None)
parser.add_argument('-i', '--batch', default=None)

args = parser.parse_args()

# Load the H5AD
import muon as mu
import numpy as np
import scvi
import pandas as pd

objpath = Path(args.h5mu)
mdata = mu.read_h5mu(objpath)
mdata.mod['atac'].obs = mdata.obs

pwd = objpath.parent.absolute()
os.chdir(pwd)
try:
    os.mkdir(args.outdir)
except FileExistsError:
    pass

pwd = pwd.joinpath(args.outdir)

# Prepare for training
scvi.model.PEAKVI.setup_anndata(
    mdata.mod['atac'],
    batch_key=args.batch
)
model = scvi.model.PEAKVI(
    mdata.mod['atac'],
    n_latent=args.nfactor
)

model.train(
    max_epochs=1000,
    plan_kwargs={
        'n_epochs_kl_warmup': 400,
        'reduce_lr_on_plateau': True,
        'lr_scheduler_metric': 'reconstruction_loss_validation',
    },
    check_val_every_n_epoch=10,
    early_stopping=False
)

# Save the model under the present working directory
if args.batch is not None:
    model_dir = os.path.join(args.outdir, f"peakVI_{args.batch}_{args.nfactor}")
else:
    model_dir = os.path.join(args.outdir, f"peakVI_{args.nfactor}_no_int")

model.save(model_dir, overwrite=True)
latent = model.get_latent_representation()
np.savetxt(os.path.join(args.outdir, f"peakVI_k{args.nfactor}.csv"), latent, delimiter=',')

# Save loss
def flatten_history(history) -> pd.DataFrame:
    cols = {}
    for k, v in history.items():
        if v is None:
            continue
        # v can be a pd.DataFrame, pd.Series, numpy array, or list
        if isinstance(v, pd.DataFrame):
            # if it’s 1-column DF, take that column; else keep as-is with multi-columns
            if v.shape[1] == 1:
                cols[k] = v.iloc[:, 0]
            else:
                # prefix multi-col metrics
                for c in v.columns:
                    cols[f"{k}.{c}"] = v[c]
        elif isinstance(v, pd.Series):
            cols[k] = v
        else:
            cols[k] = pd.Series(v)

    df = pd.concat(cols, axis=1)

    # normalize index name (epoch/step) and reset to a column
    if df.index.name is None:
        df.index.name = "index"
    df = df.reset_index()
    return df

hist = flatten_history(model.history)
hist.to_csv(os.path.join(args.outdir, "history.csv"), index=False)

# Estimate latent dim usage
v = latent.var(axis=0)
np.savetxt(os.path.join(args.outdir, "latent_var.csv"), v, delimiter=",")
print("dims var>0.5:", (v > 0.5).sum(), " / ", v.size)
print("dims var>0.1:", (v > 0.1).sum(), " / ", v.size)

# Save normalized accessibility
acc = model.get_normalized_accessibility(transform_batch='stf_2')
np.savetxt(os.path.join(args.outdir, f"peakVI_normacc_k{args.nfactor}.csv"), acc, delimiter=',')
