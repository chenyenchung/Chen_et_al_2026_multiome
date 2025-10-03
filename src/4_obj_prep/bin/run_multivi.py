import argparse
import os
from pathlib import Path

# Set commandline arguments
parser = argparse.ArgumentParser(
    prog='scverse Runner',
    description='Run MultiVI on multiome data'
)

parser.add_argument('-f', '--h5mu', required=True)
parser.add_argument('-o', '--outdir', default='multiVI')
parser.add_argument('-k', '--nfactor', type=int, default=None)
parser.add_argument('-i', '--batch', default=None)
parser.add_argument('-p', '--penalty', default='Jeffreys')

args = parser.parse_args()

# Load the H5AD
import muon as mu
import numpy as np
import scvi
import pandas as pd

objpath = Path(args.h5mu)
mdata = mu.read_h5mu(objpath)

try:
    os.mkdir(args.outdir)
except FileExistsError:
    pass

modalities_map = {
    'rna_layer': 'rna',
    'atac_layer': 'atac',
}
# Prepare for training
scvi.model.MULTIVI.setup_mudata(
    mdata,
    batch_key=args.batch,
    modalities = modalities_map,
)
model = scvi.model.MULTIVI(
    mdata,
    n_latent=args.nfactor,
    modality_weights='cell',
    modality_penalty=args.penalty,
    fully_paired=False,
)

model.train(
    max_epochs=1000,
    plan_kwargs={
        'n_epochs_kl_warmup': 400,
        'reduce_lr_on_plateau': True,
        'lr_scheduler_metric': 'reconstruction_loss_validation',
    },
    check_val_every_n_epoch=10,
    adversarial_mixing=False,
    early_stopping=False
)

# Save the model under the present working directory
# in a directory named poissonvi
if args.batch is not None:
    model_dir = os.path.join(args.outdir, f"multiVI_{args.batch}_{args.nfactor}")
else:
    model_dir = os.path.join(args.outdir, f"multiVI_{args.nfactor}_no_int")

model.save(model_dir, overwrite=True)
joint_latent = model.get_latent_representation(modality='joint')
gex_latent = model.get_latent_representation(modality='expression')
atac_latent = model.get_latent_representation(modality='accessibility')
np.savetxt(os.path.join(args.outdir, f"multiVI_k{args.nfactor}_joint.csv"), joint_latent, delimiter=',')
np.savetxt(os.path.join(args.outdir, f"multiVI_k{args.nfactor}_gex.csv"), gex_latent, delimiter=',')
np.savetxt(os.path.join(args.outdir, f"multiVI_k{args.nfactor}_atac.csv"), atac_latent, delimiter=',')

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
v_gex = gex_latent.var(axis=0)
np.savetxt(os.path.join(args.outdir, "latent_var_gex.csv"), v_gex, delimiter=",")
print("Gex dims var>0.5:", (v_gex > 0.5).sum(), " / ", v_gex.size)
print("Gex dims var>0.1:", (v_gex > 0.1).sum(), " / ", v_gex.size)

v_atac = atac_latent.var(axis=0)
np.savetxt(os.path.join(args.outdir, "latent_var_atac.csv"), v_atac, delimiter=",")
print("ATAC dims var>0.5:", (v_atac > 0.5).sum(), " / ", v_atac.size)
print("ATAC dims var>0.1:", (v_atac > 0.1).sum(), " / ", v_atac.size)

v_joint = joint_latent.var(axis=0)
np.savetxt(os.path.join(args.outdir, "latent_var_joint.csv"), v_joint, delimiter=",")
print("Joint dims var>0.5:", (v_joint > 0.5).sum(), " / ", v_joint.size)
print("Joint dims var>0.1:", (v_joint > 0.1).sum(), " / ", v_joint.size)

