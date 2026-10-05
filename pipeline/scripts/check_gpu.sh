set -eo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
init_conda
conda activate binderflow
load_tools
set -u
conda run --no-capture-output -n binder-rfd-wsl python -c 'import torch,dgl,rfdiffusion,se3_transformer; assert torch.cuda.is_available(), "RFdiffusion CUDA unavailable"; g=dgl.graph(([0,1],[1,0]),device="cuda"); g.ndata["x"]=torch.ones((2,1),device="cuda"); g.update_all(dgl.function.copy_u("x","m"),dgl.function.sum("m","h")); torch.cuda.synchronize(); assert torch.equal(g.ndata["h"],torch.ones((2,1),device="cuda")); print("RFdiffusion:",torch.__version__,dgl.__version__,torch.cuda.get_device_name())'
conda run --no-capture-output -n binder-mpnn-wsl python -c 'import torch; assert torch.cuda.is_available(), "MPNN CUDA unavailable"; a=torch.ones((32,32),device="cuda"); assert (a@a).sum().item()==32768; print("MPNN:",torch.__version__,torch.cuda.get_device_name())'
"$(dirname "$COLABFOLD_BATCH")/python" -c 'import jax,jax.numpy as jnp; import colabfold.batch; gpu=jax.devices("gpu")[0]; a=jax.device_put(jnp.ones((32,32)),gpu); assert float((a@a).sum().block_until_ready())==32768; print("ColabFold/JAX:",jax.__version__,gpu)'
"$COLABFOLD_BATCH" --help >/dev/null
printf '%s\n' 'GPU import and arithmetic checks passed. Protein inference remains a separate test.'
