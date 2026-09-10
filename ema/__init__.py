"""PeakATail package init.

IMPORTANT: this module sets BLAS / OpenMP thread caps BEFORE any submodule
(or numpy/scipy/scanpy) gets imported.  Without this, the multiprocessing
``spawn`` Pool used by peak calling leaves the parent process's BLAS thread
pool in a state where downstream ``scanpy``/``leidenalg`` calls deadlock on
a futex inside OpenBLAS — observed as a hang at the clustering step on
datasets larger than a few hundred cells.

We default to single-threaded BLAS in the parent because (a) all heavy
matrix work is already parallelised at the multiprocessing level (one
worker per dataset / per tile), and (b) the parent process orchestrates
spawn workers — over-subscribing BLAS threads in the parent only competes
with worker CPU.  Users who want multi-threaded BLAS for a specific call
can override by exporting these vars themselves before invoking ``ema``.
"""
import os as _os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_var, "1")
