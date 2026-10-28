"""Exact compressed K55-K64 local-profile screen for F2 2+1+1.

This extends the valuation DP from the K52-K54 receipt. It is a finite
profile screen only; the result cannot be promoted to whole-branch exclusion.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pyopencl as cl

from f2_75_125_k17_nonzero_first_slice import INFINITY

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_k52_k54_opencl_dp_20261024.json"
OUT = HERE / "f2_75_125_211_k55_k64_opencl_dp_20261024.json"
FIRST, LAST, MAX_N = 55, 64, 64
GRID_COUNT = 51_870
INF = 1 << 29

KERNEL = r"""
__kernel void profile_dp(
    __global const int *profiles, const int nprofiles,
    const int infinity_order, const int inf,
    __global int *main_min, __global int *main_ways,
    __global uchar *main_small, __global int *corr_min,
    __global int *corr_argmin)
{
    const int p = get_global_id(0);
    if (p >= nprofiles) return;
    int dp[65], ways[65];
    uchar small[65];
    for (int s=0; s<=64; ++s) { dp[s]=inf; ways[s]=0; small[s]=0; }
    dp[0]=0; ways[0]=1; small[0]=1;
    for (int i=1; i<=64; ++i) {
        int cost=-6;
        if (i<=5) {
            const int oi=profiles[p*5+i-1];
            if (oi>=infinity_order) continue;
            cost=oi-6;
        }
        for (int s=i; s<=64; ++s) {
            const int prev=s-i;
            if (ways[prev]==0) continue;
            const int candidate=dp[prev]+cost;
            if (candidate<dp[s]) {
                dp[s]=candidate; ways[s]=ways[prev];
                small[s]=(uchar)((i<=5) && small[prev]);
            } else if (candidate==dp[s]) {
                ways[s]=min(2,ways[s]+ways[prev]); small[s]=0;
            }
        }
    }
    for (int b=0; b<10; ++b) {
        const int n=55+b, out=p*10+b;
        main_min[out]=10+dp[n]; main_ways[out]=ways[n];
        main_small[out]=small[n];
        int best=inf, best_r=-1;
        for (int r=1; r<=n/5; ++r) {
            const int total=n-5*r;
            if (r==5) continue;
            int value;
            if (r==2) value=0;
            else {
                if (dp[total]==inf) continue;
                value=2*(5-r)+dp[total];
            }
            if (value<best) { best=value; best_r=r; }
        }
        corr_min[out]=best; corr_argmin[out]=best_r;
    }
}
"""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cpu_dp(profiles: np.ndarray):
    count = len(profiles)
    dp = np.full((count, MAX_N + 1), INF, dtype=np.int32)
    ways = np.zeros((count, MAX_N + 1), dtype=np.uint8)
    small = np.zeros((count, MAX_N + 1), dtype=np.bool_)
    dp[:, 0] = 0; ways[:, 0] = 1; small[:, 0] = True
    for i in range(1, MAX_N + 1):
        if i <= 5:
            orders = profiles[:, i-1]
            allowed = orders < INFINITY
            costs = orders.astype(np.int32) - 6
        else:
            allowed = np.ones(count, dtype=np.bool_)
            costs = np.full(count, -6, dtype=np.int32)
        for s in range(i, MAX_N + 1):
            prev_ways = ways[:, s-i]
            valid = allowed & (prev_ways != 0)
            candidate = dp[:, s-i] + costs
            improve = valid & (candidate < dp[:, s])
            tie = valid & (candidate == dp[:, s])
            dp[improve, s] = candidate[improve]
            ways[improve, s] = prev_ways[improve]
            small[improve, s] = (i <= 5) & small[improve, s-i]
            ways[tie, s] = np.minimum(
                2, ways[tie, s].astype(np.int16) + prev_ways[tie].astype(np.int16)
            ).astype(np.uint8)
            small[tie, s] = False
    shape = (count, LAST-FIRST+1)
    mm = np.empty(shape, dtype=np.int32)
    mw = np.empty(shape, dtype=np.int32)
    ms = np.empty(shape, dtype=np.uint8)
    cm = np.full(shape, INF, dtype=np.int32)
    cr = np.full(shape, -1, dtype=np.int32)
    for b, n in enumerate(range(FIRST, LAST+1)):
        mm[:, b] = 10 + dp[:, n]
        mw[:, b] = ways[:, n]
        ms[:, b] = small[:, n]
        for r in range(1, n//5+1):
            total = n-5*r
            if r == 5: continue
            if r == 2:
                values = np.zeros(count, dtype=np.int32)
                valid = np.ones(count, dtype=np.bool_)
            else:
                values = 2*(5-r)+dp[:, total]
                valid = dp[:, total] != INF
            take = valid & (values < cm[:, b])
            cm[take, b] = values[take]
            cr[take, b] = r
    return mm, mw, ms, cm, cr


def partitions(n: int) -> int:
    ways = [0]*(n+1); ways[0]=1
    for part in range(1,n+1):
        for total in range(part,n+1): ways[total]+=ways[total-part]
    return ways[n]


def main() -> None:
    started = time.monotonic()
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    if prior.get("last_complete_band") != FIRST-1:
        raise AssertionError("input receipt is not complete through K54")
    raw = prior["survivor_profiles"]
    if len(raw) != prior.get("survivor_profiles_after_last_complete_band"):
        raise AssertionError("K54 survivor count mismatch")
    if prior.get("finite_profile_grid_count") != GRID_COUNT:
        raise AssertionError("unexpected finite profile grid")
    profiles = np.asarray([
        [INFINITY if v == "infinity" else int(v) for v in row] for row in raw
    ], dtype=np.int32)
    nprofiles = len(profiles)

    devices = [d for p in cl.get_platforms() for d in p.get_devices()
               if d.type & cl.device_type.GPU]
    if not devices: raise RuntimeError("no OpenCL GPU is available")
    device = next((d for d in devices if "Iris" in d.name), devices[0])
    context = cl.Context([device])
    queue = cl.CommandQueue(context, properties=cl.command_queue_properties.PROFILING_ENABLE)
    program = cl.Program(context, KERNEL).build()
    shape = (nprofiles, LAST-FIRST+1)
    gm = np.empty(shape, np.int32); gw = np.empty(shape, np.int32)
    gs = np.empty(shape, np.uint8); gc = np.empty(shape, np.int32)
    gr = np.empty(shape, np.int32)
    pbuf = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=profiles)
    host = (gm, gw, gs, gc, gr)
    buffers = [cl.Buffer(context, cl.mem_flags.WRITE_ONLY, a.nbytes) for a in host]
    event = program.profile_dp(
        queue, (((nprofiles+127)//128)*128,), (128,), pbuf,
        np.int32(nprofiles), np.int32(INFINITY), np.int32(INF), *buffers,
    )
    event.wait()
    for a, b in zip(host, buffers): cl.enqueue_copy(queue, a, b)
    queue.finish()
    for b in buffers: b.release()
    pbuf.release()

    cpu = cpu_dp(profiles)
    for label, a, b in zip(("min", "ties", "support", "correction", "argmin"), host, cpu):
        if not np.array_equal(a, b):
            ij = np.argwhere(a != b)[0].tolist()
            raise AssertionError(f"GPU/CPU mismatch in {label} at {ij}")
    flags = (gm < 0) & (gw == 1) & (gs != 0) & (gc != INF) & (gm < gc)

    active = np.ones(nprofiles, dtype=np.bool_)
    old_excluded = int(prior["cumulative_profiles_excluded_K9_K54"])
    cumulative_new = 0
    counts = dict(prior["first_excluding_band_counts"])
    by_u2 = dict(prior["first_excluding_band_by_ord_U2"])
    term_counts = dict(prior["term_counts_by_band"])
    records = {}; examples = {}
    for b, band in enumerate(range(FIRST, LAST+1)):
        entering = active.copy()
        indices = np.flatnonzero(entering & flags[:, b])
        counts[str(band)] = int(len(indices))
        cumulative_new += int(len(indices))
        hist = {}; ex = []
        for idx in indices.tolist():
            key = "infinity" if raw[idx][1] == "infinity" else str(raw[idx][1])
            hist[key] = hist.get(key, 0)+1
            if len(ex)<8:
                ex.append({"orders_U1_to_U5":raw[idx],"main_order":int(gm[idx,b]),
                           "correction_order":int(gc[idx,b]),"correction_r":int(gr[idx,b])})
        by_u2[str(band)] = hist
        examples[str(band)] = ex
        term_counts[str(band)] = {
            "main_F_5_3": partitions(band),
            "corrections_by_r": {
                str(r): (1 if r==2 else 0 if r==5 else partitions(band-5*r))
                for r in range(1,band//5+1)
            },
        }
        records[str(band)] = {
            "profiles_entering_band":int(np.count_nonzero(entering)),
            "profiles_excluded":int(len(indices)),
            "profiles_remaining":int(np.count_nonzero(entering)-len(indices)),
            "main_min_order_range":[int(gm[entering,b].min()),int(gm[entering,b].max())],
            "unique_main_minimizers":int(np.count_nonzero(gw[entering,b]==1)),
            "unique_small_support_minimizers":int(np.count_nonzero((gw[entering,b]==1)&(gs[entering,b]!=0))),
            "profiles_with_finite_correction":int(np.count_nonzero(gc[entering,b]!=INF)),
        }
        active[indices] = False

    survivors = [raw[i] for i in np.flatnonzero(active).tolist()]
    total = old_excluded+cumulative_new
    result = {
        "scope":"exact compressed K55-K64 finite local valuation screen for F2 2+1+1 profiles",
        "field":"integer dynamic programming; exact profile orders and partition minima",
        "source_receipts":{SOURCE.name:sha256(SOURCE)},
        "script_sha256":sha256(Path(__file__)),
        "gpu":{"platform":device.platform.name,"device":device.name,"driver":device.driver_version},
        "gpu_cpu_parity":True,"gpu_cpu_profile_results_compared":int(nprofiles*(LAST-FIRST+1)),
        "gpu_kernel_seconds":float(event.profile.end-event.profile.start)/1e9,
        "wall_seconds":round(time.monotonic()-started,3),
        "finite_profile_grid_count":GRID_COUNT,"profiles_entering_K55":nprofiles,
        "last_complete_band":LAST,"first_excluding_band_counts":counts,
        "first_excluding_band_by_ord_U2":by_u2,"term_counts_by_band":term_counts,
        "band_records":records,"new_profiles_excluded_K55_K64":cumulative_new,
        "cumulative_profiles_excluded_K9_K64":total,
        "cumulative_percent_of_finite_profile_grid":round(100*total/GRID_COUNT,6),
        "survivor_profiles_after_last_complete_band":len(survivors),
        "survivor_profiles":survivors,"examples_K55_K64":examples,
        "status":"finite local-profile pruning only; whole 2+1+1 branch exclusion remains unproved",
    }
    temp=OUT.with_suffix(OUT.suffix+".tmp")
    temp.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    temp.replace(OUT)
    print(json.dumps({"bands":{str(k):counts[str(k)] for k in range(FIRST,LAST+1)},
        "newly_excluded":cumulative_new,"cumulative_excluded":total,
        "finite_grid_percent":result["cumulative_percent_of_finite_profile_grid"],
        "survivors":len(survivors),"gpu":device.name,"gpu_cpu_parity":True,
        "receipt":OUT.name}),flush=True)


if __name__ == "__main__":
    main()
