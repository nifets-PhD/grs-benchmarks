import argparse
import os
import resource
import sys
import time

import basico
import pandas


METHODS = ["stochastic", "directMethod", "tauLeap", "adaptiveSA",
           "hybridlsoda", "hybrid", "hybridode45"]

EPSILON_METHODS = {"tauLeap", "adaptiveSA"}

p = argparse.ArgumentParser()
p.add_argument("out")
p.add_argument("sbml")
p.add_argument("--seed", type=int, default=1)
p.add_argument("--method", default="stochastic", choices=METHODS)
p.add_argument("--duration", type=float, default=20000.0)
p.add_argument("--samples", type=int, default=10)
p.add_argument("--trajectories", type=int, default=3)
p.add_argument("--epsilon", type=float, default=None)
p.add_argument("--max-steps", dest="max_steps", type=int, default=10**9)
a = p.parse_args()

if a.method in EPSILON_METHODS and a.epsilon is None:
    p.error(f"--epsilon is required for {a.method}; COPASI's silent default is not a "
            f"configuration we can name or reproduce")
if a.method not in EPSILON_METHODS and a.epsilon is not None:
    p.error(f"{a.method} has no Epsilon parameter")


def peak_rss_bytes():
    n = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return n if sys.platform == "darwin" else n * 1024


def report(msg, t0=time.time()):
    print(f"[{time.time() - t0:.0f}s] {msg}", flush=True)


report("loading sbml")
t_load = time.time()
model = basico.load_model(a.sbml)
t_load = time.time() - t_load

def epsilon():
    return model.getTask("Time-Course").getMethod().getParameter("Epsilon")


def course(seed):
    basico.set_task_settings(
        basico.T.TIME_COURSE,
        {"method": {"name": a.method, "Use Random Seed": True,
                    "Random Seed": seed, "Max Internal Steps": a.max_steps}},
        model=model,
    )
    if a.epsilon is not None:
        epsilon().setDblValue(a.epsilon)
    df = basico.run_time_course(
        duration=a.duration, automatic=False, intervals=a.samples,
        method=a.method, model=model, max_steps=a.max_steps,
    )
    if len(df) != a.samples + 1:
        raise RuntimeError(
            f"TRUNCATED: seed {seed} returned {len(df)} of {a.samples + 1} samples "
            f"(last t={df.index[-1]:.1f} of {a.duration:.1f}); "
            f"method={a.method} model={a.sbml}"
        )
    return df

if a.epsilon is not None:
    basico.run_time_course(duration=0.0, automatic=False, intervals=1,
                           method=a.method, model=model, max_steps=a.max_steps)
    epsilon().setDblValue(a.epsilon)
    got = epsilon().getValue()
    if abs(got - a.epsilon) > 1e-12:
        raise RuntimeError(f"Epsilon not applied: wanted {a.epsilon}, method has {got}")

report("simulating")
frames = []
t_sim = time.time()
for r in range(a.trajectories):
    frames.append(course(a.seed + r))
    report(f"trajectory {r + 1}/{a.trajectories}")
t_sim = time.time() - t_sim
per_traj = t_sim / a.trajectories

report("writing")
os.makedirs(a.out, exist_ok=True)

rows = []
for r, df in enumerate(frames):
    df = df.reset_index().rename(columns={"Time": "t"})
    long = df.melt(id_vars="t", var_name="name", value_name="value")
    long["name"] = long["name"].str.replace(r'^var"(.*)"\(t\)$', r"\1", regex=True)
    long["name"] = long["name"].str.replace("\u208a", ".", regex=False)
    long.insert(0, "path", f"{a.seed + r}")
    long["sample"] = (long["t"] / (a.duration / a.samples)).round().astype("int32")
    rows.append(long)
counts = pandas.concat(rows, ignore_index=True)

import pyarrow
import pyarrow.feather

pyarrow.feather.write_feather(
    pyarrow.Table.from_pandas(counts[["path", "sample", "t", "name", "value"]]),
    os.path.join(a.out, "counts.arrow"),
    compression="uncompressed",
)

with open(os.path.join(a.out, "config.txt"), "w") as io:
    io.write("engine=copasi\n")
    io.write(f"model={a.sbml}\n")
    io.write(f"copasi_max_steps={a.max_steps}\n")
    for key in ("method", "seed", "duration", "samples", "trajectories", "epsilon"):
        io.write(f"{key}={getattr(a, key)}\n")

with open(os.path.join(a.out, "metrics.csv"), "w") as io:
    io.write("metric,value\n")
    io.write(f"build_seconds,{t_load}\n")
    io.write(f"simulate_seconds,{t_sim}\n")
    io.write(f"simulate_seconds_per_trajectory,{per_traj}\n")
    io.write(f"rss_peak_bytes,{peak_rss_bytes()}\n")

print(
    f"copasi {a.out}: {t_sim:.1f} s simulate "
    f"({a.trajectories} traj, {per_traj:.1f} s/traj), "
    f"{t_load:.1f} s load, {peak_rss_bytes() / 2**30:.2f} GiB peak rss, {len(counts)} count rows",
    flush=True,
)
