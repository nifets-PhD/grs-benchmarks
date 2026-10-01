GRS_PATH := $(HOME)/Code/GeneRegulatorySystems.jl
JP_PATH  := $(HOME)/Code/JumpProcesses.jl
GRS_URL  := https://github.com/nifets/GeneRegulatorySystems.jl
JP_URL   := https://github.com/nifets/JumpProcesses.jl
GRS_REV  := d8a5a364d2e61d836895db1fe7e9e33387796437
JP_REV   := cb8b9a0914eb319c91d74ceb9d6d02a9d01632c7

STAMP := .make

.DEFAULT_GOAL := help
.NOTPARALLEL:

setup: setup-julia setup-r setup-python
smoke: smoke-grs smoke-copasi smoke-dyngen smoke-experiments

setup-julia:  $(STAMP)/julia
setup-r:      $(STAMP)/r
setup-python: $(STAMP)/python

$(STAMP):
	mkdir -p $@

$(STAMP)/julia: Project.toml Manifest.toml | $(STAMP)
	julia --project=. -e 'using Pkg; Pkg.instantiate()'
	@touch $@

$(STAMP)/r: renv.lock | $(STAMP)
	Rscript -e 'renv::restore(prompt = FALSE)'
	@touch $@

$(STAMP)/python: requirements.txt | $(STAMP)
	python3 -m venv .venv
	.venv/bin/pip install -q -r requirements.txt
	@touch $@

smoke-grs: $(STAMP)/julia
	julia --project=. runners/smoke_grs.jl

smoke-copasi: $(STAMP)/python
	.venv/bin/python runners/smoke_copasi.py

smoke-dyngen: $(STAMP)/r
	Rscript runners/smoke_dyngen.R

smoke-experiments: setup
	@for e in experiments/*.sh; do echo "=== $$e"; SMOKE=1 $$e || exit 1; done

experiments:
	$(MAKE) -k kronecker-scaling kronecker-accuracy dyngen-model promoter-approximation

kronecker-scaling: setup
	experiments/$@.sh

kronecker-accuracy: setup
	experiments/$@.sh

dyngen-model: setup
	experiments/$@.sh

promoter-approximation: setup
	experiments/$@.sh

figures: figures/cross-engine-kronecker.pdf figures/scaling-kronecker-full.pdf figures/dyngen-model.tex \
    figures/promoter-approximation.pdf

figures/cross-engine-kronecker.pdf: figures/cross-engine-kronecker.jl \
    measurements/kronecker-scaling.csv measurements/kronecker-accuracy.csv \
    measurements/kronecker-accuracy-cost.csv
	julia --project=. $<

figures/scaling-kronecker-full.pdf: figures/scaling-kronecker-full.jl \
    measurements/kronecker-scaling.csv
	julia --project=. $<

figures/dyngen-model.tex: figures/dyngen-model.jl measurements/dyngen-model.csv
	julia --project=. $<

figures/promoter-approximation.pdf: figures/promoter-approximation.jl \
    measurements/promoter-approximation/costs.csv
	julia --project=. $<

dev:
	julia --project=. -e 'using Pkg; Pkg.develop([PackageSpec(path="$(GRS_PATH)"), PackageSpec(path="$(JP_PATH)")])'

pin:
	julia --project=. -e 'using Pkg; Pkg.add([PackageSpec(url="$(GRS_URL)", rev="$(GRS_REV)"), PackageSpec(url="$(JP_URL)", rev="$(JP_REV)")])'

bootstrap: dev
	Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")'
	Rscript -e 'renv::init(bare = TRUE)'
	Rscript -e 'renv::install(c("GillespieSSA2", "dyngen")); renv::snapshot(type = "all")'
	python3 -m venv .venv
	.venv/bin/pip install -q copasi-basico
	.venv/bin/pip freeze > requirements.txt

clean:
	rm -rf $(STAMP)

help:
	@grep -hE '^[a-z][a-z0-9-]*:' Makefile | cut -d: -f1 | sort -u

.PHONY: setup setup-julia setup-r setup-python smoke smoke-grs smoke-copasi \
        smoke-dyngen smoke-experiments experiments kronecker-scaling kronecker-accuracy dyngen-model promoter-approximation figures dev pin \
        bootstrap clean help
