FROM mambaorg/micromamba:1.5.10

COPY envs/qc-secundario.yml /tmp/qc-secundario.yml

RUN micromamba install -y -n base -f /tmp/qc-secundario.yml && \
    micromamba clean --all --yes

WORKDIR /workspace
