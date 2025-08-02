{ pkgs ? import <nixpkgs> {} }:
let
  fhs = pkgs.buildFHSEnv {
    name = "qforte";
    targetPkgs = pkgs: with pkgs; [
      micromamba
      gnumake
      gcc
    ];

    runScript = "bash --login";

    profile = ''
      set -e

      export MAMBA_ROOT_PREFIX="$PWD/.mamba"

      if [ "$0" = zsh ]; then
        eval "$(micromamba shell hook --shell=zsh)"
      elif [ "$0" = bash ]; then
        eval "$(micromamba shell hook --shell=bash)"
      else
        eval "$(micromamba shell hook --shell=posix)"
      fi

      if [ ! -d "$MAMBA_ROOT_PREFIX/envs/qforte-qiskit-env" ]; then
        micromamba create -n qforte-default-env -y -c conda-forge python=3.10 openblas psi4 cmake pytest \
        qiskit qiskit-ibm-runtime qiskit-aer
      fi

      micromamba activate qforte-default-env

      set +e
    '';
  };

in fhs.env
