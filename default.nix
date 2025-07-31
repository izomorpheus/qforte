{ pkgs ? import <nixpkgs> {} }:
let
  fhs = pkgs.buildFHSEnv {
    name = "qforte";
    targetPkgs = pkgs: with pkgs; [
      micromamba
      which
      vim
      gcc
      gdb
      glibc
      starship
      gnumake
    ];

    runScript = "bash --login";

    profile = ''
      set -e
      export SHELL=/usr/bin/bash
      export MAMBA_ROOT_PREFIX="$PWD/.mamba"
      if [ "$0" = bash ]; then
        eval "$(micromamba shell hook --shell=bash)"
        eval "$(starship init bash)"
        micromamba activate qforte-default-env
      else
        eval "$(micromamba shell hook --shell=posix)"
      fi
      if [ ! -d "$MAMBA_ROOT_PREFIX/envs/qforte-default-env" ]; then
        micromamba create -n qforte-default-env -y
        micromamba install -n qforte-default-env -y -c conda-forge python=3.8 openblas psi4 cmake pytest
      fi
      echo DONE: source /etc/profile
      set +e
    '';
  };

in fhs.env
