# flake-qiskit.nix
{
  description = "qforte (Qiskit): dev shell + build-and-exit";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:

    let
      pkgs = import nixpkgs { inherit system; };

      projectRoot = toString ../..;

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

          if [ ! -d "$MAMBA_ROOT_PREFIX/envs/qforte-default-env" ]; then
            micromamba create -n qforte-default-env -y -c conda-forge python=3.8 openblas psi4 cmake pytest
          fi

          micromamba activate qforte-default-env

          set +e
        '';
      };

    in {
      ## 1) Interactive dev shell
      devShells.default = pkgs.mkShell {
        shellHook  = ''
          starting qforte dev shell...
          exec ${fhs.out}/bin/qforte
        '';
      };

      devShells.build = pkgs.mkShell {
        shellHook  = ''
          echo building and installing qforte...
          cd ${projectRoot} || exit 1
          exec ${fhs.out}/bin/qforte ./scripts/sh/build-mamba.sh
        '';
      };

    });
}

