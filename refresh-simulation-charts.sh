#!/usr/bin/env bash
set -euo pipefail

julia --project="$(dirname "$0")" -e '
using DDEExamples
run_pid_simulation(; verbose=false)
run_smith_simulation(; verbose=false)
run_corrdenom_simulation(; verbose=false)
run_mpc_simulation(; verbose=false)
'
