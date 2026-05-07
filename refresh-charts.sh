#!/usr/bin/env bash
set -euo pipefail

julia --project="$(dirname "$0")" -e '
using DDEExamples
demo()
demo_zero_delay()
demo_budget_delay()
demo_budget_delay_with_noise()
demo_budget_controllers()
demo_budget_controllers_noise()
demo_pid_pacer()
demo_pid_pacer_noise()
demo_demand_spike()
demo_demand_two_spikes()
demo_demand_spike_then_drop()
demo_demand_random()
'
