function plot_simulation(name::String, data::Vector{PlotDataPoint},
                        start_time::Float64, end_time::Float64;
                        output_dir::String="src/pid/plots")
    isempty(data) && return nothing

    mkpath(output_dir)

    times = [dp.time for dp in data]
    tick_deltas = [dp.tick_delta for dp in data]
    on_target = [dp.linear_on_target_ratio for dp in data]
    expected_util = [dp.expected_budget_utilization for dp in data]
    budget_util = [dp.budget_utilization for dp in data]
    throttles = [dp.throttle for dp in data]
    eligible_rates = [dp.eligible_impression_rate for dp in data]

    # Compute -log10(eligible_rate) for display (negative axis)
    log_rates = [r > 0 ? -log10(Float64(r)) : NaN for r in eligible_rates]

    # Compute x-axis limits
    t_min = times[1]
    t_max = max(times[end], end_time)
    margin = (t_max - t_min) * 0.02
    xl = (t_min - margin, t_max + margin)

    # Trim data to visible range for plotting
    visible_mask = [t_min - margin .<= t .<= t_max + margin for t in times]
    vis_times = times[visible_mask]
    vis_tick_deltas = tick_deltas[visible_mask]
    vis_log_rates = log_rates[visible_mask]
    vis_on_target = on_target[visible_mask]
    vis_expected_util = expected_util[visible_mask]
    vis_budget_util = budget_util[visible_mask]
    vis_throttles = throttles[visible_mask]

    # Left Y-axis: tick delta and -log10(rate)
    p = Plots.plot(vis_times, vis_tick_deltas;
        label="Tick Delta (s)", color=:lightblue, linewidth=1, alpha=0.6,
        xlabel="Time (seconds)", ylabel="Tick Δ (s) / -log₁₀(rate)",
        title=replace(name, "_" => " "),
        legend=:topright, size=(1200, 600), xlims=xl)

    Plots.plot!(p, vis_times, vis_log_rates;
        label="-log₁₀(Eligible Rate)", color=:cyan, linewidth=1, alpha=0.7)

    # Right Y-axis: ratios and throttle (use twinx)
    p2 = Plots.twinx(p)
    Plots.plot!(p2, vis_times, vis_on_target;
        label="On-Target Ratio", color=:red, linewidth=2, ylabel="Ratio", xlims=xl)
    Plots.plot!(p2, vis_times, vis_expected_util;
        label="Expected Utilization", color=:orange, linewidth=2)
    Plots.plot!(p2, vis_times, vis_budget_util;
        label="Budget Utilization", color=:green, linewidth=2)
    Plots.plot!(p2, vis_times, vis_throttles;
        label="Throttle", color=:purple, linewidth=2)

    # Reference lines
    Plots.hline!(p2, [1.0]; label="", linestyle=:dash, color=:gray, linewidth=1)

    # Campaign window markers on p2 (twinx has correct x coordinate mapping)
    if start_time >= t_min - margin
        Plots.vline!(p2, [start_time]; label="Campaign Start", linestyle=:dash, color=:goldenrod, linewidth=1.5)
    end
    Plots.vline!(p2, [end_time]; label="Campaign End", linestyle=:dash, color=:red, linewidth=1.5)

    # Force xlims on both axes after all data is plotted
    Plots.xlims!(p, xl)
    Plots.xlims!(p2, xl)

    safe_name = replace(name, "/" => "_", " " => "_")
    filepath = joinpath(output_dir, "$(safe_name).png")
    Plots.savefig(p, filepath)
    println("Plot saved to $filepath")
    return filepath
end
