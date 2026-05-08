function create_default_smith_simulation_use_cases()::Vector{SimulationUseCase}
    # Time convention: campaign start_time = 0.0 seconds
    # Simulation may begin before campaign start (negative current_time)

    return SimulationUseCase[
        # 1. no_supply
        SimulationUseCase(
            "no_supply",
            CampaignMetadata(86400 * MICRODOLLAR, 0.0, 3600.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(-240.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                0, 0, 0.0, default_eligible_impression_generator,
                0.45, 0.03, default_win_rate_generator,
                60.0, 0.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.0, Int64(0), 0.0)
        ),

        # 2. exactly_1cent_per_second_spend_no_win_delay
        SimulationUseCase(
            "exactly_1cent_per_second_spend_no_win_delay",
            CampaignMetadata(36 * MICRODOLLAR, 0.0, 3600.0, MICROCENT),
            CampaignState(-60.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                1, 0, 0.0, default_eligible_impression_generator,
                1.0, 0.0, default_win_rate_generator,
                1.0, 0.0, default_tick_interval_generator,
                0, true,
                nothing, nothing
            ),
            SimulationExpects(false, 0.05, Int64(3600), 0.051)
        ),

        # 3. exactly_1cent_per_second_spend
        SimulationUseCase(
            "exactly_1cent_per_second_spend",
            CampaignMetadata(36 * MICRODOLLAR, 0.0, 3600.0, MICROCENT),
            CampaignState(-60.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                1, 0, 0.0, default_eligible_impression_generator,
                1.0, 0.0, default_win_rate_generator,
                1.0, 0.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(false, 0.05, Int64(3600), 0.05)
        ),

        # 4. pace_new_campaign_86400_budget
        SimulationUseCase(
            "pace_new_campaign_86400_budget",
            CampaignMetadata(86400 * MICRODOLLAR, 0.0, 3600.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(-240.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                9000, 100, 0.0, default_eligible_impression_generator,
                0.25, 0.03, default_win_rate_generator,
                60.0, 0.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.05, Int64(3_456_000), 0.05)
        ),

        # 5. pace_with_fuzzy_subsecond_intervals
        SimulationUseCase(
            "pace_with_fuzzy_subsecond_intervals",
            CampaignMetadata(86400 * MICRODOLLAR, 0.0, 3600.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(-240.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                9000, 100, 0.0, default_eligible_impression_generator,
                0.25, 0.03, default_win_rate_generator,
                0.5, 0.5, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.05, Int64(3_456_000), 0.05)
        ),

        # 6. pace_with_fuzzy_15s_tick
        SimulationUseCase(
            "pace_with_fuzzy_15s_tick",
            CampaignMetadata(86400 * MICRODOLLAR, 0.0, 3600.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(-240.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                9000, 100, 0.0, default_eligible_impression_generator,
                0.25, 0.03, default_win_rate_generator,
                15.0, 1.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.05, Int64(3_456_000), 0.05)
        ),

        # 7. pace_with_fuzzy_multi_minute_intervals
        SimulationUseCase(
            "pace_with_fuzzy_multi_minute_intervals",
            CampaignMetadata(86400 * MICRODOLLAR, 0.0, 3600.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(-240.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                9000, 100, 0.0, default_eligible_impression_generator,
                0.25, 0.03, default_win_rate_generator,
                15.0, 0.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.05, Int64(3_456_000), 0.05)
        ),

        # 8. hard_catch_up (25-hour campaign, simulation starts near end)
        SimulationUseCase(
            "hard_catch_up",
            CampaignMetadata(86400 * MICRODOLLAR, -86400.0, 3600.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(-240.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                9000, 100, 0.0, default_eligible_impression_generator,
                0.25, 0.03, default_win_rate_generator,
                15.0, 0.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.05, Int64(3_456_000), 0.05)
        ),

        # 9. resumed_near_end_of_flight (pre-loaded 90% spent)
        SimulationUseCase(
            "resumed_near_end_of_flight",
            CampaignMetadata(86400 * MICRODOLLAR, -86400.0, 3600.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(
                -240.0, 0, 0, 0.0, 0.0,
                round(Int64, 3_456_000 * 0.9),
                round(Int64, 3_456_000 * 0.9 * 25 * MICRODOLLAR / MICROMILLI),
                PendingWinQueue(), PendingRateQueue()
            ),
            SimulationConfig(
                9000, 100, 0.0, default_eligible_impression_generator,
                0.25, 0.03, default_win_rate_generator,
                15.0, 0.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.05, Int64(3_456_000), 0.05)
        ),

        # 10. resumed_over_paced_near_end_of_flight (pre-loaded 95% spent)
        SimulationUseCase(
            "resumed_over_paced_near_end_of_flight",
            CampaignMetadata(86400 * MICRODOLLAR, -86400.0, 7200.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(
                -240.0, 0, 0, 0.0, 0.0,
                round(Int64, 3_456_000 * 0.95),
                round(Int64, 3_456_000 * 0.95 * 25 * MICRODOLLAR / MICROMILLI),
                PendingWinQueue(), PendingRateQueue()
            ),
            SimulationConfig(
                9000, 100, 0.0, default_eligible_impression_generator,
                0.25, 0.03, default_win_rate_generator,
                15.0, 0.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.05, Int64(3_456_000), 0.05)
        ),

        # 11. compensate_for_spikes (3-day campaign)
        SimulationUseCase(
            "compensate_for_spikes",
            CampaignMetadata(86400 * MICRODOLLAR, 0.0, 259200.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(-240.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                9000, 100, 0.002, default_eligible_impression_generator,
                0.25, 0.03, default_win_rate_generator,
                15.0, 0.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.05, Int64(3_456_000), 0.05)
        ),

        # 12. unstable_eligible_impressions_with_spikes
        SimulationUseCase(
            "unstable_eligible_impressions_with_spikes",
            CampaignMetadata(86400 * MICRODOLLAR, 0.0, 86400.0, div(25 * MICRODOLLAR, MICROMILLI)),
            CampaignState(-240.0, 0, 0, 0.0, 0.0, 0, 0, PendingWinQueue(), PendingRateQueue()),
            SimulationConfig(
                10, 200, 0.004, spike_and_trickle_generator,
                0.80, 0.03, default_win_rate_generator,
                15.0, 0.0, default_tick_interval_generator,
                0, true,
                impression_event_latency_generator,
                default_eligible_impression_latency_generator
            ),
            SimulationExpects(true, 0.0, Int64(3_456_000), 0.50)
        ),
    ]
end
