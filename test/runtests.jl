using Test
using DDEExamples

@testset "PID Discrete Simulation" begin
    use_cases = create_default_pid_pacer_simulation_use_cases()

    for (i, uc) in enumerate(use_cases)
        @testset "$(uc.name)" begin
            hard_v, soft_v, iters = run_single_scenario!(uc; plot=true, verbose=false, index=i-1)

            if !uc.expects.ignore_pacing_violations_while_in_flight
                @test hard_v == 0
                @test soft_v == 0
            end

            if !isnothing(uc.expects.expected_impressions)
                expected = uc.expects.expected_impressions
                if expected > 0
                    delta = round(Int64, expected * uc.expects.expected_impressions_delta_percentage)
                    @test abs(uc.state.spent_impressions - expected) <= delta
                end
            end
        end
    end
end

@testset "Smith Discrete Simulation" begin
    use_cases = create_default_smith_simulation_use_cases()

    for (i, uc) in enumerate(use_cases)
        @testset "$(uc.name)" begin
            hard_v, soft_v, iters = run_single_smith_scenario!(uc; plot=true, verbose=false, index=i-1)

            if !uc.expects.ignore_pacing_violations_while_in_flight
                @test hard_v == 0
                @test soft_v == 0
            end

            if !isnothing(uc.expects.expected_impressions)
                expected = uc.expects.expected_impressions
                if expected > 0
                    delta = round(Int64, expected * uc.expects.expected_impressions_delta_percentage)
                    @test abs(uc.state.spent_impressions - expected) <= delta
                end
            end
        end
    end
end

@testset "CorrDenom Discrete Simulation" begin
    use_cases = create_default_corrdenom_simulation_use_cases()

    for (i, uc) in enumerate(use_cases)
        @testset "$(uc.name)" begin
            hard_v, soft_v, iters = run_single_corrdenom_scenario!(uc; plot=true, verbose=false, index=i-1)

            if !uc.expects.ignore_pacing_violations_while_in_flight
                @test hard_v == 0
                @test soft_v == 0
            end

            if !isnothing(uc.expects.expected_impressions)
                expected = uc.expects.expected_impressions
                if expected > 0
                    delta = round(Int64, expected * uc.expects.expected_impressions_delta_percentage)
                    @test abs(uc.state.spent_impressions - expected) <= delta
                end
            end
        end
    end
end

@testset "MPC Discrete Simulation" begin
    use_cases = create_default_mpc_simulation_use_cases()

    for (i, uc) in enumerate(use_cases)
        @testset "$(uc.name)" begin
            hard_v, soft_v, iters = run_single_mpc_scenario!(uc; plot=true, verbose=false, index=i-1)

            if !uc.expects.ignore_pacing_violations_while_in_flight
                @test hard_v == 0
                @test soft_v == 0
            end

            if !isnothing(uc.expects.expected_impressions)
                expected = uc.expects.expected_impressions
                if expected > 0
                    delta = round(Int64, expected * uc.expects.expected_impressions_delta_percentage)
                    @test abs(uc.state.spent_impressions - expected) <= delta
                end
            end
        end
    end
end

@testset "Artstein Discrete Simulation" begin
    use_cases = create_default_artstein_simulation_use_cases()

    for (i, uc) in enumerate(use_cases)
        @testset "$(uc.name)" begin
            hard_v, soft_v, iters = run_single_artstein_scenario!(uc; plot=true, verbose=false, index=i-1)

            if !uc.expects.ignore_pacing_violations_while_in_flight
                @test hard_v == 0
                @test soft_v == 0
            end

            if !isnothing(uc.expects.expected_impressions)
                expected = uc.expects.expected_impressions
                if expected > 0
                    delta = round(Int64, expected * uc.expects.expected_impressions_delta_percentage)
                    @test abs(uc.state.spent_impressions - expected) <= delta
                end
            end
        end
    end
end
