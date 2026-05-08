using Test
using DDEExamples

@testset "PID Discrete Simulation" begin
    use_cases = create_default_pid_pacer_simulation_use_cases()

    for uc in use_cases
        @testset "$(uc.name)" begin
            hard_v, soft_v, iters = run_single_scenario!(uc; plot=true, verbose=false)

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
