struct SpendEvent
    impressions::Int64
    budget::Int64
    triggered_at::Float64
    delivered_at::Float64
end

struct RateObservation
    rate::Int64
    observed_at::Float64
    delivered_at::Float64
end

mutable struct PendingWinQueue
    events::Vector{SpendEvent}
end

PendingWinQueue() = PendingWinQueue(SpendEvent[])

mutable struct PendingRateQueue
    events::Vector{RateObservation}
end

PendingRateQueue() = PendingRateQueue(RateObservation[])

# ── PendingWinQueue Operations ──────────────────────────────────────────────

function enqueue!(q::PendingWinQueue, event::SpendEvent)
    idx = searchsortedfirst(q.events, event; by=e -> e.delivered_at)
    insert!(q.events, idx, event)
end

function dequeue_until!(q::PendingWinQueue, now::Float64)::Vector{SpendEvent}
    split_idx = 0
    for i in eachindex(q.events)
        if q.events[i].delivered_at <= now
            split_idx = i
        else
            break
        end
    end
    split_idx == 0 && return SpendEvent[]
    arrived = q.events[1:split_idx]
    deleteat!(q.events, 1:split_idx)
    return arrived
end

function total_budget(q::PendingWinQueue)::Int64
    isempty(q.events) && return Int64(0)
    sum(e.budget for e in q.events)
end

function total_impressions(q::PendingWinQueue)::Int64
    isempty(q.events) && return Int64(0)
    sum(e.impressions for e in q.events)
end

Base.length(q::PendingWinQueue) = length(q.events)

# ── PendingRateQueue Operations ─────────────────────────────────────────────

function enqueue!(q::PendingRateQueue, obs::RateObservation)
    idx = searchsortedfirst(q.events, obs; by=e -> e.delivered_at)
    insert!(q.events, idx, obs)
end

function deliver_latest!(q::PendingRateQueue, now::Float64)::Tuple{Int64,Bool}
    split_idx = 0
    for i in eachindex(q.events)
        if q.events[i].delivered_at <= now
            split_idx = i
        else
            break
        end
    end
    split_idx == 0 && return (Int64(0), false)
    rate = q.events[split_idx].rate
    deleteat!(q.events, 1:split_idx)
    return (rate, true)
end

Base.length(q::PendingRateQueue) = length(q.events)
