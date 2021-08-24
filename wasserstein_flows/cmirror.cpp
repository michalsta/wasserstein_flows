#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>


class CMirror
{
public:
    CMirror(const std::vector<double>&& _masses, const std::vector<double>&& _probs, const std::vector<int32_t>&& _idxes, double _exp_ab_cost, double _the_ab_cost);

    const std::vector<double> masses;
    const std::vector<double> probs;
    const std::vector<int32_t> idxes;

    const size_t nconfs;

    const double exp_ab_cost;
    const double the_ab_cost;
    const double tot_ab_cost;

    class Graph{
    public:
        std::vector<double> inline_flows;
        std::vector<double> costs;
        std::vector<double> directed_probs;
        std::vector<double> from_abyss_flows;

        void send(size_t src, size_t dst, double howmuch);
    } G;

    double yankable_flow(size_t idx) const;
    double stuffable_flow(size_t idx) const;

    double get_abyss_cost(size_t idx) const;

    std::pair<double, double> stuffing_into_abyss_cost(size_t idx) const;
    std::pair<double, double> yanking_out_of_abyss_cost(size_t idx) const;


    std::pair<double, double> delta_cost(size_t src, size_t tgt);

    inline bool pushout_once(size_t idx);
};


CMirror::CMirror(const std::vector<double>&& _masses, const std::vector<double>&& _probs, const std::vector<int32_t>&& _idxes, double _exp_ab_cost, double _the_ab_cost) :
    masses(std::move(_masses)),
    probs(std::move(_probs)),
    idxes(std::move(_idxes)),
    nconfs(idxes.size()),
    exp_ab_cost(_exp_ab_cost),
    the_ab_cost(_the_ab_cost),
    tot_ab_cost(_exp_ab_cost + _the_ab_cost)
{
    G.inline_flows.resize(nconfs-1);

    G.costs.reserve(nconfs-1);
    for(size_t ii=0; ii < nconfs-1; ii++)
        G.costs.push_back(masses[ii+1] - masses[ii]);

    G.directed_probs.reserve(nconfs);
    G.from_abyss_flows.reserve(nconfs);
    for(size_t ii=0; ii < nconfs; ii++)
    {
        const double dp = idxes[ii] >= 0 ? probs[ii] : -probs[ii];
        G.directed_probs.push_back(dp);
        G.from_abyss_flows.push_back(-dp);
    }
}


inline double CMirror::yankable_flow(size_t idx) const
{
    if(G.directed_probs[idx] > 0.0)
        return -G.from_abyss_flows[idx];
    return -G.directed_probs[idx] - G.from_abyss_flows[idx];
}

inline double CMirror::stuffable_flow(size_t idx) const
{
    if(G.directed_probs[idx] > 0.0)
        return G.directed_probs[idx] + G.from_abyss_flows[idx];
    return G.from_abyss_flows[idx];
}

inline double CMirror::get_abyss_cost(size_t idx) const
{
    return idxes[idx] == -1 ? exp_ab_cost : the_ab_cost;
}


inline std::pair<double, double> CMirror::stuffing_into_abyss_cost(size_t idx) const
{
    if(G.from_abyss_flows[idx] > 0.0)
        return {-get_abyss_cost(idx), G.from_abyss_flows[idx]};
    return {get_abyss_cost(idx), probs[idx] + G.from_abyss_flows[idx]};
}


inline std::pair<double, double> CMirror::yanking_out_of_abyss_cost(size_t idx) const
{
    if(G.from_abyss_flows[idx] < 0.0)
        return {-get_abyss_cost(idx), -G.from_abyss_flows[idx]};
    return {get_abyss_cost(idx), probs[idx] - G.from_abyss_flows[idx]};
}


std::pair<double, double> CMirror::delta_cost(size_t src, size_t tgt)
{
    if(src == tgt)
        return {0.0, 0.0};

    auto [cost, flow] = yanking_out_of_abyss_cost(src);
    auto [c, f] = stuffing_into_abyss_cost(tgt);
    cost += c;
    flow = (std::min)(f, flow);

    if(flow == 0.0)
        return {0.0, 0.0};

    if(src < tgt)
    {
        for(size_t ii = src; ii < tgt; ii++)
            if(G.inline_flows[ii] >= 0.0)
                cost += G.costs[ii];
            else
            {
                cost -= G.costs[ii];
                flow = (std::min)(flow, -G.inline_flows[ii]);
            }
    }
    else
    {
        for(size_t ii = tgt; ii < src; ii++)
            if(G.inline_flows[ii] <= 0.0)
                cost += G.costs[ii];
            else
            {
                cost -= G.costs[ii];
                flow = (std::min)(flow, G.inline_flows[ii]);
            }
    }

    return {cost, flow};
}

inline void CMirror::Graph::send(size_t src, size_t tgt, double howmuch)
{
    if(src == tgt)
        return;
    if(tgt < src)
    {
        std::swap(src, tgt);
        howmuch = -howmuch;
    }

    from_abyss_flows[src] += howmuch;
    from_abyss_flows[tgt] -= howmuch;

    for(size_t ii = src; ii < tgt; ii++)
        inline_flows[ii] += howmuch;
}

/*
std::tuple<size_t, double, double> CMirror::closest_left_improving(size_t src)
{
    auto [cost, flow] = yanking_out_of_abyss_cost(src);
    if(flow == 0.0)
        return {0, 0.0, 0.0};

    for(size_t ii = src-1; ii != std::numeric_limits<decltype(ii)>::max(); ii--)
    {
        if(G.inline_flows[ii] <= 0.0)
            cost += G.costs[ii];
        else
        {
            cost -= G.costs[ii];
            flow = (std::min)(flow, -G.inline_flows[ii]);
        }
        auto [c, f] = stuffing_into_abyss_cost(ii);
        if(f > 0.0 and cost + c < 0.0)
            return {ii, cost+c, (std::min)(flow, f)};
    }
    return {0, 0.0, 0.0};
}*/


inline bool CMirror::pushout_once(size_t idx)
{
        if(yankable_flow(idx) == 0.0)
            return false;
        double dcost_left = 0.0;
        double fl_left = 0.0;
        ssize_t left_idx = idx-1;
        while(left_idx >= 0)
        {
            if(stuffable_flow(left_idx) > 0.0)
            {
                std::tie(dcost_left, fl_left) = delta_cost(idx, left_idx);
                if((dcost_left < 0.0 and fl_left > 0.0) or dcost_left > tot_ab_cost)
                    break;
            }
            left_idx -= 1;
        }
        double dcost_right = 0.0;
        double fl_right = 0.0;
        size_t right_idx = idx+1;
        while(right_idx < nconfs)
        {
            if(stuffable_flow(right_idx) > 0.0)
            {
                std::tie(dcost_right, fl_right) = delta_cost(idx, right_idx);
                if((dcost_right < 0.0 and fl_right > 0.0) or dcost_right > tot_ab_cost)
                    break;
            }
            right_idx += 1;
        }

        if(dcost_left >= 0.0 and dcost_right >= 0.0)
            return false;

        if(dcost_left < dcost_right)
            G.send(idx, left_idx, fl_left);
        else
            G.send(idx, right_idx, fl_right);
        return true;
}
