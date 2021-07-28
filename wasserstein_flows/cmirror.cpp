#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>


class CMirror
{
public:
    CMirror(const std::vector<double>&& _masses, const std::vector<double>&& _probs, const std::vector<int32_t>&& _idxes, double _exp_ab_cost, double _the_ab_cost);

    std::vector<double> masses;
    std::vector<double> probs;
    std::vector<int32_t> idxes;

    const size_t nconfs;

    const double exp_ab_cost;
    const double the_ab_cost;

    std::vector<double> inline_flows;
    std::vector<double> costs;
    std::vector<double> directed_probs;
    std::vector<double> from_abyss_flows;

    double yankable_flow(size_t idx) const;
    double stuffable_flow(size_t idx) const;

    double get_abyss_cost(size_t idx) const;

    std::pair<double, double> stuffing_into_abyss_cost(size_t idx) const;
    std::pair<double, double> yanking_out_of_abyss_cost(size_t idx) const;


    std::pair<double, double> delta_cost(size_t src, size_t tgt);
};


CMirror::CMirror(const std::vector<double>&& _masses, const std::vector<double>&& _probs, const std::vector<int32_t>&& _idxes, double _exp_ab_cost, double _the_ab_cost) :
    masses(std::move(_masses)),
    probs(std::move(_probs)),
    idxes(std::move(_idxes)),
    nconfs(idxes.size()),
    exp_ab_cost(_exp_ab_cost),
    the_ab_cost(_the_ab_cost)
{
    inline_flows.resize(nconfs-1);

    costs.reserve(nconfs-1);
    for(size_t ii=0; ii < nconfs-1; ii++)
        costs.push_back(masses[ii+1] - masses[ii]);

    directed_probs.reserve(nconfs);
    from_abyss_flows.reserve(nconfs);
    for(size_t ii=0; ii < nconfs; ii++)
    {
        const double dp = idxes[ii] >= 0 ? probs[ii] : -probs[ii];
        directed_probs.push_back(dp);
        from_abyss_flows.push_back(-dp);
    }
}


inline double CMirror::yankable_flow(size_t idx) const
{
    if(directed_probs[idx] > 0.0)
        return -from_abyss_flows[idx];
    return -directed_probs[idx] - from_abyss_flows[idx];
}

inline double CMirror::stuffable_flow(size_t idx) const
{
    if(directed_probs[idx] > 0.0)
        return directed_probs[idx] + from_abyss_flows[idx];
    return from_abyss_flows[idx];
}

inline double CMirror::get_abyss_cost(size_t idx) const
{
    return idxes[idx] == -1 ? exp_ab_cost : the_ab_cost;
}


inline std::pair<double, double> CMirror::stuffing_into_abyss_cost(size_t idx) const
{
    if(from_abyss_flows[idx] > 0.0)
        return {-get_abyss_cost(idx), from_abyss_flows[idx]};
    return {get_abyss_cost(idx), probs[idx] + from_abyss_flows[idx]};
}


inline std::pair<double, double> CMirror::yanking_out_of_abyss_cost(size_t idx) const
{
    if(from_abyss_flows[idx] < 0.0)
        return {-get_abyss_cost(idx), -from_abyss_flows[idx]};
    return {get_abyss_cost(idx), probs[idx] - from_abyss_flows[idx]};
}


std::pair<double, double> CMirror::delta_cost(size_t src, size_t tgt)
{
    if(src == tgt)
        return {0.0, 0.0};

    auto [cost, flow] = yanking_out_of_abyss_cost(src);
    auto [c, f] = stuffing_into_abyss_cost(tgt);
    cost += c;
    flow = (std::min)(f, flow);

    if(src < tgt)
    {
        for(size_t ii = src; ii < tgt; ii++)
            if(inline_flows[ii] >= 0.0)
                cost += costs[ii];
            else
            {
                cost -= costs[ii];
                flow = (std::min)(flow, -inline_flows[ii]);
            }
    }
    else
    {
        for(size_t ii = tgt; ii < src; ii++)
            if(inline_flows[ii] <= 0.0)
                cost += costs[ii];
            else
            {
                cost -= costs[ii];
                flow = (std::min)(flow, inline_flows[ii]);
            }
    }

    return {cost, flow};
}
