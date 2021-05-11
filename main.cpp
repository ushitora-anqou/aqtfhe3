#include "aqtfhe3.hpp"

#include "benchmark.hpp"
#include "test.hpp"

#include <execution>
#include <queue>
#include <set>
#include <sstream>

template <class P>
typename trlwe_lvl1<P>::poly_torus uint2weight(uint64_t n)
{
    typename trlwe_lvl1<P>::poly_torus w;
    const torus mu = 1u << 31;  // 1/2
    for (size_t i = 0; i < P::N; i++)
        w[i] = ((n >> i) & 1u) ? mu : 0;
    return w;
}

template <class P>
std::string weight2bitstring(const typename trlwe_lvl1<P>::poly_torus &weight)
{
    const torus mu25 = 1u << 30, mu75 = (1u << 30) + (1u << 31);
    std::stringstream ss;
    for (size_t i = 0; i < P::N; i++) {
        if (i % 32 == 0)
            ss << "\n";
        else if (i % 8 == 0)
            ss << " ";
        if (mu25 <= weight[i] && weight[i] <= mu75)
            ss << 1;
        else
            ss << 0;
    }
    return ss.str();
}

template <class P>
uint64_t weight2uint(const typename trlwe_lvl1<P>::poly_torus &weight)
{
    const torus mu25 = 1u << 30, mu75 = (1u << 30) + (1u << 31);
    uint64_t n = 0;
    for (size_t i = 0; i < P::N; i++)
        n |= (((mu25 <= weight[i] && weight[i] <= mu75) ? 1 : 0) << i);
    return n;
}

class Graph {
public:
    using State = int;

private:
    struct TableItem {
        State index, child0, child1;
        uint64_t marker;
    };
    std::vector<TableItem> table_;
    std::vector<std::vector<State>> states_at_depth_;

public:
    Graph()
    {
        table_ = {
#include "test_table.inc"
        };
    }

    size_t size() const
    {
        return table_.size();
    }

    State next_state(State state, bool input) const
    {
        auto &t = table_.at(state);
        return input ? t.child1 : t.child0;
    }

    State initial_state() const
    {
        return 0;
    }

    uint64_t marker_of_state(State state) const
    {
        return table_.at(state).marker;
    }

    std::set<uint64_t> markers() const
    {
        std::set<uint64_t> ret;
        for (auto &&t : table_)
            ret.insert(t.marker);
        return ret;
    }

    void reserve_states_at_depth(size_t depth)
    {
        states_at_depth_.clear();
        states_at_depth_.shrink_to_fit();
        states_at_depth_.reserve(depth);

        std::set<int> sts0, sts1;
        sts0.insert(initial_state());
        for (size_t i = 0; i < depth; i++) {
            states_at_depth_.emplace_back(sts0.begin(), sts0.end());

            for (State st : sts0) {
                sts1.insert(next_state(st, false));
                sts1.insert(next_state(st, true));
            }
            {
                using std::swap;
                swap(sts0, sts1);
            }
        }
    }

    std::vector<State> states_at_depth(size_t depth) const
    {
        return states_at_depth_.at(depth);
    }
};

template <class P>
trlwe_lvl1<P> trlwe_lvl1_zero()
{
    typename trlwe_lvl1<P>::poly_torus m = {};
    return trlwe_lvl1<P>::trivial_encrypt_poly_torus(m);
}

template <class P>
trlwe_lvl1<P> eval_det_wfa(const Graph &gr,
                           const std::vector<trgsw_lvl1_fft<P>> &input)
{
    auto get_w = [&gr](Graph::State s) {
        auto m = gr.marker_of_state(s);
        auto w_s = trlwe_lvl1<P>::trivial_encrypt_poly_torus(uint2weight<P>(m));
        return w_s;
    };
    auto mult_X_1 = [](const trlwe_lvl1<P> &src) {
        trlwe_lvl1<P> out;
        poly_mult_by_X_k(out.a, src.a, 1);
        poly_mult_by_X_k(out.b, src.b, 1);
        return out;
    };

    size_t total_cnt_cmux = 0;
    std::vector<trlwe_lvl1<P>> ci(gr.size(), trlwe_lvl1_zero<P>()),
        co(gr.size());
    int d = input.size();
    for (int j = d - 1; j >= 0; --j) {
        auto states = gr.states_at_depth(j);
        std::for_each(std::execution::par, states.begin(), states.end(),
                      [&](auto &&q) {
                          Graph::State q1 = gr.next_state(q, true),
                                       q0 = gr.next_state(q, false);
                          auto thn = mult_X_1(ci.at(q1));
                          thn += get_w(q1);
                          auto els = mult_X_1(ci.at(q0));
                          els += get_w(q0);
                          cmux(co.at(q), input.at(j), thn, els);
                      });
        {
            using std::swap;
            swap(ci, co);
        }
        std::cerr << "[" << j << "] #CMUX : " << states.size() << "\n";
        total_cnt_cmux += states.size();
    }
    std::cerr << "Total #CMUX : " << total_cnt_cmux << "\n";

    /*
    // Cancel initial state when it's a final state
    Graph::State qi = gr.initial_state();
    auto ret = ci.at(qi);
    ret += get_w(qi);
    */
    return ci.at(gr.initial_state());
}

void det_wfa()
{
    using P = params::CGGI19;
    unsigned int seed = std::random_device{}();
    std::default_random_engine prng{seed};

    auto skey = secret_key<P>{prng};
    // auto bkey = bootstrapping_key<P>::make_ptr(prng, skey);
    // auto kskey = key_switching_key<P>::make_ptr(prng, skey);

    // "zabcde"
    std::vector<bool> plain_input = {
        false, true,  false, true,  true,  true,  true,  false, true,  false,
        false, false, false, true,  true,  false, false, true,  false, false,
        false, true,  true,  false, true,  true,  false, false, false, true,
        true,  false, false, false, true,  false, false, true,  true,  false,
        true,  false, true,  false, false, true,  true,  false
        //#include "rand_input_valid.inc"
        //#include "rand_input_invalid.inc"
    };
    std::vector<trgsw_lvl1_fft<P>> input;
    for (bool b : plain_input)
        input.emplace_back(trgsw_lvl1<P>::encrypt_bool(prng, skey, b));

    Graph gr;
    gr.reserve_states_at_depth(input.size());

    std::cout << "Input size: " << input.size() << "\n"
              << "State size: " << gr.size() << "\n"
              << "=====\n";
    trlwe_lvl1<P> enc_res = eval_det_wfa(gr, input);
    trlwe_lvl1<P>::poly_torus res = enc_res.decrypt_poly_torus(skey);
    std::cout << "Result (uint):\t" << weight2uint<P>(res) << "\n";
    std::cout << "Result (bitstr):\n" << weight2bitstring<P>(res) << "\n";
}
//////////////////////////

//////////////////////////////
//// MAIN
//////////////////////////////
int main()
{
    using namespace std::chrono;

    // Test
    auto elapsed = timeit([] { test(1, 1); });
    debug_log("Test passed. (", duration_cast<milliseconds>(elapsed).count(),
              " ms)");

    // Bench
    auto ms_per_gate = bench_hom_nand();
    debug_log("Benchmark result: ", ms_per_gate, " ms/gate");

    det_wfa();

    return 0;
}
