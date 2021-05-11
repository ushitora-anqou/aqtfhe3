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
    // When on state `index`,
    //   input 0 -> next state is `child0`,
    //              next weight is (w_{org} * X^k + gain0)
    //   input 1 -> next state is `child1`
    //              next weight is (w_{org} * X^k + gain1)
    struct TableItem {
        State index, child0, child1;
        uint64_t gain0, gain1;
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

    uint64_t gain(State from, bool input) const
    {
        auto &s = table_.at(from);
        return input ? s.gain1 : s.gain0;
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
class DetWFARunner {
private:
    Graph graph_;
    std::vector<trgsw_lvl1_fft<P>> input_;
    std::vector<trlwe_lvl1<P>> weight_;
    bool hasEvaluated_;
    int shiftWidth_, shiftInterval_;

public:
    DetWFARunner(Graph graph, std::vector<trgsw_lvl1_fft<P>> input)
        : graph_(std::move(graph)),
          input_(std::move(input)),
          weight_(graph_.size(), trlwe_lvl1_zero<P>()),
          hasEvaluated_(false),
          shiftWidth_(1),
          shiftInterval_(8)
    {
    }

    const trlwe_lvl1<P> &result() const
    {
        assert(hasEvaluated_);
        return weight_.at(graph_.initial_state());
    }

    void eval()
    {
        assert(!hasEvaluated_);
        hasEvaluated_ = true;

        size_t total_cnt_cmux = 0;
        std::vector<trlwe_lvl1<P>> out(graph_.size());
        for (int j = input_.size() - 1; j >= 0; --j) {
            auto states = graph_.states_at_depth(j);
            std::for_each(std::execution::par, states.begin(), states.end(),
                          [&](auto &&q) {
                              trlwe_lvl1<P> w0, w1;
                              next_weight(w1, j, q, true);
                              next_weight(w0, j, q, false);
                              cmux(out.at(q), input_.at(j), w1, w0);
                          });
            {
                using std::swap;
                swap(out, weight_);
            }
            std::cerr << "[" << j << "] #CMUX : " << states.size() << "\n";
            total_cnt_cmux += states.size();
        }
        std::cerr << "Total #CMUX : " << total_cnt_cmux << "\n";
    }

private:
    void mult_X_k(trlwe_lvl1<P> &out, const trlwe_lvl1<P> &src, size_t k) const
    {
        poly_mult_by_X_k(out.a, src.a, k);
        poly_mult_by_X_k(out.b, src.b, k);
    }

    void next_weight(trlwe_lvl1<P> &out, int j, Graph::State from,
                     bool input) const
    {
        Graph::State to = graph_.next_state(from, input);
        uint64_t w = graph_.gain(from, input);

        // Return Enc(
        //   w_{to} * X^{shiftWidth} +
        //   w_0 * X^0 +...+ w_{N-1} * X^{N-1}
        // )
        if (j % shiftInterval_ == shiftInterval_ - 1)
            mult_X_k(out, weight_.at(to), shiftWidth_);
        else
            out = weight_.at(to);
        out += trlwe_lvl1<P>::trivial_encrypt_poly_torus(uint2weight<P>(w));
    }
};

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
    // trlwe_lvl1<P> enc_res = eval_det_wfa(gr, input);

    DetWFARunner<P> runner{gr, input};
    runner.eval();
    trlwe_lvl1<P> enc_res = runner.result();

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

    /*
    // Test
    auto elapsed = timeit([] { test(1, 1); });
    debug_log("Test passed. (", duration_cast<milliseconds>(elapsed).count(),
              " ms)");

    // Bench
    auto ms_per_gate = bench_hom_nand();
    debug_log("Benchmark result: ", ms_per_gate, " ms/gate");
    */

    det_wfa();

    return 0;
}
