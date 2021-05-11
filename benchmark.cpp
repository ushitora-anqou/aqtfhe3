#include "benchmark.hpp"

//////////////////////////////
//// BENCHMARK
//////////////////////////////
double bench_hom_nand()
{
    using P = params::CGGI19;
    unsigned int seed = std::random_device{}();
    std::default_random_engine prng{seed};

    auto skey = secret_key<P>{prng};
    auto bkey = bootstrapping_key<P>::make_ptr(prng, skey);
    auto kskey = key_switching_key<P>::make_ptr(prng, skey);

    // Prepare inputs and answers
    constexpr size_t N = 1000;
    std::array<bool, N> lhs_plain = random_bool_array<N>(prng),
                        rhs_plain = random_bool_array<N>(prng), ans_plain;
    for (size_t i = 0; i < N; i++)
        ans_plain[i] = !(lhs_plain[i] && rhs_plain[i]);
    std::array<tlwe_lvl0<P>, N> lhs_enc, rhs_enc, res_enc;
    for (size_t i = 0; i < N; i++) {
        lhs_enc[i] = tlwe_lvl0<P>::encrypt_bool(prng, skey, lhs_plain[i]);
        rhs_enc[i] = tlwe_lvl0<P>::encrypt_bool(prng, skey, rhs_plain[i]);
    }

    // Let's start calculation!
    debug_log("Benchmark started...");
    auto elapsed = timeit([&] {
        for (size_t i = 0; i < N; i++)
            hom_nand(res_enc[i], lhs_enc[i], rhs_enc[i], *bkey, *kskey);
    });
    debug_log("done.");

    // Check if the result is correct.
    std::array<bool, N> res_plain;
    for (size_t i = 0; i < N; i++)
        res_plain[i] = res_enc[i].decrypt_bool(skey);
    TEST_ASSERT(res_plain == ans_plain);

    using namespace std::chrono;
    return duration_cast<milliseconds>(elapsed).count() /
           static_cast<double>(N);
}
