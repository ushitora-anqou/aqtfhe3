#include "test.hpp"

//////////////////////////////
//// TEST
//////////////////////////////
template <class P>
void test_tlwe(unsigned int seed, const secret_key<P> &s)
{
    debug_log("test_tlwe:");

    // !!! CAVEAT !!!: std::default_random_engine is NOT
    // cryptographically secure. DO NOT use it in production!
    std::default_random_engine prng{seed};

    {
        torus m = torus_uniform(prng);
        auto tlwe = tlwe_lvl0<P>::encrypt_torus(prng, s, m);
        torus m_ = tlwe.decrypt_torus(s);
        debug_log("\t", m, " == ", m_);
        TEST_ASSERT(m == m_);
    }

    {
        std::binomial_distribution<int> dist;
        bool m = dist(prng);
        auto tlwe = tlwe_lvl0<P>::encrypt_bool(prng, s, m);
        bool m_ = tlwe.decrypt_bool(s);
        debug_log("\t", m, " == ", m_);
        TEST_ASSERT(m == m_);
    }
}

template <class P>
void test_trlwe(unsigned int seed, const secret_key<P> &s)
{
    debug_log("test_trlwe:");

    std::default_random_engine prng{seed};

    {
        auto m = random_bool_array<P::N>(prng);
        auto trlwe = trlwe_lvl1<P>::encrypt_poly_bool(prng, s, m);
        typename trlwe_lvl1<P>::poly_bool m_ = trlwe.decrypt_poly_bool(s);
        debug_log("\t[0] ", m[0], " == ", m_[0]);
        TEST_ASSERT(m == m_);
    }

    {
        auto trlwe = trlwe_lvl1<P>::encrypt_zero(prng, s);
        auto m = trlwe.decrypt_poly_torus(s);
        debug_log("\t[0] ", m[0], " == ", 0);
        for (size_t i = 0; i < P::N; i++)
            TEST_ASSERT(m[i] == 0);
    }
}

/*
template <class P>
void test_external_product(unsigned int seed, const secret_key<P> &s)
{
    debug_log("test_external_product:");
    std::default_random_engine prng{seed};

    auto m = random_bool_array<P::N>(prng);
    auto c = random_bool_value(prng);
    auto trlwe = trlwe_lvl1<P>::encrypt_poly_bool(prng, s, m);
    auto trgsw = trgsw_lvl1<P>::encrypt_bool(prng, s, c);
    trlwe_lvl1<P> res;
    external_product(res, trgsw, trlwe);
    typename trlwe_lvl1<P>::poly_bool res_plain = res.decrypt_poly_bool(s);
    debug_log("\t[0] ", c ? m[0] : 0, " == ", res_plain[0]);
    if (c)
        TEST_ASSERT(m == res_plain);
    else
        for (size_t i = 0; i < P::N; i++)
            TEST_ASSERT(res_plain[i] == 0);
}
*/

template <class P>
void test_cmux(unsigned int seed, const secret_key<P> &s)
{
    debug_log("test_cmux:");
    std::default_random_engine prng{seed};

    auto t_plain = random_bool_array<P::N>(prng);
    auto f_plain = random_bool_array<P::N>(prng);
    auto c_plain = random_bool_value(prng);
    auto t = trlwe_lvl1<P>::encrypt_poly_bool(prng, s, t_plain);
    auto f = trlwe_lvl1<P>::encrypt_poly_bool(prng, s, f_plain);
    auto c = trgsw_lvl1_fft<P>{trgsw_lvl1<P>::encrypt_bool(prng, s, c_plain)};

    trlwe_lvl1<P> res;
    cmux(res, c, t, f);
    typename trlwe_lvl1<P>::poly_bool res_plain = res.decrypt_poly_bool(s);

    debug_log("\t[0] ", c_plain, " ? ", t_plain[0], " : ", f_plain[0],
              " == ", res_plain[0]);
    TEST_ASSERT((c_plain ? t_plain : f_plain) == res_plain);
}

template <class P>
void test_blind_rotate(unsigned int seed, const secret_key<P> &s,
                       const bootstrapping_key<P> &b)
{
    debug_log("test_blind_rotate:");
    std::default_random_engine prng{seed};

    auto plain = random_bool_value(prng);
    auto tlwe = tlwe_lvl0<P>::encrypt_bool(prng, s, plain);
    constexpr auto testvec = gate_bootstrapping_test_vector<P>();

    trlwe_lvl1<P> res_trlwe;
    tlwe_lvl1<P> res_tlwe;
    blind_rotate(res_trlwe, tlwe, testvec, b);
    sample_extract_index(res_tlwe, res_trlwe, 0);
    auto res_plain = res_tlwe.decrypt_bool(s);

    debug_log("\t", plain, " == ", res_plain);
    TEST_ASSERT(plain == res_plain);
}

template <class P>
void test_identity_key_switch(unsigned int seed, const secret_key<P> &s)
{
    debug_log("test_identity_key_switch:");
    std::default_random_engine prng{seed};

    auto ks = key_switching_key<P>::make_ptr(prng, s);
    bool plain = random_bool_value(prng);
    auto tlwe1 = tlwe_lvl1<P>::encrypt_bool(prng, s, plain);
    tlwe_lvl0<P> tlwe0;
    identity_key_switch(tlwe0, tlwe1, *ks);
    bool res_plain = tlwe0.decrypt_bool(s);

    debug_log("\t", plain, " == ", res_plain);
    TEST_ASSERT(plain == res_plain);
}

template <class P>
void test_hom_nand(unsigned int seed, const secret_key<P> &s,
                   const bootstrapping_key<P> &b)
{
    debug_log("test_hom_nand:");
    std::default_random_engine prng{seed};

    auto ks = key_switching_key<P>::make_ptr(prng, s);
    bool lhs = random_bool_value(prng), rhs = random_bool_value(prng);
    auto lhs_tlwe = tlwe_lvl0<P>::encrypt_bool(prng, s, lhs),
         rhs_tlwe = tlwe_lvl0<P>::encrypt_bool(prng, s, rhs);

    tlwe_lvl0<P> res_tlwe;
    hom_nand(/* out */ res_tlwe, lhs_tlwe, rhs_tlwe, b, *ks);
    bool res = res_tlwe.decrypt_bool(s);

    debug_log("\t", "nand(", lhs, ", ", rhs, ") = ", res);
    TEST_ASSERT((!(lhs & rhs)) == res);
}

void test(size_t N, size_t M)
{
    using P = params::Test1;

    for (size_t j = 0; j < N; j++) {
        // Generate secret key and corresponding bootstrapping key
        const auto [skey, bkey] = [j] {
            unsigned int seed = std::random_device{}();
            // !!! CAVEAT !!!: std::default_random_engine is NOT
            // cryptographically secure. DO NOT use it in production!
            std::default_random_engine prng{seed};

            debug_log("Generating secret key (", j, ") (SEED: ", seed, ")");
            auto skey = secret_key<P>{prng};
            auto bkey = bootstrapping_key<P>::make_ptr(prng, skey);
            return std::make_pair(skey, bkey);
        }();

        for (size_t i = 0; i < M; i++) {
            unsigned int seed = std::random_device{}();
            debug_log("==============================");
            debug_log(j * M + i, " (SEED: ", seed, ")");

            test_tlwe<P>(seed, skey);
            test_trlwe<P>(seed, skey);
            // test_external_product<P>(seed, skey);
            test_cmux<P>(seed, skey);
            test_blind_rotate<P>(seed, skey, *bkey);
            test_identity_key_switch<P>(seed, skey);
            test_hom_nand<P>(seed, skey, *bkey);

            debug_log("==============================");
            debug_log("");
        }
    }
}
