#include "params.hpp"

#include <array>
#include <cassert>
#include <chrono>
#include <cstring>
#include <iostream>
#include <random>

template <class G>
concept RandGen = std::uniform_random_bit_generator<G>;

using torus = uint32_t;
static_assert(std::numeric_limits<torus>::digits == 32);

template <class T, size_t Size>
using poly = std::array<T, Size>;

template <class... Args>
void debug_log(Args &&... args)
{
    const char *envvar = std::getenv("AQTFHE3_VERBOSE");
    if (envvar != nullptr && std::strcmp(envvar, "1") == 0)
        (std::cerr << ... << args) << std::endl;
}

torus double2torus(double src)
{
    return static_cast<torus>(
        std::fmod(src, 1.0) *
        std::pow(2.0, std::numeric_limits<torus>::digits));
}

template <RandGen RG>
torus torus_uniform(RG &rng)
{
    std::uniform_int_distribution<torus> dist{
        0, std::numeric_limits<torus>::max()};
    return dist(rng);
}

template <RandGen RG>
torus torus_modular_normal(RG &rng, double alpha)
{
    std::normal_distribution<> dist{0.0, alpha};
    return double2torus(dist(rng));
}

template <class P>
struct secret_key {
    // Use type torus instead of bool for FFT
    std::array<torus, P::n> lvl0;
    std::array<torus, P::N> lvl1;

    secret_key()
    {
    }

    secret_key(const secret_key &that)
    {
        this->lvl0 = that.lvl0;
        this->lvl1 = that.lvl1;
    }

    template <RandGen RG>
    secret_key(RG &rng)
    {
        std::binomial_distribution<int> dist;
        dist(rng);
        for (torus &t : lvl0)  // s <- B_n
            t = dist(rng);
        for (torus &t : lvl1)  // s[X] <- B_N[X]
            t = dist(rng);
    }
};

namespace impl {
template <class P, int lvl>
struct tlwe {
    static_assert(lvl == 0 || lvl == 1);

    static constexpr size_t N() noexcept
    {
        if constexpr (lvl == 0)
            return P::n;
        else
            return P::N;
    }

    static constexpr double stddev() noexcept
    {
        if constexpr (lvl == 0)
            return P::alpha;
        else
            return P::alpha_bk;
    }

    static const std::array<torus, N()> &key(const secret_key<P> &skey) noexcept
    {
        if constexpr (lvl == 0)
            return skey.lvl0;
        else
            return skey.lvl1;
    }

    std::array<torus, N() + 1> data;

    torus &operator[](size_t i) noexcept
    {
        return data[i];
    }

    torus operator[](size_t i) const noexcept
    {
        return data[i];
    }

    torus &a(size_t i) noexcept
    {
        assert(i < N());
        return (*this)[i];
    }

    torus a(size_t i) const noexcept
    {
        assert(i < N());
        return (*this)[i];
    }

    torus &b() noexcept
    {
        return (*this)[N()];
    }

    torus b() const noexcept
    {
        return (*this)[N()];
    }

    template <RandGen RG>
    static tlwe encrypt_torus(RG &rng, const secret_key<P> &skey, torus m)
    {
        tlwe tlwe;
        const auto &s = key(skey);

        // a[i] <- U_{T}
        for (size_t i = 0; i < N(); i++)
            tlwe.a(i) = torus_uniform(rng);

        // b = m + e + a * s
        torus e = torus_modular_normal(rng, stddev());  // e <- D_{T, alpha}
        tlwe.b() = m + e;
        for (size_t i = 0; i < N(); i++)
            tlwe.b() += tlwe.a(i) * s[i];

        return tlwe;
    }

    template <RandGen RG>
    static tlwe encrypt_bool(RG &rng, const secret_key<P> &skey, bool m)
    {
        const torus mu = 1u << 29;  // 1/8
        return encrypt_torus(rng, skey, m ? mu : -mu);
    }

    torus decrypt_torus(const secret_key<P> &skey) const
    {
        const auto &s = key(skey);

        // m = b - a * s - e
        torus m = b();
        for (size_t i = 0; i < N(); i++)
            m -= a(i) * s[i];

        return m;
    }

    bool decrypt_bool(const secret_key<P> &s) const
    {
        return static_cast<int32_t>(decrypt_torus(s)) > 0;
    }
};
}  // namespace impl

template <class P>
using tlwe_lvl0 = impl::tlwe<P, 0>;

//////////////////////////////
//// TEST
//////////////////////////////
template <size_t N, RandGen RG>
poly<bool, N> random_bool_array(RG &rng)
{
    poly<bool, N> ret;
    std::binomial_distribution<int> dist;
    for (bool &v : ret)
        v = dist(rng);
    return ret;
}

template <RandGen RG>
bool random_bool_value(RG &rng)
{
    return random_bool_array<1>(rng)[0];
}

void failwith(bool cond, const char *expr, unsigned long line)
{
    if (cond)
        return;
    debug_log("\e[1;31mASSERTION FAILED! (L.", line, ") ", expr, "\e[0m");
    exit(1);
}

#define TEST_ASSERT(expr) failwith(expr, #expr, __LINE__)

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

template <class Proc>
auto timeit(Proc &&proc)
{
    using namespace std::chrono;
    auto begin = high_resolution_clock::now();
    proc();
    auto end = high_resolution_clock::now();
    return end - begin;
}

void test(size_t N, size_t M)
{
    using P = params::Test1;

    for (size_t j = 0; j < N; j++) {
        // Generate secret key and corresponding bootstrapping key
        const auto skey = [j] {
            unsigned int seed = std::random_device{}();
            // !!! CAVEAT !!!: std::default_random_engine is NOT
            // cryptographically secure. DO NOT use it in production!
            std::default_random_engine prng{seed};

            debug_log("Generating secret key (", j, ") (SEED: ", seed, ")");
            auto skey = secret_key<P>{prng};
            return skey;
        }();

        for (size_t i = 0; i < M; i++) {
            unsigned int seed = std::random_device{}();
            debug_log("==============================");
            debug_log(j * M + i, " (SEED: ", seed, ")");

            test_tlwe<P>(seed, skey);

            debug_log("==============================");
            debug_log("");
        }
    }
}

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

    return 0;
}
