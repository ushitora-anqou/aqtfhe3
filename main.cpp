#include "params.hpp"

#include "spxlios.hpp"

#include <array>
#include <cassert>
#include <chrono>
#include <cstring>
#include <iostream>
#include <memory>
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

template <size_t N>
class fft_processor {
    static_assert(N >= 16, "N must be >=16");
    static_assert((N & (N - 1)) == 0, "N must be a power of 2");

private:
    spxlios::processor<N> spxlios_;

public:
    fft_processor()
    {
    }

    void twist_fft_lvl1(std::array<uint32_t, N> &out,
                        const std::array<double, N> &src)
    {
        spxlios_.execute_direct_torus32(out, src);
    }

    void twist_ifft_lvl1(std::array<double, N> &out,
                         const std::array<uint32_t, N> &src)
    {
        spxlios_.execute_reverse_torus32(out, src);
    }
};

namespace {
template <size_t N>
fft_processor<N> fftproc;
}

// Calculate out = (X^k * poly) mod (X^N + 1),
// where 0 <= k < 2 * N and length(poly) = N
template <class T, class T1, size_t N>
void poly_mult_by_X_k(poly<T, N> &out, const poly<T1, N> &poly,
                      size_t k) noexcept
{
    if (k < N) {
        for (size_t i = 0; i < N - k; i++)
            out[i + k] = poly[i];
        for (size_t i = N - k; i < N; i++)
            out[i + k - N] = -poly[i];
    }
    else {
        const size_t l = k - N;
        for (size_t i = 0; i < N - l; i++)
            out[i + l] = -poly[i];
        for (size_t i = N - l; i < N; i++)
            out[i + l - N] = poly[i];
    }
}

template <size_t N>
void mul_in_fd(std::array<double, N> &res, const std::array<double, N> &a,
               const std::array<double, N> &b)
{
    for (size_t i = 0; i < N / 2; i++) {
        double aimbim = a[i + N / 2] * b[i + N / 2];
        double arebim = a[i] * b[i + N / 2];
        res[i] = a[i] * b[i] - aimbim;
        res[i + N / 2] = a[i + N / 2] * b[i] + arebim;
    }
}

template <size_t N>
void fma_in_fd(std::array<double, N> &res, const std::array<double, N> &a,
               const std::array<double, N> &b)
{
    for (size_t i = 0; i < N / 2; i++) {
        res[i] = a[i + N / 2] * b[i + N / 2] - res[i];
        res[i] = a[i] * b[i] - res[i];
        res[i + N / 2] += a[i] * b[i + N / 2];
        res[i + N / 2] += a[i + N / 2] * b[i];
    }
}

// Calculate out = (lhs * rhs) mod (X^N + 1) using FFT,
// where length(lhs) = length(rhs) = N
template <size_t N>
void poly_mult(poly<uint32_t, N> &out, const poly<uint32_t, N> &lhs,
               const poly<uint32_t, N> &rhs)
{
    std::array<double, N> lhs_fft, rhs_fft, tmp;
    fftproc<N>.twist_ifft_lvl1(lhs_fft, lhs);
    fftproc<N>.twist_ifft_lvl1(rhs_fft, rhs);
    mul_in_fd(tmp, lhs_fft, rhs_fft);
    fftproc<N>.twist_fft_lvl1(out, tmp);
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
template <class P>
using tlwe_lvl1 = impl::tlwe<P, 1>;

template <class P>
struct trlwe_lvl1 {
    using poly_torus = poly<torus, P::N>;
    using poly_bool = poly<bool, P::N>;

    poly_torus a, b;

    template <RandGen RG>
    static trlwe_lvl1 encrypt_poly_torus(RG &rng, const secret_key<P> &skey,
                                         const poly_torus &m)
    {
        trlwe_lvl1 trlwe;
        poly_torus &a = trlwe.a, &b = trlwe.b;
        const poly_torus &s = skey.lvl1;

        // a[X] <- U_{T_N[X]}
        for (torus &t : a)
            t = torus_uniform(rng);

        // e[X] <- D_{T_N[X], alpha_bk}
        poly_torus e;
        for (torus &t : e)
            t = torus_modular_normal(rng, P::alpha_bk);

        // b[X] = a[X] * s[X] + m[X] + e[X]
        poly_mult(/* out */ b, a, s);
        for (size_t i = 0; i < P::N; i++)
            b[i] += m[i] + e[i];

        return trlwe;
    }

    template <RandGen RG>
    static trlwe_lvl1 encrypt_zero(RG &rng, const secret_key<P> &skey)
    {
        poly_torus m = {};  // Initialize with 0
        return encrypt_poly_torus(rng, skey, m);
    }

    template <RandGen RG>
    static trlwe_lvl1 encrypt_poly_bool(RG &rng, const secret_key<P> &skey,
                                        const poly_bool &m)
    {
        const torus mu = 1u << 29;  // 1/8
        poly_torus t;
        for (size_t i = 0; i < P::N; i++)
            t[i] = m[i] ? mu : -mu;
        return encrypt_poly_torus(rng, skey, t);
    }

    // Generate trivial ciphertext
    static constexpr trlwe_lvl1 trivial_encrypt_poly_torus(const poly_torus &m)
    {
        trlwe_lvl1 trlwe;
        for (size_t i = 0; i < P::N; i++)
            trlwe.a[i] = 0;
        trlwe.b = m;
        return trlwe;
    }

    poly_torus decrypt_poly_torus(const secret_key<P> &skey) const noexcept
    {
        const poly_torus &a = this->a, &b = this->b, &s = skey.lvl1;

        // m[X] = b[X] - a[X] * s[X] - e[X]
        poly_torus m, as;
        poly_mult(/* out */ as, a, s);
        for (size_t i = 0; i < P::N; i++)
            m[i] = b[i] - as[i];

        return m;
    }

    poly_bool decrypt_poly_bool(const secret_key<P> &skey) const noexcept
    {
        poly_torus t = decrypt_poly_torus(skey);
        poly_bool m;
        for (size_t i = 0; i < P::N; i++)
            m[i] = static_cast<int32_t>(t[i]) > 0 ? true : false;
        return m;
    }
};

template <class P>
struct trgsw_lvl1 {
    std::array<trlwe_lvl1<P>, 2 * P::l> data;

    trlwe_lvl1<P> &operator[](size_t i)
    {
        return data[i];
    }

    const trlwe_lvl1<P> &operator[](size_t i) const
    {
        return data[i];
    }

    template <RandGen RG>
    static trgsw_lvl1 encrypt_poly_int(RG &rng, const secret_key<P> &skey,
                                       const poly<int, P::N> &mu)
    {
        // Assume Bg = (1 << Bgbit)
        constexpr size_t N = P::N, l = P::l, Bgbit = P::Bgbit;

        trgsw_lvl1 trgsw;
        // trgsw[i] = (a[i][X], b[i][X])
        for (auto &trlwe : trgsw.data)
            trlwe = trlwe_lvl1<P>::encrypt_zero(rng, skey);
        // trgsw[i].a += mu[X] / (Bg^i)  (0 <= i < l)
        // trgsw[l + i].b += mu[X] / (Bg^i)  (0 <= i < l)
        for (size_t i = 0; i < l; i++) {
            for (size_t j = 0; j < N; j++) {
                torus t =
                    static_cast<torus>(mu[j]) * (1u << (32 - (i + 1) * Bgbit));
                trgsw[i].a[j] += t;
                trgsw[l + i].b[j] += t;
            }
        }

        return trgsw;
    }

    template <RandGen RG>
    static trgsw_lvl1 encrypt_bool(RG &rng, const secret_key<P> &skey, bool mu)
    {
        poly<int, P::N> t = {};  // Initialize with 0
        t[0] = mu ? 1 : 0;
        return encrypt_poly_int(rng, skey, t);
    }
};

template <class P>
struct trgsw_lvl1_fft {
    struct body {
        std::array<double, P::N> a, b;
    };
    std::array<body, 2 * P::l> data;

    trgsw_lvl1_fft()
    {
    }

    trgsw_lvl1_fft(const trgsw_lvl1<P> &src)
    {
        for (size_t i = 0; i < 2 * P::l; i++) {
            fftproc<P::N>.twist_ifft_lvl1(data[i].a, src[i].a);
            fftproc<P::N>.twist_ifft_lvl1(data[i].b, src[i].b);
        }
    }

    const body &operator[](size_t i) const noexcept
    {
        return data[i];
    }
};

template <class P>
struct bootstrapping_key {
    std::array<trgsw_lvl1_fft<P>, P::n> data;

    bootstrapping_key()
    {
    }

    bootstrapping_key(const bootstrapping_key &that)
    {
        this->data = that.data;
    }

    template <RandGen RG>
    bootstrapping_key(RG &rng, const secret_key<P> &skey)
    {
        // Encrypt every bit of secret key as TRGSWlvl1
        for (size_t i = 0; i < P::n; i++) {
            auto trgsw = trgsw_lvl1<P>::encrypt_bool(rng, skey, skey.lvl0[i]);
            data[i] = trgsw_lvl1_fft<P>{trgsw};
        }
    }

    // NOTE: struct bootstrapping_key needs large space of memory. Allocating it
    // on stack may cause segmentaion fault.
    template <RandGen RG>
    static std::shared_ptr<bootstrapping_key> make_ptr(
        RG &rng, const secret_key<P> &skey)
    {
        return std::make_shared<bootstrapping_key>(rng, skey);
    }

    const trgsw_lvl1_fft<P> &operator[](size_t i) const noexcept
    {
        return data[i];
    }
};

template <class P>
struct key_switching_key {
    std::array<tlwe_lvl0<P>, P::N * P::t *((1 << P::basebit) - 1)> data;

    tlwe_lvl0<P> &operator()(size_t i, size_t j, size_t k)
    {
        constexpr size_t t = P::t, base = (1 << P::basebit);
        return data[k + (base - 1) * (j + t * i)];
    }

    const tlwe_lvl0<P> &operator()(size_t i, size_t j, size_t k) const
    {
        constexpr size_t t = P::t, base = (1 << P::basebit);
        return data[k + (base - 1) * (j + t * i)];
    }

    template <RandGen RG>
    static std::shared_ptr<key_switching_key> make_ptr(
        RG &rng, const secret_key<P> &skey)
    {
        constexpr size_t N = P::N, t = P::t, basebit = P::basebit,
                         base = (1 << basebit);

        auto ks = std::make_shared<key_switching_key>();
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < t; j++) {
                for (size_t k = 1; k <= base - 1; k++) {  // Ignore k = 0
                    torus p =
                        k * skey.lvl1[i] * (1u << (32 - basebit * (j + 1)));
                    (*ks)(i, j, k - 1) =
                        tlwe_lvl0<P>::encrypt_torus(rng, skey, p);
                }
            }
        }

        return ks;
    }
};

template <class P>
void sample_extract_index(tlwe_lvl1<P> &out, const trlwe_lvl1<P> &trlwe,
                          size_t k) noexcept
{
    constexpr size_t N = P::N;
    assert(k < N);

    out.b() = trlwe.b[k];
    for (size_t i = 0; i < N; i++)
        if (i <= k)
            out.a(i) = trlwe.a[k - i];
        else
            out.a(i) = -trlwe.a[N + k - i];
}

template <class P>
void decompose(std::array<poly<torus, P::N>, P::l> &out,
               const poly<torus, P::N> &a)
{
    constexpr torus Bg = P::Bg;
    constexpr size_t l = P::l, N = P::N, Bgbit = P::Bgbit;

    torus offset = 0;
    for (size_t i = 0; i < l; i++)
        offset += Bg / 2 * (1u << (32 - (i + 1) * Bgbit));

    torus round_offset = 1ULL << (32 - l * Bgbit - 1);

    poly<torus, P::N> a_tilde;
    for (size_t i = 0; i < N; i++)
        a_tilde[i] = a[i] + offset + round_offset;

    for (size_t i = 0; i < l; i++)
        for (size_t j = 0; j < N; j++)
            out[i][j] =
                ((a_tilde[j] >> (32 - Bgbit * (i + 1))) & (Bg - 1)) - Bg / 2;
}

template <class P>
void external_product(trlwe_lvl1<P> &out, const trgsw_lvl1_fft<P> &trgsw_fft,
                      const trlwe_lvl1<P> &trlwe)
{
    constexpr size_t l = P::l, N = P::N;

    // Decompose trlwe
    std::array<poly<torus, N>, l> dec_a, dec_b;
    decompose<P>(dec_a, trlwe.a);
    decompose<P>(dec_b, trlwe.b);

    // Apply inverse FFT to decomposed trlwe
    std::array<std::array<double, N>, l> dec_a_fft, dec_b_fft;
    for (size_t i = 0; i < P::l; i++) {
        fftproc<N>.twist_ifft_lvl1(dec_a_fft[i], dec_a[i]);
        fftproc<N>.twist_ifft_lvl1(dec_b_fft[i], dec_b[i]);
    }

    // Do multiplication
    std::array<double, N> out_a_fft = {}, out_b_fft = {};  // Initialize with 0
    for (size_t i = 0; i < l; i++) {
        fma_in_fd(out_a_fft, dec_a_fft[i], trgsw_fft[i].a);
        fma_in_fd(out_a_fft, dec_b_fft[i], trgsw_fft[l + i].a);
        fma_in_fd(out_b_fft, dec_a_fft[i], trgsw_fft[i].b);
        fma_in_fd(out_b_fft, dec_b_fft[i], trgsw_fft[l + i].b);
    }

    // Apply FFT
    fftproc<N>.twist_fft_lvl1(out.a, out_a_fft);
    fftproc<N>.twist_fft_lvl1(out.b, out_b_fft);
}

template <class P>
void cmux(trlwe_lvl1<P> &out, const trgsw_lvl1_fft<P> &cond,
          const trlwe_lvl1<P> &thn, const trlwe_lvl1<P> &els)
{
    const trlwe_lvl1<P> &trlwe0 = els, &trlwe1 = thn;

    // tmp0 = (a1[X], b1[X]) - (a0[X], b0[X])
    trlwe_lvl1<P> tmp0;
    for (size_t i = 0; i < P::N; i++) {
        tmp0.a[i] = trlwe1.a[i] - trlwe0.a[i];
        tmp0.b[i] = trlwe1.b[i] - trlwe0.b[i];
    }

    // tmp1 = external_product(cond, tmp0)
    trlwe_lvl1<P> tmp1;
    external_product(/* out */ tmp1, cond, tmp0);

    // out = tmp1 + trlwe0
    for (size_t i = 0; i < P::N; i++) {
        out.a[i] = tmp1.a[i] + trlwe0.a[i];
        out.b[i] = tmp1.b[i] + trlwe0.b[i];
    }
}

template <class P>
void blind_rotate(trlwe_lvl1<P> &out, const tlwe_lvl0<P> &src,
                  const trlwe_lvl1<P> &testvec,
                  const bootstrapping_key<P> &bkey)
{
    constexpr size_t n = P::n, N = P::N, Nbit = P::Nbit;
    const size_t b_tilda =
        2 * N - ((src.b() + (1u << (31 - Nbit - 1))) >> (32 - Nbit - 1));

    // Initialize out = X^{b_tilda} * (a[X], b[X])
    for (size_t i = 0; i < N; i++)
        out.a[i] = out.b[i] = 0;
    poly_mult_by_X_k(out.a, testvec.a, b_tilda);
    poly_mult_by_X_k(out.b, testvec.b, b_tilda);

    for (size_t i = 0; i < n; i++) {
        const size_t a_tilda =
            (src.a(i) + (1 << (31 - Nbit - 1))) >> (32 - Nbit - 1);
        const trlwe_lvl1<P> trlwe0 = out;

        // Let trlwe1 = X^{a_tilda} * trlwe0
        trlwe_lvl1<P> trlwe1 = {};  // Initialize with 0.
        poly_mult_by_X_k(trlwe1.a, trlwe0.a, a_tilda);
        poly_mult_by_X_k(trlwe1.b, trlwe0.b, a_tilda);

        // out = dec(bkey[i]) ? (X^{a_tilda} * out) : out
        cmux(out, bkey[i], trlwe1, trlwe0);
    }
}

template <class P>
constexpr trlwe_lvl1<P> gate_bootstrapping_test_vector()
{
    constexpr torus mu = 1u << 29;  // 1/8
    typename trlwe_lvl1<P>::poly_torus m;
    for (size_t i = 0; i < P::N; i++)
        m[i] = mu;
    return trlwe_lvl1<P>::trivial_encrypt_poly_torus(m);
}

template <class P>
void identity_key_switch(tlwe_lvl0<P> &out, const tlwe_lvl1<P> &src,
                         const key_switching_key<P> &ks)
{
    constexpr size_t n = P::n, N = P::N, bb = P::basebit, t = P::t;

    out.b() = src.b();
    for (size_t i = 0; i < n; i++)
        out.a(i) = 0;

    size_t prec_offset = 1u << (32 - (1 + bb * t));
    for (size_t i = 0; i < N; i++) {
        size_t abar = src.a(i) + prec_offset;
        for (size_t j = 0; j < t; j++) {
            size_t k = (abar >> (32 - (j + 1) * bb)) & ((1u << bb) - 1);
            if (k == 0)
                continue;

            // (a_tilda, b_tilda) -= KS_ijk
            auto &ks_ijk = ks(i, j, k - 1);
            for (size_t l = 0; l < n; l++)
                out.a(l) -= ks_ijk.a(l);
            out.b() -= ks_ijk.b();
        }
    }
}

template <class P>
void gate_bootstrapping_tlwe_to_tlwe(tlwe_lvl1<P> &out, const tlwe_lvl0<P> &src,
                                     const bootstrapping_key<P> &bkey)
{
    constexpr auto testvec = gate_bootstrapping_test_vector<P>();

    trlwe_lvl1<P> trlwe;
    blind_rotate(/* out */ trlwe, src, testvec, bkey);
    sample_extract_index(/* out */ out, trlwe, 0);
}

template <class P>
void hom_nand(tlwe_lvl0<P> &out, const tlwe_lvl0<P> &lhs,
              const tlwe_lvl0<P> &rhs, const bootstrapping_key<P> &bkey,
              const key_switching_key<P> &ks)
{
    // out = ((0, 1/8) - lhs - rhs) |> gate bootstrapping |> identity key switch
    constexpr size_t n = P::n;
    const torus mu = 1u << 29;  // 1/8

    tlwe_lvl0<P> tlwe;
    for (size_t i = 0; i < n; i++)
        tlwe.a(i) = -lhs.a(i) - rhs.a(i);
    tlwe.b() = mu - lhs.b() - rhs.b();

    tlwe_lvl1<P> tlwe1;
    gate_bootstrapping_tlwe_to_tlwe(/* out */ tlwe1, tlwe, bkey);
    identity_key_switch(/* out */ out, tlwe1, ks);
}

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

    return 0;
}
