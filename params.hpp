#include <cstdint>

namespace params {
struct CGGI19 {
    static constexpr uint32_t n = 630;
    // static constexpr double alpha = 3.0517578125e-05;  // 2^(-15)
    static constexpr double alpha = 0;
    static constexpr uint32_t Nbit = 10;
    static constexpr uint32_t N = 1 << Nbit;
    static constexpr uint32_t l = 3;
    static constexpr uint32_t Bgbit = 6;
    static constexpr uint32_t Bg = 1 << Bgbit;
    // static constexpr double alpha_bk = 2.98023223876953125e-08;  // 2^(-25)
    static constexpr double alpha_bk = 0;
    static constexpr uint32_t t = 8;
    static constexpr uint32_t basebit = 2;
    static constexpr double alphaks = alpha;
    static constexpr uint32_t mu = 1U << 29;
};

struct Test1 {
    static constexpr uint32_t n = 630;
    static constexpr double alpha = 0.0;
    static constexpr uint32_t Nbit = 10;
    static constexpr uint32_t N = 1 << Nbit;
    static constexpr uint32_t l = 3;
    static constexpr double alpha_bk = 0.0;
    static constexpr uint32_t Bgbit = 6;
    static constexpr uint32_t Bg = 1 << Bgbit;
    static constexpr uint32_t t = 8;
    static constexpr uint32_t basebit = 2;
};
}  // namespace params
