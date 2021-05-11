#ifndef AQTFHE3_UTILITY_HPP
#define AQTFHE3_UTILITY_HPP

#include <array>
#include <cassert>
#include <chrono>
#include <cstring>
#include <iostream>
#include <memory>
#include <random>

template <class G>
concept RandGen = std::uniform_random_bit_generator<G>;

template <class... Args>
void debug_log(Args &&... args)
{
    const char *envvar = std::getenv("AQTFHE3_VERBOSE");
    if (envvar != nullptr && std::strcmp(envvar, "1") == 0)
        (std::cerr << ... << args) << std::endl;
}

inline void failwith(bool cond, const char *expr, unsigned long line)
{
    if (cond)
        return;
    debug_log("\e[1;31mASSERTION FAILED! (L.", line, ") ", expr, "\e[0m");
    exit(1);
}

#define TEST_ASSERT(expr) failwith(expr, #expr, __LINE__)

template <class Proc>
auto timeit(Proc &&proc)
{
    using namespace std::chrono;
    auto begin = high_resolution_clock::now();
    proc();
    auto end = high_resolution_clock::now();
    return end - begin;
}

#endif
