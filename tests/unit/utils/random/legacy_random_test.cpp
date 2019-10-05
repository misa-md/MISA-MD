//
// Created by genshen on 2019/10/1.
//

#include <gtest/gtest.h>
#include <utils/random/legacy_random.hpp>
#include <xoshiro_gen.h>

class LegacyRandTest : public md_rand::LegacyRand {
    FRIEND_TEST(legacy_random_compile_test, factorial_test);

    FRIEND_TEST(legacy_random_test, factorial_test);

    friend double uniform(int &_random_seed);
};

double uniform(int &_random_seed) {
    int k = _random_seed / 127773;
    _random_seed = 16807 * (_random_seed - k * 127773) - 2836 * k;
    if (_random_seed < 0) {
        _random_seed += 0x7fffffff;
    }
    double ans = (1.0/0x7fffffff) * _random_seed;
    return ans;
}

TEST(legacy_random_compile_test, factorial_test) {
    LegacyRandTest rand1;
    rand1.seed(466953);
    std::uniform_int_distribution<> dist{0, 6};
    uint32_t num = dist(rand1);
    EXPECT_TRUE(num <= 6 && num >= 0);
}

TEST(legacy_random_test, factorial_test) {
    LegacyRandTest rand1;
    rand1.seed(466953);

    int _random_seed = 466953;

    for (int i = 0; i < 100000; i++) {
        EXPECT_EQ(rand1() * LegacyRandTest::AM, uniform(_random_seed));
    }
}
