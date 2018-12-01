//
// Created by genshen on 2018-3-12.
//

#include <gtest/gtest.h>

unsigned int Factorial(unsigned int number) {
    return number <= 1 ? number : Factorial(number - 1) * number;
}

TEST(Factorials_are_computed, factorial_test) {
    EXPECT_EQ(Factorial(1), 1);
    EXPECT_EQ(Factorial(2), 2);
    EXPECT_EQ(Factorial(3), 6);
    EXPECT_EQ(Factorial(10), 3628800);
}
