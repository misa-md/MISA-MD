//
// Created by genshen on 2019/10/1.
//

#ifndef CRYSTAL_MD_LEGACY_RANDOM_HPP
#define CRYSTAL_MD_LEGACY_RANDOM_HPP

namespace md_rand {
    class LegacyRand {
    public:
        void seed(const uint32_t seed) {
            _random_seed = seed;
        }

        uint32_t rand() {
            int k = _random_seed / IQ;
            _random_seed = IA * (_random_seed - k * IQ) - IR * k;
            if (_random_seed < 0) {
                _random_seed += IM;
            }
            double ans = AM * _random_seed;
            return ans;
        }

    private:
        uint32_t _random_seed;
        static const uint32_t IA = 16807;
        static const uint32_t IM = 2147483647;
        static const uint32_t AM = (1.0 / IM);
        static const uint32_t IQ = 127773;
        static const uint32_t IR = 2836;
    };
}

#endif //CRYSTAL_MD_LEGACY_RANDOM_HPP
