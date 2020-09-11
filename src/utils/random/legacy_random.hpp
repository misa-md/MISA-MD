//
// Created by genshen on 2019/10/1.
//

#ifndef MISA_MD_LEGACY_RANDOM_HPP
#define MISA_MD_LEGACY_RANDOM_HPP

namespace md_rand {
    class LegacyRand {
    public:
        typedef uint32_t result_type;
        typedef int32_t seed_type;

        void seed(const seed_type seed) {
            _random_seed = seed;
        }

        result_type operator()() {
            seed_type k = _random_seed / IQ;
            _random_seed = IA * (_random_seed - k * IQ) - IR * k;
            if (_random_seed < 0) {
                _random_seed += IM;
            }
            return _random_seed;
        }

        static inline constexpr result_type max() {
            return IM;
        }

        static inline constexpr result_type min() {
            return 0;
        }

    protected:
        seed_type _random_seed;
        static const seed_type IA = 16807;
        static const seed_type IM = 0x7fffffff;
        static constexpr double AM = (1.0 / IM);
        static const seed_type IQ = 127773;
        static const seed_type IR = 2836;
    };
}

#endif //MISA_MD_LEGACY_RANDOM_HPP
