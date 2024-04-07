//    firpm
//    Copyright (C) 2015 - 2024  S. Filip

#include "firpm/band.h"
#include "firpm/pmmath.h"

namespace pm {

    template<typename T>
    void bandconv(std::vector<band_t<T>>& out, std::vector<band_t<T>>& in,
            convdir_t direction)
    {
        out.resize(in.size());
        std::size_t n = in.size() - 1u;
        for (std::size_t i{0u}; i < in.size(); ++i)
        {
            out[i].weight    = in[n - i].weight;
            out[i].amplitude = in[n - i].amplitude;
            out[i].xs        = in[n - i].xs;
            out[i].part      = in[n - i].part;
            if (direction == convdir_t::FROMFREQ)
            {
                out[i].start = pmmath::cos(in[n - i].stop);
                out[i].stop  = pmmath::cos(in[n - i].start);
                out[i].space = space_t::CHEBY;
                for(std::size_t j{0}; j < out[i].part.size(); ++j)
                    out[i].part[j] = pmmath::cos(out[i].part[j]);
                std::sort(begin(out[i].part), end(out[i].part));
            } else {
                out[i].start = pmmath::acos(in[n - i].stop);
                out[i].stop  = pmmath::acos(in[n - i].start);
                out[i].space = space_t::FREQ;
                for(std::size_t j{0}; j < out[i].part.size(); ++j)
                    out[i].part[j] = pmmath::acos(out[i].part[j]);
                std::sort(begin(out[i].part), end(out[i].part));
            }
        }
    }

    /* Template instantiation */
    template void bandconv<double>(
        std::vector<band_t<double>>& out,
        std::vector<band_t<double>>& in,
            convdir_t direction);

    template void bandconv<long double>(
        std::vector<band_t<long double>>& out,
        std::vector<band_t<long double>>& in,
            convdir_t direction);

#ifdef HAVE_MPFR
    template void bandconv<mpfr::mpreal>(
        std::vector<band_t<mpfr::mpreal>>& out,
        std::vector<band_t<mpfr::mpreal>>& in,
            convdir_t direction);
#endif

} // namespace pm