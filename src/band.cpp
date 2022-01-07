//    firpm_d
//    Copyright (C) 2015 - 2019  S. Filip
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "firpm/band.h"
#include "firpm/pmmath.h"

namespace pm {

    template<typename T>
    std::vector<band_t<T>> bandconv(std::vector<band_t<T>>& in,
            convdir_t direction)
    {
        std::vector<band_t<T>> out(in.size());
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
        return out;
    }

    /* Template instantiation */
    template std::vector<band_t<double>> bandconv<double>(
        std::vector<band_t<double>>& in,
            convdir_t direction);

    template std::vector<band_t<long double>> bandconv<long double>(
        std::vector<band_t<long double>>& in,
            convdir_t direction);

#ifdef HAVE_MPFR
    template std::vector<band_t<mpfr::mpreal>> bandconv<mpfr::mpreal>(
        std::vector<band_t<mpfr::mpreal>>& in,
            convdir_t direction);
#endif

} // namespace pm
