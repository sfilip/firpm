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

template<typename T>
void bandconv(std::vector<band_t<T>> &out, std::vector<band_t<T>> &in,
        convdir_t direction)
{
    out.resize(in.size());
    std::size_t n = in.size() - 1u;
    for (std::size_t i{0u}; i < in.size(); ++i)
    {
        out[i].weight    = in[n - i].weight;
        out[i].amplitude = in[n - i].amplitude;
        out[i].xs        = in[n - i].xs;
        if (direction == convdir_t::FROMFREQ)
        {
            out[i].start = cosl(in[n - i].stop);
            out[i].stop  = cosl(in[n - i].start);
            out[i].space = space_t::CHEBY;
        } else {
            out[i].start = acosl(in[n - i].stop);
            out[i].stop  = acosl(in[n - i].start);
            out[i].space = space_t::FREQ;
        }
    }
}

/* Template instantiation */
template void bandconv<double>(
	std::vector<band_t<double>> &out,
	std::vector<band_t<double>> &in,
        convdir_t direction);

template void bandconv<long double>(
	std::vector<band_t<long double>> &out,
	std::vector<band_t<long double>> &in,
        convdir_t direction);
