//    firpm_d
//    Copyright (C) 2015  S. Filip
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

void bandConversion(std::vector<Band> &out, std::vector<Band> &in,
        ConversionDirection direction)
{
    out.resize(in.size());
    int n = in.size() - 1;
    for (std::size_t i = 0u; i < in.size(); ++i)
    {
        out[i].weight    = in[n - i].weight;
        out[i].amplitude = in[n - i].amplitude;
        out[i].extremas  = in[n - i].extremas;
        if (direction == ConversionDirection::FROMFREQ)
        {
            out[i].start = cos(in[n - i].stop);
            out[i].stop  = cos(in[n - i].start);
            out[i].space = BandSpace::CHEBY;
        } else {
            out[i].start = acos(in[n - i].stop);
            out[i].stop  = acos(in[n - i].start);
            out[i].space = BandSpace::FREQ;
        }
    }
}
