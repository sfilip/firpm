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
            out[i].start = cosl(in[n - i].stop);
            out[i].stop  = cosl(in[n - i].start);
            out[i].space = BandSpace::CHEBY;
        } else {
            out[i].start = acosl(in[n - i].stop);
            out[i].stop  = acosl(in[n - i].start);
            out[i].space = BandSpace::FREQ;
        }
    }
}
