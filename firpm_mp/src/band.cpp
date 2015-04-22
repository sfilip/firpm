#include "firpm/band.h"


void bandConversion(std::vector<Band> &out, std::vector<Band> &in,
        ConversionDirection direction, mp_prec_t prec)
{
    using mpfr::mpreal;
    mp_prec_t prevPrec = mpreal::get_default_prec();
    mpreal::set_default_prec(prec);

    out.resize(in.size());
    int n = in.size() - 1;
    for (std::size_t i = 0u; i < in.size(); ++i)
    {
        out[i].weight    = in[n - i].weight;
        out[i].amplitude = in[n - i].amplitude;
        out[i].extremas  = in[n - i].extremas;
        if (direction == ConversionDirection::FROMFREQ)
        {
            out[i].start = mpfr::cos(in[n - i].stop);
            out[i].stop  = mpfr::cos(in[n - i].start);
            out[i].space = BandSpace::CHEBY;
        } else {
            out[i].start = mpfr::acos(in[n - i].stop);
            out[i].stop  = mpfr::acos(in[n - i].start);
            out[i].space = BandSpace::FREQ;
        }
    }

    mpreal::set_default_prec(prevPrec);
}