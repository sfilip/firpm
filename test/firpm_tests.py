import numpy as np
import pytest

import firpm


@pytest.mark.parametrize("numtaps, fpass, fstop, ripple, attenuation", [
    # Specification 4
    (1001, 0.8, 0.81, 0.1, -90),
    (2002, 0.8, 0.81, 0.1, -90),

    # Specification 5
    (2003, 0.1, 0.105, 0.1, -90),
])
def test_lpf_spectrum(numtaps, fpass, fstop, ripple, attenuation):
    # Around 2000 taps the uniform initial reference leaves the leveled error
    # below double-precision noise, so whether the exchange converges depends
    # on libm rounding. Use reference scaling there, as the C++ tests do.
    kwargs = dict(strategy=firpm.Strategy.SCALING, depth=1) if numtaps > 2000 else {}
    h = firpm.firpm(numtaps - 1, [0.0, fpass, fstop, 1.0],
                    [1.0, 1.0, 0.0, 0.0], [1.0, 10.0], **kwargs)

    # Spectrum in dB, against frequency normalized to Nyquist
    hdb = 20 * np.log10(np.abs(np.fft.fft(h)))
    f = np.linspace(0, 2, len(hdb))

    assert np.max(np.abs(hdb[f < fpass])) <= ripple
    assert np.max(hdb[(f > fstop) & (f < 1)]) <= attenuation


FLAT = ([0.0, 0.2, 0.25, 0.5], [1.0, 0.0])
FLAT_WEIGHT = [1.0, 10.0]


def test_matches_scipy_remez():
    signal = pytest.importorskip("scipy.signal")
    for numtaps in [51, 52]:
        h = firpm.remez(numtaps, *FLAT, weight=FLAT_WEIGHT)
        h_scipy = signal.remez(numtaps, *FLAT, weight=FLAT_WEIGHT)
        assert len(h) == numtaps
        assert np.allclose(h, h_scipy, atol=1e-4)


def test_returns_ndarray():
    for h in [firpm.remez(51, *FLAT),
              firpm.remez(51, [(0.0, 0.2, lambda f: 1.0), (0.25, 0.5, 0.0)]),
              firpm.firpm(50, [0.0, 0.4, 0.5, 1.0], [1.0, 1.0, 0.0, 0.0]),
              firpm.firpm(30, [0.1, 0.9], [1.0, 1.0], "hilbert")]:
        assert isinstance(h, np.ndarray) and h.dtype == np.float64


def test_version():
    assert firpm.__version__.count(".") == 2


# MATLAB's firpm(11, [0 .3 .5 1], [1 1 0 0]) and firpm(20, [0 .8 .9 1],
# [0 0 1 1]), from scipy's TestRemez.test_compare. MATLAB uses a frequency
# grid, so its designs differ from the grid-free optimum by about 6e-5.
MATLAB_12 = [0.024590270518440, -0.041314581814658, -0.075943803756711,
             -0.003530911231040, 0.193140296954975, 0.373400753484939]
MATLAB_21 = [-0.038976016082299, 0.018704846485491, -0.014644062687875,
             0.002879152556419, 0.016849978528150, -0.043276706138248,
             0.073641298245579, -0.103908158578635, 0.129770906801075,
             -0.147163447297124, 0.153302248456347]


def peak_error(h, f, a):
    # The largest error over the bands (edges normalized to Nyquist, unit
    # weights), on a dense grid.
    H = np.abs(np.fft.rfft(h, 2**18))
    x = np.linspace(0, 1, len(H))
    peak = 0.0
    for i in range(0, len(f), 2):
        m = (x >= f[i]) & (x <= f[i + 1])
        D = a[i] + (a[i + 1] - a[i]) * (x[m] - f[i]) / (f[i + 1] - f[i])
        peak = max(peak, np.max(np.abs(H[m] - D)))
    return peak


@pytest.mark.parametrize("n, f, a, matlab", [
    (11, [0, 0.3, 0.5, 1], [1, 1, 0, 0], MATLAB_12 + MATLAB_12[::-1]),
    (20, [0, 0.8, 0.9, 1], [0, 0, 1, 1], MATLAB_21 + MATLAB_21[-2::-1]),
])
def test_firpm_no_worse_than_matlab(n, f, a, matlab):
    # Close to MATLAB's design, and at least as close to the minimax optimum.
    h = firpm.firpm(n, f, a)
    assert np.allclose(h, matlab, atol=1e-4)
    assert peak_error(h, f, a) <= peak_error(np.array(matlab), f, a)


@pytest.mark.parametrize("numtaps", [51, 52])
def test_firpm_matches_remez(numtaps):
    # The same design in both conventions: order vs. taps, edges normalized
    # to Nyquist vs. in units of fs, amplitudes per edge vs. per band.
    h_firpm = firpm.firpm(numtaps - 1, [0.0, 0.4, 0.5, 1.0],
                          [1.0, 1.0, 0.0, 0.0], [1.0, 10.0])
    h_remez = firpm.remez(numtaps, [0.0, 0.2, 0.25, 0.5], [1.0, 0.0],
                          weight=[1.0, 10.0])
    assert np.array_equal(h_firpm, h_remez)


@pytest.mark.parametrize("numtaps", [31, 32])
def test_firpm_ftype_matches_remez_type(numtaps):
    h_firpm = firpm.firpm(numtaps - 1, [0.1, 0.9], [1.0, 1.0], "hilbert")
    h_remez = firpm.remez(numtaps, [0.05, 0.45], [1.0], type="hilbert")
    assert np.array_equal(h_firpm, h_remez)
    # A remez differentiator's slope is desired / fs per unit frequency;
    # firpm's amplitudes are given at the edges.
    h_firpm = firpm.firpm(numtaps - 1, [0.0, 0.9], [0.0, 0.45],
                          "differentiator")
    h_remez = firpm.remez(numtaps, [0.0, 0.45], [1.0],
                          type="differentiator")
    assert np.allclose(h_firpm, h_remez, atol=1e-12)


def test_firpm_ramp():
    # Amplitudes are linear between a band's edges.
    h = firpm.firpm(80, [0.0, 0.3, 0.4, 1.0], [1.0, 0.5, 0.0, 0.0])
    H = np.abs(np.fft.rfft(h, 4096))
    f = np.linspace(0, 1, len(H))
    m = f <= 0.3
    assert np.max(np.abs(H[m] - (1.0 - f[m] / 0.6))) < 1e-2


def test_firpm_bands_normalized_frequency():
    # firpm's Bands, and the frequencies its callables see, are normalized
    # to Nyquist = 1.
    seen = []
    def amplitude(f):
        seen.append(f)
        return 1.0
    h = firpm.firpm(50, [(0.0, 0.4, amplitude), (0.5, 1.0, 0.0, 10.0)])
    assert 0.0 <= min(seen) and max(seen) <= 0.4 + 1e-12
    h_remez = firpm.remez(51, [(0.0, 0.2, lambda f: 1.0),
                               (0.25, 0.5, 0.0, 10.0)])
    assert np.allclose(h, h_remez, atol=1e-12)


def test_firpm_bad_arguments():
    with pytest.raises(ValueError):
        firpm.firpm(20, [0.0, 0.4, 0.5, 1.0], [1.0, 1.0, 0.0])
    with pytest.raises(ValueError):
        firpm.firpm(20, [0.0, 0.4, 0.5, 1.0], [1.0, 1.0, 0.0, 0.0], [1.0])
    with pytest.raises(ValueError):
        firpm.firpm(20, [0.1, 0.9], [1.0, 1.0], "lowpass")
    with pytest.raises(TypeError):
        firpm.firpm(20, [0.1, 0.9], [1.0, 1.0], "hilbert", "hilbert")
    with pytest.raises(ValueError):
        firpm.firpm(20, [(0.1, 0.9, 1.0)], ftype="hilbert")
    with pytest.raises(ValueError):
        firpm.firpm(-1, [0.0, 1.0], [1.0, 1.0])


def test_no_module_level_strategy_aliases():
    assert not hasattr(firpm, "SCALING")
    assert not hasattr(firpm._firpm, "SCALING")
    assert firpm.Strategy.SCALING.name == "SCALING"


@pytest.mark.parametrize("numtaps", [31, 32])
@pytest.mark.parametrize("type, bands, desired, weight", [
    # A wideband differentiator, and one with a stopband.
    ("differentiator", [0.0, 0.45], [1.0], [1.0]),
    ("differentiator", [0.0, 0.2, 0.3, 0.5], [1.0, 0.0], [1.0, 1.0]),
    ("hilbert", [0.05, 0.45], [1.0], [1.0]),
])
def test_type_matches_scipy_remez(numtaps, type, bands, desired, weight):
    signal = pytest.importorskip("scipy.signal")
    if type == "hilbert" and numtaps % 2:
        # Type III is zero at Nyquist as well as at 0.
        bands = [0.05, 0.4]
    h = firpm.remez(numtaps, bands, desired, weight=weight, type=type)
    # remez's default grid is coarse enough to differ by about 5e-4.
    h_scipy = signal.remez(numtaps, bands, desired, weight=weight, type=type,
                           grid_density=256)
    assert len(h) == numtaps
    assert np.allclose(h, -h[::-1], atol=1e-12)
    assert np.allclose(h, h_scipy, atol=1e-4)


def test_differentiator_response():
    # Amplitude desired * f / fs, here with fs in Hz.
    fs = 1000.0
    h = firpm.remez(32, [0.0, 400.0], [2.0], type="differentiator", fs=fs)
    H = np.fft.rfft(h, 4096)
    f = np.linspace(0, fs / 2, len(H))
    m = f <= 400.0
    assert np.max(np.abs(np.abs(H[m]) - 2.0 * f[m] / fs)) < 1e-2


def test_type_needs_band_edges():
    with pytest.raises(ValueError):
        firpm.remez(31, [(0.05, 0.45, 1.0)], type="hilbert")
    with pytest.raises(ValueError):
        firpm.remez(31, [0.05, 0.45], [1.0], type="lowpass")


def test_fs_scales_band_edges():
    h = firpm.remez(51, *FLAT, weight=FLAT_WEIGHT)
    h_hz = firpm.remez(51, [0.0, 9600.0, 12000.0, 24000.0], FLAT[1],
                       weight=FLAT_WEIGHT, fs=48000.0)
    assert np.allclose(h, h_hz, atol=1e-12)


def test_default_weight_is_one():
    assert np.array_equal(firpm.remez(51, *FLAT),
                          firpm.remez(51, *FLAT, weight=[1.0, 1.0]))


@pytest.mark.parametrize("numtaps", [51, 52])
def test_numeric_bands_match_flat(numtaps):
    # Numeric bands are the flat specification, so the design is identical.
    h_flat = firpm.remez(numtaps, *FLAT, weight=FLAT_WEIGHT)
    h_bands = firpm.remez(numtaps, [firpm.Band(0.0, 0.2, 1.0),
                                    (0.25, 0.5, 0.0, 10.0)])
    assert np.array_equal(h_bands, h_flat)


@pytest.mark.parametrize("numtaps", [51, 52])
def test_callable_bands_match_flat(numtaps):
    # Constant callables must reproduce the flat specification, for both
    # type I (odd numtaps) and type II (even numtaps).
    h_flat = firpm.remez(numtaps, *FLAT, weight=FLAT_WEIGHT)
    h_bands = firpm.remez(numtaps, [(0.0, 0.2, lambda f: 1.0, lambda f: 1.0),
                                    (0.25, 0.5, lambda f: 0.0, lambda f: 10.0)])
    assert len(h_bands) == numtaps
    assert np.allclose(h_flat, h_bands, atol=1e-12)


def test_mixed_bands_match_flat():
    # In a design with a callable, numeric amplitudes and weights are
    # evaluated natively and must still match the flat form.
    h_flat = firpm.remez(61, *FLAT, weight=FLAT_WEIGHT)
    h_bands = firpm.remez(61, [(0.0, 0.2, 1.0),
                               (0.25, 0.5, lambda f: 0.0, 10.0)])
    assert np.allclose(h_flat, h_bands, atol=1e-12)


def test_ramp_callable():
    # A ramp is a callable; check the passband follows it.
    def ramp(f):
        return 1.0 - f / 0.3
    h = firpm.remez(81, [(0.0, 0.15, ramp), (0.2, 0.5, 0.0)])
    H = np.abs(np.fft.rfft(h, 4096))
    f = np.linspace(0, 0.5, len(H))
    passband = f <= 0.15
    assert np.max(np.abs(H[passband] - ramp(f[passband]))) < 1e-2


@pytest.mark.parametrize("fs", [None, 48000.0])
def test_callables_see_band_frequency(fs):
    # Callables receive frequency in the units of the band edges.
    scale = 1.0 if fs is None else fs
    edges = [(0.0, 0.2 * scale), (0.25 * scale, 0.5 * scale)]
    seen = {0: [], 1: []}
    def amplitude(band, value):
        def fn(f):
            seen[band].append(f)
            return value
        return fn
    firpm.remez(51, [(*edges[0], amplitude(0, 1.0)),
                     (*edges[1], amplitude(1, 0.0))], fs=fs)
    for band, (f0, f1) in enumerate(edges):
        assert seen[band]
        tol = 1e-12 * scale
        assert f0 - tol <= min(seen[band]) and max(seen[band]) <= f1 + tol


def test_bands_cic_droop_correction():
    # The motivating use case: pre-compensate the passband droop of a
    # third-order CIC decimator (R=8) with a non-flat amplitude target.
    R = 8
    def cic_inv(f):
        om = 2 * np.pi * f
        return 1.0 if om == 0 else (np.sin(om/2) / (R*np.sin(om/(2*R))))**-3
    numtaps = 202
    h = firpm.remez(numtaps, [(0.0, 0.04, cic_inv),
                              (0.06, 0.5, 0.0)])
    assert len(h) == numtaps
    H = np.abs(np.fft.fft(h, 8192))[:4097]
    f = np.linspace(0, 0.5, 4097)
    want = np.array([cic_inv(fi) for fi in f])
    passband = f < 0.04
    assert np.max(np.abs(H[passband] - want[passband])) < 1e-3
    assert np.max(H[f > 0.06]) < 1e-3


def test_bands_callback_exception_propagates():
    # Raised immediately (serial code path) and after many evaluations (from
    # inside the OpenMP worker threads of exchange()); both must surface as the
    # original Python exception rather than aborting or being swallowed.
    def boom(f):
        raise ValueError("boom")
    with pytest.raises(ValueError, match="boom"):
        firpm.remez(51, [(0.0, 0.2, boom), (0.25, 0.5, 0.0)])

    # A well-posed design (about 2700 calls to converge), so the exception
    # is raised mid-run whichever path the exchange takes.
    calls = [0]
    def late_boom(f):
        calls[0] += 1
        if calls[0] > 1000:
            raise KeyError("late")
        return 1.0
    with pytest.raises(KeyError):
        firpm.remez(201, [(0.0, 0.2, late_boom, lambda f: 1.0),
                          (0.25, 0.5, lambda f: 0.0, lambda f: 1.0)])


def test_bands_invalid_specification():
    # Unsorted bands, edges above Nyquist and nonpositive weights are
    # rejected rather than producing a filter with STATUS_SUCCESS.
    with pytest.raises(ValueError):
        firpm.remez(21, [(0.25, 0.5, lambda f: 0.0),
                         (0.0, 0.2, lambda f: 1.0)])
    with pytest.raises(ValueError):
        firpm.remez(21, [(0.0, 0.6, lambda f: 1.0)])
    with pytest.raises(ValueError):
        firpm.remez(21, [(0.0, 0.2, 1.0, -1.0), (0.25, 0.5, lambda f: 0.0)])


@pytest.mark.parametrize("bands, desired, weight", [
    ([0.25, 0.5, 0.0, 0.2], [0.0, 1.0], [1.0, 1.0]),   # unsorted edges
    ([0.0, 0.6], [1.0], [1.0]),                        # edge above Nyquist
    ([0.0, 0.2, 0.25, 0.5], [1.0, 0.0], [1.0, -1.0]),  # negative weight
])
def test_flat_invalid_specification(bands, desired, weight):
    # These used to come back as an empty tap vector.
    with pytest.raises(ValueError):
        firpm.remez(21, bands, desired, weight=weight)


def test_bad_arguments():
    with pytest.raises(ValueError):
        firpm.remez(21, FLAT[0], [1.0, 0.0, 0.0])
    with pytest.raises(ValueError):
        firpm.remez(21, FLAT[0], FLAT[1], weight=[1.0])
    with pytest.raises(ValueError):
        firpm.remez(21, [0.0, 0.2, 0.25], FLAT[1])
    with pytest.raises(TypeError):
        firpm.remez(21, [(0.0, 0.2, 1.0)], weight=[1.0])
    with pytest.raises(ValueError):
        firpm.remez(21, *FLAT, precision="single")
    with pytest.raises(ValueError):
        firpm.remez(0, *FLAT)
    with pytest.raises(ValueError, match="Sampling.*single scalar"):
        firpm.remez(11, [0.1, 0.4], [1.0], fs=np.array([10, 20]))


HAVE_MPFR = hasattr(firpm._firpm, "firpm_mp")


@pytest.mark.skipif(not HAVE_MPFR, reason="built without MPFR")
@pytest.mark.parametrize("numtaps", [51, 52])
def test_mpfr_matches_double(numtaps):
    callable_bands = [(0.0, 0.2, lambda f: 1.0),
                      (0.25, 0.5, lambda f: 0.0, 10.0)]
    for args, kwargs in [(FLAT, dict(weight=FLAT_WEIGHT)),
                         ((callable_bands,), {})]:
        h = firpm.remez(numtaps, *args, **kwargs)
        h_mp = firpm.remez(numtaps, *args, **kwargs, precision="mpfr")
        assert len(h_mp) == numtaps
        assert np.allclose(h, h_mp, atol=1e-10)


@pytest.mark.skipif(HAVE_MPFR, reason="built with MPFR")
def test_mpfr_unavailable():
    with pytest.raises(ValueError, match="MPFR"):
        firpm.remez(21, *FLAT, precision="mpfr")


def test_bands_even_numtaps_nonzero_nyquist_promotes_to_type_i():
    # Even numtaps with nonzero gain at Nyquist cannot be a type II filter;
    # the band form must promote to type I exactly as the flat form does.
    h_flat = firpm.remez(32, [0.0, 0.2, 0.25, 0.5], [0.0, 1.0])
    h_bands = firpm.remez(32, [(0.0, 0.2, lambda f: 0.0),
                               (0.25, 0.5, lambda f: 1.0)])
    assert len(h_flat) == len(h_bands) == 33
    assert np.allclose(h_flat, h_bands, atol=1e-12)


def test_narrow_band_reference_allocation():
    # Four bands whose widths round to 6 + 6 + 4 initial references used to
    # leave the narrow last band with none and corrupt the reference set.
    f = [0.0, 0.23, 0.3, 0.53, 0.6, 0.74, 0.8, 0.8000001]
    a = [1, 1, 0, 0, 1, 1, 0, 0]
    w = [1.0, 1.0, 1.0, 1.0]
    h = firpm.firpm(28, f, a, w)
    assert len(h) == 29 and np.all(np.isfinite(h))
    H = np.abs(np.fft.fft(h, 8192))[:4097]
    fr = np.linspace(0, 1, 4097)
    for i in range(4):
        m = (fr >= f[2*i]) & (fr <= f[2*i+1])
        if m.any():
            assert np.max(np.abs(H[m] - a[2*i])) < 0.1

