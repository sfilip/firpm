import pytest
import pyfirpm

import numpy as np
import matplotlib.pyplot as plt


@pytest.mark.parametrize("n, fpass, fstop, ripple, attenuation", [
    # Specification 4
    (1000, 0.8, 0.81, 0.1, -90),
    (2001, 0.8, 0.81, 0.1, -90),

    # Specification 5
    (2002, 0.1, 0.105, 0.1, -90),
])
def test_lpf_spectrum(request, n, fpass, fstop, ripple, attenuation):
    success = True

    NFFT=1024

    h = pyfirpm.firpm(n, [0.0, fpass, fstop, 1.0], [1.0, 1.0, 0.0, 0.0], [1.0, 10.0])

    # Form spectrum in dB
    hdb = 20 * np.log10(np.abs(np.fft.fft(h)))
    f = np.linspace(0, 2, len(hdb))

    # Plot spectrum in dB
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(f, hdb)

    # Verify passband ripple
    passband_exceptions = np.argwhere((f < fpass) * (abs(hdb) > ripple) * (f < 1))
    if len(passband_exceptions):
        ax.plot(f[passband_exceptions], hdb[passband_exceptions],
                'ro', fillstyle='none')
        success = False

    # Verify stopband attenuation
    amp_stopband = np.argwhere((f > fstop) * (hdb > attenuation) * (f < 1))
    if len(amp_stopband):
        ax.plot(f[amp_stopband], hdb[amp_stopband],
                'ro', fillstyle='none')
        success = False

    # Produce plot
    ax.set_xlim([0, 1])
    ax.set_xlabel('Frequency')
    ax.set_ylabel('Amplitude (dB)')
    ax.set_title(request.node.name)
    ax.grid(True)
    fig.savefig(f'{request.node.name}.pdf')

    # Raise deferred exception if we failed
    assert success


@pytest.mark.parametrize("n", [50, 51])
def test_bands_matches_flat(n):
    # Constant per-band callables must reproduce the flat-vector API exactly,
    # for both type I (even n) and type II (odd n).
    f, a, w = [0.0, 0.4, 0.5, 1.0], [1.0, 1.0, 0.0, 0.0], [1.0, 10.0]
    h_flat = pyfirpm.firpm(n, f, a, w)
    h_bands = pyfirpm.firpm_bands(n, [(0.0, 0.4, lambda om: 1.0, lambda om: 1.0),
                                      (0.5, 1.0, lambda om: 0.0, lambda om: 10.0)])
    assert len(h_bands) == n + 1
    assert np.allclose(h_flat, h_bands, atol=1e-12)


def test_bands_cic_droop_correction():
    # The motivating use case: pre-compensate the passband droop of a
    # third-order CIC decimator (R=8) with a non-flat amplitude target.
    R = 8
    def cic_inv(om):
        return 1.0 if om == 0 else (np.sin(om/2) / (R*np.sin(om/(2*R))))**-3
    n = 201
    h = pyfirpm.firpm_bands(n, [(0.0, 0.08, cic_inv, lambda om: 1.0),
                                (0.12, 1.0, lambda om: 0.0, lambda om: 1.0)])
    assert len(h) == n + 1
    H = np.abs(np.fft.fft(h, 8192))[:4097]
    om = np.linspace(0, np.pi, 4097)
    want = np.array([cic_inv(o) for o in om])
    passband = om < 0.08*np.pi
    assert np.max(np.abs(H[passband] - want[passband])) < 1e-3
    assert np.max(H[om > 0.12*np.pi]) < 1e-3


def test_bands_callback_exception_propagates():
    # Raised immediately (serial code path) and after many evaluations (from
    # inside the OpenMP worker threads of exchange()); both must surface as the
    # original Python exception rather than aborting or being swallowed.
    def boom(om):
        raise ValueError("boom")
    with pytest.raises(ValueError, match="boom"):
        pyfirpm.firpm_bands(50, [(0.0, 0.4, boom, lambda om: 1.0),
                                 (0.5, 1.0, lambda om: 0.0, lambda om: 1.0)])

    calls = [0]
    def late_boom(om):
        calls[0] += 1
        if calls[0] > 5000:
            raise KeyError("late")
        return 1.0
    with pytest.raises(KeyError):
        pyfirpm.firpm_bands(400, [(0.0, 0.4, late_boom, lambda om: 1.0),
                                  (0.5, 1.0, lambda om: 0.0, lambda om: 1.0)])


def test_bands_invalid_specification():
    # Unsorted bands and edges outside [0, 1] are rejected rather than
    # producing a filter with STATUS_SUCCESS.
    with pytest.raises(RuntimeError):
        pyfirpm.firpm_bands(20, [(0.5, 1.0, lambda om: 0.0, lambda om: 1.0),
                                 (0.0, 0.4, lambda om: 1.0, lambda om: 1.0)])
    with pytest.raises(RuntimeError):
        pyfirpm.firpm_bands(20, [(0.0, 1.2, lambda om: 1.0, lambda om: 1.0)])


def test_bands_odd_n_nonzero_nyquist_promotes_to_type_i():
    # Odd n with nonzero gain at Nyquist cannot be a type II filter; the band
    # overload must promote to type I exactly as the flat-vector API does.
    f, a, w = [0.0, 0.4, 0.5, 1.0], [0.0, 0.0, 1.0, 1.0], [1.0, 1.0]
    h_flat = pyfirpm.firpm(31, f, a, w)
    h_bands = pyfirpm.firpm_bands(31, [(0.0, 0.4, lambda om: 0.0, lambda om: 1.0),
                                       (0.5, 1.0, lambda om: 1.0, lambda om: 1.0)])
    assert len(h_flat) == len(h_bands) == 33
    assert np.allclose(h_flat, h_bands, atol=1e-12)
