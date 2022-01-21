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
