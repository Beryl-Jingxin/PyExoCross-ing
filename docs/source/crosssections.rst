Cross sections
==============

Please provide the line lists, ``Temperature``, ``Pressure``, 
wavenumber ``Range``, number of the points ``Npoints`` or bin size 
``BinSize`` and uncertainty filter ``UncFilter``.

``Npoints`` is the number of the points in grid.

``BinSize`` is the interval size of the grid.

``UncFilter`` means uncertainty filter. 
If you don't use the uncertainty filter, write ``N`` here. 
You don't need to change the number behind it.


*Example*

.. code:: bash

    # Calculate stick spectra or cross-sections #
    Temperature                             300
    Pressure                                1
    Range                                   0          30000
    Npoints/BinSize                         Npoints    30001

    Cutoff(Y/N)                             Y          25             # Default value 25
    Threshold(Y/N)                          Y          1e-30          # Default value 1e-30
    UncFilter(Y/N)                          Y          0.001
    QNsFilter(Y/N)                          Y          par[+]   e/f[e]   v[0,1,2,3]  
    DopplerHWHM(Y/N)                        Y          3              # Set Doppler HWHM as a constant
    LorentzianHWHM(Y/N)                     Y          0.7            # Set Lorentzian HWHM as a constant

    Broadeners                              Default    Air    Self    H2    He    CO2
    Ratios                                  1.0        0.0    0.0     0.0   0.0   0.0

    Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
    Profile                                 Gaussian        
    Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'

.. code:: bash

    # Calculate stick spectra or cross-sections #
    Temperature                             1000
    Pressure                                1.2
    Range                                   1000       5000
    Npoints/BinSize                         BinSize    0.1

    Cutoff(Y/N)                             N          25             # Default value 25
    Threshold(Y/N)                          N          1e-30          # Default value 1e-30
    UncFilter(Y/N)                          N          0.001
    QNsFilter(Y/N)                          N          par[+]   e/f[e]   v[0,1,2,3]  
    DopplerHWHM(Y/N)                        N          3              # Set Doppler HWHM as a constant
    LorentzianHWHM(Y/N)                     N          0.7            # Set Lorentzian HWHM as a constant

    Broadeners                              Default    Air    Self    H2    He    CO2
    Ratios                                  0.0        0.7    0.3     0.0   0.0   0.0

    Absorption/Emission                     Emission                # 'Absorption' or 'Emission'
    Profile                                 SciPyVoigt        
    Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'


**Note**

If the line lists data is not in the ExoMol format, 
please convert your data format into the ExoMol format at first 
and then compute partition functions with *PyExoCross*.
