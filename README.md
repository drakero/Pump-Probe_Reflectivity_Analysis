## Synopsis
<p>
Code used for analyzing reflectivity data taken in pump-probe experiments. In these experiments, a "pump" laser pulse was focused 
onto a semiconductor such as GaAs. After some delay (femtoseconds to picoseconds), the "probe" laser pulse was focused onto the
same spot, interacting with region that had been excited by the pump pulse. The probe reflects off the surface and onto a
photodiode, giving a measurement of the surface reflectivity when paired with a measurement of the probe pulse energy before
reflection. This is repeated many times for each delay and for many delays, giving the time evolution of the surface's reflectivity
after excitation from the pump pulse.
</p>

<p> This code reads in the data, adjusts for drift in the photodiode signal, bins the data based on energy, and produces plots of
the reflectivity as a function of time and pump energy. Some rudimentary code for calculated the change in refractive index from
the Fresnel equations was also added in the Jupyter notebook. This value is then used to estimate the electron density in the
conduction band by using the Drude model. This addition is still a work in progress.</p>

## License
MIT license
