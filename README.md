# Mueller matrix calculator, based on a JVIE-T-matrix code
The calculator allows one to compute the T-matrices for arbitrary geometries, constructed as tetgen-compatible tetrahedral meshes. 

## Usage
Key parameters are the particle size and incident EM-field (wavelength band, currently supports logarithmic scaling). Also, for the Mueller (scattering) matrix the angular resolution (N<sub>theta</sub>, N<sub>phi</sub>) for scattering angle and azimuthal angle must be given.

### Mueller matrix indexing
The resulting scattering (and extinction) matrices are written in a single file containing indices N<sub>size</sub>, N<sub>ia</sub>, and N<sub>pts</sub> per each row. A row also contains the Mueller matrix elements S<sub>11</sub>, S<sub>12</sub>, ..., S<sub>34</sub>, S<sub>44</sub>, and the scattering (or extinction) cross section c<sub>sca</sub> (c<sub>ext</sub>).

* N<sub>size</sub>: 1--bars (WILL run from the smallest size parameter to largest, i.e. opposite to the order of wavelenghts printed)
* N<sub>ia</sub>: 0--2 (0=randomly oriented case, will limit N<sub>pts</sub> to N<sub>theta</sub>, 1=incident angle 90 deg, 2=incident angle 180 deg)
* N<sub>pts</sub>: 1--N<sub>theta</sub>N<sub>phi</sub>. The current indices i<sub>theta</sub> and i<sub>phi</sub> for theta and phi, respectively, given by i<sub>theta</sub> = (N<sub>pts</sub>-1)%N<sub>theta</sub>, i<sub>phi</sub> = floor(N<sub>pts</sub>/N<sub>theta</sub>). The angles theta and phi are then given by theta = pi*(i<sub>theta</sub>+0.5)/N<sub>theta</sub>, phi = 2pi i<sub>phi</sub>/N<sub>phi</sub>. NOTE: i<sub>phi</sub> runs from 0 to N<sub>phi</sub>-1, i<sub>theta</sub> runs from 0 to N<sub>theta</sub>-1. 

