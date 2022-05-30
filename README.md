# Paraxial-Optics

A MATLAB toolkit for paraxial (Gaussian) optics and refractive sequential ray-tracing with matrix methods in the meridian plane. Includes arbitrary number of optical surfaces (lenses) and apertures, and calculates cardinal points, image formation, and entrance/exit pupils & windows. 

Mostly meant for education purposes, e.g., for undergraduate physics or engineering courses.

## Brief description

Paraxial or Gaussian optics deal with ray-tracing within the geometric optics (high-frequency) approximation, where rays are refracted (or reflected) by shapes (called optical surfaces) with dimensions much greater than the wavelength. Optical surfaces are defined by their optical power (or curvature and refractive index), i.e., how much they "bend" the light rays hitting on them, and their aperture size, i.e., how big (in diameter) a lens is. Clear apertures (those of zero optical power) can be used, for more complex compound systems, e.g., to define a variable F-number iris or an image plane, controlling depth of focus or the sensor-size in a camera. 

The present 2D implementation considers only rays within the meridional plane, i.e., the plane defined by the optical axis and an off-axis point-object. Also, frequency and polarization are irrelevant here, but one can implement chromatic dispersion (via wavelength-dependent optical power) and/or Fresnel transmission coefficients assuming a given linear polarization. Finally, note that this presently supports only sequential ray-tracing in refractive (dioptric) optical systems, i.e., composed of glass lens; development for reflective (catoptric) or catadioptric systems is straightforward, but non-sequential ray tracing and/or 3D is much more complicated.

## Examples

Consider the following includes example scripts to understand the syntax, use and visualization offered by this code-pack. I am also giving the figures that MATLAB should produce.

### 0.myScript0_OneLens_OneRay.m
   - See how an OS and a ray are defined.
   - See how how ray-tracing is done and plotted
   - Change the parameters of the OS or the ray, and see what changes
   - Add mores lenses (rows in "OS" matrix) or rays (rows in "rays" matrix), and see what change

![myScript0_OneLens_OneRay_MATLAB](https://user-images.githubusercontent.com/97299585/171019676-a416c69c-30e2-48cf-99a1-19178da415fa.png)
Fig. 0: One lens system, with a single ray. Marked also are the focal points (F1/2), principal planes (H1/2), and nodal points (N1/2).

### 1.myScript1_Magnifier.m
   - Describes how to define an optical system (as ABCD or as "optical surface" [OS] matrix)
   - Defines a couple of test-case optical systems and point objects
     * TestCase = 1 is a diverging thin-lens and a virtual object.
     * TestCase = 2 is the converging ABCD system used as a magnifier
   - Calls top-level image-formation/visualization function "ImageFormation_with_CardinalRayTracing"
   - Observe the cardinal points and three cardinal rays traced to form a real (or virtual object)
   - Change the parameters of the system (e.g. the optical power C=-P, converging or diverging) and/or position of object, and observe how it affects image formation

![myScript1_Magnifier_MATLAB](https://user-images.githubusercontent.com/97299585/171019703-02329528-dd90-41b8-abdb-1659a35615ac.png)
Fig. 1: Three main rays in two more complicated optical systems. In the left, a diverging lens and a virtual object. In the right, a magnifier (a close-by real object producing an erect magnified virtual image)
   
### 2.myScript2_TelescopeMicroscope.m
   - Defines "custom" ray tracing for a Keplerian astronomical Telescope and a Compound Microscope
   - Produces two figures for the corresponding systems   - 

   - See how custom ray-sets (bundles) can be defined and simulated
   - Microscope: Calculate microscope magnification; check how object-plane position (objS) 
                 affects image-plane position (we want VIRTUAL image!) and magnification
   - Telescope : Measure width and angle of outgoing ray-bundles, and compare them to the 
                 corresponding incoming ray-bundle parameters (width, angle); how does this compare
                 to the telescope magnification?
                 
![myScript2_TelescopeMicroscope_Octave](https://user-images.githubusercontent.com/97299585/171019783-706dcbd4-2839-4573-bc7a-d1d8531a614a.png)
Fig. 2: Two two-lens systems, with multiple rays traced. In the left, a Kepler telescope, with on-axis (red) and off-axis (cyan) pencil beams. In the right, a compound microscope (magnified inverted virtual image).
     
### 3.myScript3_Apertures.m
   - Produces three figures
     * One test figure with image-formation rays, 
     * A figure with test rays for Aperture-Stop (AS) and Pupils calculations 
     * A figure with test rays for Field-Stop (FS) and Windows calculations 
   - Read through the script comments/guidelines to understand how these are defined and simulated
  
![myScript3_Apertures_MATLAB](https://user-images.githubusercontent.com/97299585/171019888-459ce761-a4a0-478c-9853-f33fd895d24a.png)
Fig. 3: Study of a complex optical system with mulitple surfaces/lenses and two apertures. In the left, image formation from an object in infinity (black pencil beam) or a close-by point object. In the middle, a diverging cone of rays from an on-axis object to identify the aperture-stop and the pupils of the system. In the right, a converging cone of rays to identify the windows and field-stop of the system.

## List of all functions (core and utility) included:
* Convert_ABCD2OS.m (converts ABCD to optical surface matrix)
* DoParaxialAnalysis.m (calculates ABCD from optical surfaces, cardinals points etc)
* DoParaxialImageFormation.m (calculates position & magnification of image)
* DoParaxialRayTracing.m (traces given set of rays from entrance to exit)
* ImageFormation_with_CardinalRayTracing.m (top-level image-formation/visualization function)
* PlotParaxialOS.m (plots the optical system)
* calculateVisibleSpectrumColor.m (colors from wavelength)
