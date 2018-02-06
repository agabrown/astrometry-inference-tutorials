# CBJ November 2017 (Coryn Bailer-Jones, calj@mpia.de)
# Various mathematical and astronomical constants

# mathematical constants
conv <- pi/180 # degrees -> radians
km <- 1e3      # km in SI (m)

# physical constants
#G <- 6.67428e-11 # gravitational constant

# solar system constants
sgp <- 1.32712440018e20   # standard gravitational constant in SI (G*Msun)
yr  <- 365.25*24*3600     # length of Julian year in SI (s). Tropical year is 365.242190 days
au  <- 1.495978701e11     # AU in SI (m) (Hipparcos catalogue documentation p. 25, Table 1.2.2)
pc  <- (3600*180/pi)*au   # parsec in SI (m) [3.085678e+16]
kf  <- au/(yr*km)         # 1 AU/yr in km/s. Units: km s^-1 yr [4.740471]
                          # converts proper motion of 1as/yr at 1pc to velocity in km/s 
                          # year unit here must be same as in PMs

# Galactic coordinate system constants (l is Galactic longitude)
# theta0, raNGP, decNGP uniquely define the equatorial to Galactic coordinate transformation.
# With the Hipparcos values for these, the matrix trans in eq.to.gal() is equal the transpose 
# of the matrix A_G in the Hipparcos catalogue documentation p. 92, eqn. 1.5.11) to 10 dps.]
theta0 <- 90+32.93192  # PA of NGP = l of north celestial pole = 90 + l_omega 
                       # l_omega = l of ascending node of Galactic plane on the equator of ICRS
                       # (Hipparcos catalogue documentation p. 91, eqn. 1.5.10) (theta0=123 in B1950; Johnson & Soderblom 1987) 
raNGP  <- 192.85948 # (Hipparcos catalogue documentation p. 91, eqn. 1.5.9)
decNGP <- 27.12825  # (Hipparcos catalogue documentation p. 91, eqn. 1.5.9)
