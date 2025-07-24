# Introduction to TSMP-PDAF

TSMP-PDAF describes the build commands of
[TSMP](https://github.com/HPSCTerrSys/TSMP2/) that can introduce data
assimilation for an ensemble of TSMP simulations using the Parallel
Data Assimilation Framework ([PDAF](http://pdaf.awi.de/trac/wiki)).

Note that TSMP-PDAF does not follow the multiple program multiple data
(MPMD) paradigm that TSMP does. Instead component models are loaded as
libraries and compiled into a single executable (@Kurtz2016).

## Issues

If you encounter errors or issues regarding TSMP-PDAF, please add them
in one of the following places:

- TSMP-PDAF issues (related to source code)
  https://github.com/HPSCTerrSys/pdaf/issues
- Issues related to building TSMP-PDAF using TSMP2
  https://github.com/HPSCTerrSys/TSMP2/issues
- TSMP-PDAF issues (internal)
  https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/tsmp-pdaf-issues


## Citing TSMP-PDAF

Please cite the following when using TSMP-PDAF in a publication

* Kurtz, W., He, G., Kollet, S. J., Maxwell, R. M., Vereecken, H., &
  Hendricks Franssen, H. J. (2016). TerrSysMP–PDAF (version 1.0): a
  modular high-performance data assimilation framework for an
  integrated land surface–subsurface model. Geoscientific Model
  Development, 9(4), 1341-1360. doi:
  [10.5194/gmd-9-1341-2016](http://dx.doi.org/10.5194/gmd-9-1341-2016)

* Nerger, L., & Hiller, W. (2013). Software for ensemble-based data
  assimilation systems - Implementation strategies and
  scalability. Computers & Geosciences, 55(), 110–118. doi:
  [j.cageo.2012.03.026](http://dx.doi.org/10.1016/j.cageo.2012.03.026)

## Contributors

In alphabetic order (to be extended):

* Yorck Ewerdwalbesloh
* Guowei He
* Johannes Keller
* Wolfgang Kurtz
* Stefan Poll
* Mukund Pondkule
* Prabhakar Shrestha
* Lukas Strebel
