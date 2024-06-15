# Assorted header-only C++ libraries for Bioconductor

This package vendors an assortment of header-only C++ libraries for use in Bioconductor packages. 
The use of a central repository avoids duplicate vendoring of libraries across multiple R packages,
and enables better coordination of version updates across cohorts of interdependent C++ libraries.
This package is very minimalistic by design, to ensure that downstream packages are not burdened with more transitive dependencies.

All vendored libraries are stored in [`inst/include`](inst/include) and can be used via the usual `LinkingTo` mechanism.
The corresponding license for each library is stored in [`inst/licenses`](inst/include).
Developers should edit and run [`inst/fetch.sh`](inst/fetch.sh) to update the version of any library.
