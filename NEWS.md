# chromePlus 2.0.0

* Added `plotChromeplus()` for density plots of MCMC results with HPD intervals.
* Added `simChrom()` for simulating chromosome number evolution under multiple models.
* Added `makeSSEchrom()` for simulating chromosome data under binary SSE models.
* Added `startVals()` helper for generating starting values from uniform or normal distributions.
* Added sex chromosome-autosome fusion (SAF) model support in `constrainMkn()`.
* Added `sym.hyperstates` and `oneway` constraint options to `constrainMkn()`.
* Improved vignette with detailed model descriptions and worked examples.

# chromePlus 0.0.1

* Initial release with `constrainMkn()`, `constrainMuSSE()`, and `datatoMatrix()`.
* Implements chromEvol-style models using the diversitree framework.
