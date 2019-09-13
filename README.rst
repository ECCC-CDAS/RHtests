RHtests
========

Data homogenization
-------------------

*  This RHtestsV4 software package can be used to detect, and adjust for, multiple changepoints (shifts) that could exist in a data series that may have first order autoregressive errors but excluding daily precipitation data series. It is based on the penalized maximal t test (Wang et al. 2007) and the penalized maximal F test (Wang 2008b), which are embedded in a recursive testing algorithm (Wang 2008a), with the lag-1 autocorrelation (if any) of the time series being empirically accounted for. The problem of uneven distribution of false alarm rate and detection power is also greatly alleviated by using empirical penalty functions (Wang et al. 2007, Wang 2008b). The time series being tested may have zero-trend or a linear trend throughout the whole period of record. A homogenous time series that is well correlated with the base series may be used as a reference series. However, detection of changepoints is also possible with the RHtestsV3 package when a homogenous reference series is not available.

*  The RHtestsV4 is the RHtestsV3 with the addition of provision of QM-adjustments that are estumated with the use of a reference series. The RHtestsV3 package is an extended version of the RHtestV2 package. The extension includes: (1) provision of Quantile-Matching (QM) adjustments ( Wang et al. 2010, section 5 ) in addition to the mean-adjustments that were provided in the RHtestV2; (2) choice of the segment to which the base series is to be adjusted (referred to as the base segment); (3) choices of the nominal level of confidence at which to conduct the test; and (4) all functions are now available in the GUI mode.  It also conducts simple quality control on the input daily data.

Links
-----

* `Wang et al. 2010, section 5`_
* `R statistical programming language`_
* `RHtestsV3 Example One Data`_

.. _Wang et al. 2010, section 5: http://etccdi.pacificclimate.org/RHtest/transformTPRs.pdf
.. _R statistical programming language: http://www.r-project.org/
.. _RHtestsV3 Example One Data: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/Example1.dat

Using RHtests
==============

Running RHtests
----------------

* Run the following command to start ::

    > source('RHtestsV4_20180301.r')

Issues
------

* Please check the `issue page`_ and check if the issue is already reported and its current status.
* If the issue is not reported yet, please kindly submit a `new issue`_, tag the issue as bug and leave it unassigned. Please describe your issue in as much detail as possible and to include your output.

.. _issue page: https://github.com/ECCC-CDAS/RHtests/issues
.. _new issue: https://github.com/ECCC-CDAS/RHtests/issues/new

Contact Us
----------

* I am Rodney Chan and together with Yang Feng at `Climate Data and Analysis Section`_ of Environment and Climate Change Canada are the current maintainer of RHtests. You can contact me at rodney.chan@canada.ca or Yang at yang.feng@canada.ca

.. _Climate Data and Analysis Section: https://github.com/ECCC-CDAS
