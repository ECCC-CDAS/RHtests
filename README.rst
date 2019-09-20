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
* `RHtestsV4 last version`_
* `RHtestsV4 User Manual`_
* `RHtestsV4 User Manual - French Version`_
* `RHtests_dlyPrcp4 Code`_
* `RHtests_dlyPrcp User Manual`_
* `RHtestsV3 Example One Data`_
* `RHtestsV3 Example Two Base Data`_
* `RHtestsV3 Example Two Ref Data`_
* `Description of QM adjustment algorithm, Wang et al. 2010, section 5`_
* `Table of empirical percentiles of the PMT test statistic, included in the codes`_
* `Table of empirical percentiles of the PMF test statistic, included in the codes`_

.. _Wang et al. 2010, section 5: http://etccdi.pacificclimate.org/RHtest/transformTPRs.pdf
.. _R statistical programming language: http://www.r-project.org/
.. _RHtestsV4 last version: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/RHtestsV4_20190301.r
.. _RHtestsV4 User Manual: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/RHtestsV4_UserManual_10Dec2014.pdf
.. _RHtestsV4 User Manual - French Version: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/RHtestsV4_UserManual_10Dec2014_French.pdf
.. _RHtests_dlyPrcp4 Code: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/RHtests_dlyPrcp_20130719.r
.. _RHtests_dlyPrcp User Manual: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/RHtests_dlyPrcp_UserManual_10Dec2014.pdf
.. _RHtestsV3 Example One Data: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/Example1.dat
.. _RHtestsV3 Example Two Base Data: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/Example2.dat
.. _RHtestsV3 Example Two Ref Data: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/Example2_Ref.dat
.. _Description of QM adjustment algorithm, Wang et al. 2010, section 5: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/transformTPRs.pdf
.. _Table of empirical percentiles of the PMT test statistic, included in the codes: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/PTmaxRed_Nmin5_6CVs.txt
.. _Table of empirical percentiles of the PMF test statistic, included in the codes: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/PFmax31red_Nmin10_6CVs.txt

*  Users who wish to use the RClimDex software package (see below) to calculate climate indices, and who wish to detect and adjust for artificial shifts in the daily data series to be used as input to the RClimDex (which is recommended), should read the Quick Guide to RClimDex and RHtests Users, and should also download the following: 

* `Quick Guide to RClimDex and RHtests Users`_
* `Quick Guide to RClimDex and RHtests Users - French version`_
* `RClimDex-RHtests data format conversion software`_
* `Software for homogenization of daily precipitation data series, Wang et al. 2010`_
* `Example for daily precipitation data`_

.. _Quick Guide to RClimDex and RHtests Users: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/QuickGuide_to_RClimDex_and_RHtests.doc
.. _Quick Guide to RClimDex and RHtests Users - French version: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/QuickGuide_to_RClimDex_and_RHtests.French.doc
.. _RClimDex-RHtests data format conversion software: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/RClimDex_RHtest.r
.. _Software for homogenization of daily precipitation data series, Wang et al. 2010: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/RHtests_dlyPrcp.r
.. _Example for daily precipitation data: https://github.com/ECCC-CDAS/RHtests/blob/master/V4_files/RHtests_dlyPrcp_ExampleData.txt

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
