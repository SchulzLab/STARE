============
Affinity summarisation
============

To get TF scores for a gene, STARE summarises the TF-affinities from regions that were assigned to that gene. Since STARE offers different approaches to get the region-gene mapping, it also has different ways to summarise the TF affinities. 

The basic idea is always the same. We sum up the *tf* affinities of the regions *r* that were linked to a gene *g*:

.. image:: /Figures/STARE_BaseEquation.png
  :alt: STARE_BaseEquation

Now the *scaler* is where it changes depending on how the region-gene mapping was done. The different approaches make use of different types of data, which allows to integrate more data modalities.

Gene window
******************

In case that all regions within a window around the TSS are considered relevant, the equation would change to the following: 

.. image:: /Figures/STARE_WindowEquation.png
  :alt: STARE_WindowEquation

The scaling with the activity of the region is only done when activity columns were specific (-n). The scaling with the exponential decay function can be turned off with *-e False*, for example for very small window sizes. The distance decay can only be turned off for the Gene window-version, for the ABC-mode it is always on.

ABC-scoring
******************

If you did the region-gene mapping with the ABC-score, we change the scaling and also distinguish between regions close to the promoter and more distant ones. The reason is that the measurement of contact frequency of a region with itself is likely flawed. Note that the -e flag for the distance decay does not change the scaling here.

.. image:: /Figures/STARE_AdaptedABCEquation.png
  :alt: STARE_AdaptedABCEquation

Not-adapted ABC-scoring
******************

If you chose the run the ABC-scoring without the adapted score, we can't use the adapted activity for scaling. Thus, the affinity scaling changes to:

.. image:: /Figures/STARE_UnadaptedABCEquation.png
  :alt: STARE_UnadaptedABCEquation

The -e flag also doesn't change the scaling in this setup.



