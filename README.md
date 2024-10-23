Match and join astronomical catalogs, replicating the functionality of TOPCAT.

Astronomy research often involves matching catalogs of sources. Matching allows one for example to associate sources 
detected at different wavelengths, allowing the construction of a spectral energy distribution.

TOPCAT is a commonly used piece of software for this purpose. However, its interactive/GUI nature limits its usefulness in 
building a data reduction/analysis pipeline.

This small tool wraps astropy.coordinates matching functions and applies them to catalogs, implementing 
additional match and join types. This can replicate the matching functionality of TOPCAT.

Supports matching using 3 methods:
* Finding the best match in catalog B for each entry in catalog A within the specified radius
* Finding the best *symmetric* match, meaning a match is only valid when element x in A is the 
closest match to element y in B, and vice versa.
* Finding all matches within the specified radius.

And supports joining using 3 methods:
* Joining catalog A **and** B where a match was found.
* Joining A and B, preserving entries in A which did **not** have a match in B.
* Selecting those entries in A which did **not** have a match in B.