This repository is associated with the paper "Effective Pairings in Isogeny-based Cryptography" (also known as EPIC).
The algorithms described in this work can be found in the different subfolders of `/src`

This work is in no way production-ready. Instead, it should be treated as a proof-of-concept for using pairings in isogeny-based cryptography,
and may contain bugs and/or typo's.

EDIT (November, 2024).
Using biextensions, Robert (2024/517) has made pairing computations on Montgomery curves much easier. One should use such pairings
in combinations with the algorithms given above, which presumably gives better results. A preliminary implementation of such pairings
is straightforward, I have made one in Magma which can be requested by email (krijn (at) q (one) q (one) (dot) nl).
