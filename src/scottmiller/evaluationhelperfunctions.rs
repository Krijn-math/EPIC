//helper functions for eval

use num_bigint::BigUint;

use crate::{bigint::*, bigconstants::*};
use super::constructionhelperfunctions::*;



pub fn apply(
    step: Line,
    rx: &BigUint,
    ry: &BigUint,
) -> Result {
    //receives a tuple (lx, ly, l0)
    //computes the multiplication
    let mut re = fp_mul(&step.lx, &rx);
    let im = fp_mul(&step.ly, &ry);

    re = fp_sub(&step.l0, &re);

    Result { re, im }
}


pub fn square(
    f: Result
) -> Result {
    let (re, im) = fp2sqr(&f.re, &f.im);
    Result{ re, im }
}

pub fn absorb(
    f: Result,
    eval: Result
) -> Result {
    let (re, im) = fp2mul(&f.re, &f.im, &eval.re, &eval.im);
    Result{ re, im }
}