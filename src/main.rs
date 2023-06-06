#![allow(non_snake_case)]
#![allow(dead_code)]

mod bigconstants;
mod biggenverif;
mod bigint;

mod multimiller {
    pub mod bigmultimillerhelper;
    pub mod bigmultimiller;
    pub mod multi_window;
}

mod nafmiller {
    pub mod bignafmiller;
}

mod supersingularity {
    pub mod doliskani;
    pub mod supersingularity1;
    pub mod supersingularity2;
}

mod tradmiller {
    pub mod bigmiller;
    pub mod bigmillerhelper;
}

mod windowmiller {
    pub mod bigger_window;
    pub mod bigwindowmiller;
}


mod bigfastfinding;
use std::env;
mod elligator;
mod precomputed;

use biggenverif::check_if_primitive;
use bigint::*;
use num_bigint::{BigUint};

mod bigcurve;

mod bigprodtree;

mod bignaive;

mod biglucaspairing;
use biglucaspairing::*;

use crate::{bigconstants::*, bigprodtree::*};
use crate::windowmiller::bigwindowmiller::*;

pub fn fast_verif_with_gen(tpx: &BigUint, tpy: &BigUint, tmx: &BigUint, tmy: &BigUint, mont_A: &BigUint, log: &BigUint) -> bool {

    let r: BigUint = (*CSIDH512).clone() + &(*BIGONE);
    let f = windowmiller(&r, tpx, tpy, tmx, tmy, mont_A);
    print_big_cost();

    let (asq, bsq) = (fp_sq(&f.re), fp_sq(&f.im));
    let abinv = non_constant_inversion(&fp_add(&asq, &bsq));


    let zetare = fp_mul(&fp_sub(&asq, &bsq), &abinv);
    print_big_cost();

    check_if_primitive(&zetare, log)
}


fn main() {
    env::set_var("RUST_BACKTRACE", "1");
    println!("Hello, curves!");


    //our torsion basis for E3 over CSIDH-512
    let mont_A: BigUint = BigUint::parse_bytes(b"53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340", 16).unwrap();
    let px: BigUint = fp_neg(&BigUint::parse_bytes(b"17", 10).unwrap());
    let qx: BigUint = BigUint::parse_bytes(b"16", 10).unwrap();
    let py: BigUint = BigUint::parse_bytes(b"B2EF33A6FACF4ED15A1AC03EA6C9473E08DAAC754D9CCC272C775D3B89A4D42DA1B8BB396430D3198FC61F45D5927F183424A5848FB23A374AB8C57B72E19F", 16).unwrap();
    let qy: BigUint = BigUint::parse_bytes(b"28302633848ABDDFEC053E481685A2F4BF36C4405D9A5FB226160A1076AB579EC5C08E118E02E4DB1D5DF0F5C2E5DAEA84F3690086295C19A9A4EC3DE615A505", 16).unwrap();


    assert!(lucaspairing_algorithm(&px, &py, &qx, &qy, &mont_A));
    reset_big_cost();
    assert!(product_tree_algorithm(&px, &mont_A));
    assert!(product_tree_algorithm(&qx, &mont_A));
    print_big_cost();
}
