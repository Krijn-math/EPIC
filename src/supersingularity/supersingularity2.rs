//testing supersingularity by first computing pairing, then checking order is > 4sqrt(p)
//computes the pairing of degree sqrt(p) instead of degree p:
//may still many bugs, but in terms of cost estimation it should work

use num_bigint::{BigUint, ToBigUint, ToBigInt};
use num_traits::One;
use lazy_static::lazy_static;


use crate::{bigconstants::*, bigint::*, nafmiller::bignafmiller::*, bigcurve::*, elligator::*};

lazy_static! {
    pub static ref UUU: Vec<i128> = [ 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241,
    251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
    349, 353, 359, 367, 373, 587 ].to_vec();
    pub static ref NUUU: usize = UUU.len();

        
    pub static ref pUUU: BigUint = BigUint::parse_bytes(b"10A59877399FC365649BD0A396FF5CCF84EAAEB92E2E947EA2113B6E1D850F5F7C67", 16).unwrap();
    pub static ref cofacU: BigUint = BigUint::parse_bytes(b"61C0D66BDC739EEA6B062B4B8773657BA5C8DC4BFA21428704A6F45CE2624", 16).unwrap();

    pub static ref BOUND: BigUint = BigUint::parse_bytes(b"2856F1399D91D6592142B95420000000000000000000000000000000000000000", 16).unwrap();

    pub static ref NAFUUU: Vec<i8> = naf(&(*pUUU).clone().to_bigint().unwrap());

}

pub fn lucasPowOrderRecSQRT(
    are: &BigUint,
    low: usize,
    upp: usize,
    mut m: BigUint,
) -> BigUint {

    if upp - low == 0 {
        if is_one(are) {
            return m;
        }
        let k: BigUint = UUU[upp].to_biguint().unwrap();
        let zre = lucasonepow(are, &k);
        if zre {
            return &k * &m;
        }
        panic!(
            "THIS CANNOT LEGALLY HAPPEN upp = {}, low = {}, and a is {} ",
            upp, low, &are
        );
    }

    let mid: usize = mid(upp, low);
    let k: BigUint = ProdVec(&UUU, low, mid);
    let lre: BigUint;
    lre = lucaspow(are, &k);
    m = lucasPowOrderRecSQRT(&lre, mid + 1, upp, m);

    if BigUint::ge(&m, &*BOUND)
    {
        return m;
    }

    let k: BigUint = ProdVec(&UUU, mid + 1, upp);
    let rre: BigUint;
    rre = lucaspow(are, &k);

    lucasPowOrderRecSQRT(&rre, low, mid, m)
}

pub fn has_large_order(zetare: &BigUint) -> bool {
    let ord = lucasPowOrderRecSQRT(&zetare, 0, *NUUU - 1, BigUint::one());
    BigUint::ge(&ord, &*BOUND)
}


pub fn nafumiller(
    _n: &BigUint,
    px: &BigUint,
    py: &BigUint,
    rx: &BigUint,
    ry: &BigUint,
    mont_A: &BigUint,
) -> (BigUint, BigUint) {
    //cleanmiller but implemented using a NAF, using the fact that we can easily add or substract

    let mut tx2: BigUint = fp_sq(px);
    let mut txtz: BigUint = px.clone();
    let mut tz2: BigUint = (*BIGONE).clone();
    let mut tytz: BigUint = py.clone();

    //our to be return value, in Fp2 but can ignore denominator, its fp invariant
    let mut fxre: BigUint = (*BIGONE).clone();
    let mut fxim: BigUint = (*BIGZERO).clone();

    //we loop per bit over n
    for i in 1..NAFUUU.len() {
        ((tx2, txtz, tz2, tytz), (fxre, fxim)) =
            double(&tx2, &txtz, &tz2, &tytz, &fxre, &fxim, &rx, &ry, mont_A);

        //add or subtract depending on the bit
        if NAFUUU[i] == 1 {
            ((tx2, txtz, tz2, tytz), (fxre, fxim)) = add(
                &tx2, &txtz, &tz2, &tytz, &fxre, &fxim, &px, &py, &rx, &ry, mont_A,
            );
        } else if NAFUUU[i] == -1 {
            ((tx2, txtz, tz2, tytz), (fxre, fxim)) = subtract(
                &tx2, &txtz, &tz2, &tytz, &fxre, &fxim, &px, &py, &rx, &ry, mont_A,
            );
        }
        //println!("bit i: {} f is {}", i, to_aff(&fxre, &fxim));
    }

    (fxre, fxim)
}


pub fn is_supersingular(mont_A: &BigUint) -> bool{
    let (px, _py, qx, qy) = shitty_elligator(5, &mont_A); 
    print_big_cost();

    //we can adapt this to also keep track of py, so do arithmetic not just x-only
    let (tx, tz) = spe_ladder(&*cofacU, &px, mont_A);
    print_big_cost();


    let (tx, ty) = recover(&tx, &tz, &mont_A);

    print_big_cost();

    reset_big_cost();
    let (fxre, fxim) = nafumiller(&*pUUU, &tx, &ty, &qx, &qy, &mont_A);

    print_big_cost();


    //reduce by (p2 - 1)/pUUU = (p-1)*cofacUUU which can be done by Step 1: power (p-1) and Step 2: result to cofacUUU

    //Step 1
    let (asq, bsq) = (fp_sq(&fxre), fp_sq(&fxim));
    let abinv = non_constant_inversion(&fp_add(&asq, &bsq));
    let mut zetare = fp_mul(&fp_sub(&asq, &bsq), &abinv);

    //print_big_cost();


    //Step 2:
    //this can be done with lucaspow because zeta now has norm 1
    zetare = lucaspow(&zetare, &*cofacU);

    //print_big_cost();

    //compute order of zeta, if larger than bound we are good to go
    has_large_order(&zetare)
}