//algorithm to find full torsion points, using the multipairing and tricks to quickly
//compute what torsion the randomly found points are missing, and to "fill it up"
//using an adaption of the gauss algorithm

//compares against the naive method of randomly sampling points on E and checking their
//order individually with prod tree

use std::collections::VecDeque;

use num_bigint::{BigUint, ToBigUint};
use num_traits::ToPrimitive;

use crate::{elligator::shitty_elligator, multimiller::bigmultimiller::multimiller, bigconstants::*, bigint::*, biglucaspairing::*};

pub fn prodvecmiss(vec: &VecDeque<i128>, miss: i128) -> BigUint 
{
    let prod = vec.iter().fold(1, |acc, &x| acc * x);
    return (prod/miss).to_biguint().unwrap();
}

pub fn order_to_vec(order: &BigUint) -> VecDeque<i128>{
    //for now a fast function that, given an order, computes what the missing torsion is

    let mut rem =  (&*PPLUSONE / order).to_u32().unwrap() as i128;
    let mut res= VecDeque::new();

    let mut i:i128 = 2;
    while rem > 1 {
        while rem % i == 0 {
            rem = rem/i;
            res.push_back(i);
        }
        i += 1;
    }

    while res[0] == 2 {
        res.pop_front();
    }
    println!("res: {:?}", res);
    res

}


pub fn fastfinder(mont_A: &BigUint) -> (BigUint, BigUint, BigUint, BigUint) {
    //given a montgomery coefficient, finds with probabilistic method
    //two points P and Q that form a basis for E[Elkies], one on the curve
    //and one of the twist. This is equivalent to finding a point T+ and T- as
    //used in two=point strategies in CTIDH and CSIDH-variants.

    //one can do two or three pairings to start off with, which is a trade-off between
    //speed and probability of success

    //on the other hand, if memory is not an issue, one can apply the Scott-Miller approach
    //of seperating the computation of the lines function, from their evaluation in the points Q1, Q1, ...

    //first sample quasi-random points on E and twist using elligator
    let (p1x, p1y, q1x, q1y) = shitty_elligator(1, &mont_A);
    let (_p2x, _p2y, q2x, q2y) = shitty_elligator(2, &mont_A);

    let qx = vec![ &q1x, &q2x];
    let qy = vec![ &q1y, &q2y];

    //we compute the tate pairing
    let f = multimiller(&*PPLUSONE, &p1x, &p1y, qx, qy, &mont_A);

    //we reduce the result with the trick
    let (asq1, bsq1) = (fp_sq(&f[0].re), fp_sq(&f[0].im));
    let ab1 = fp_add(&asq1, &bsq1);

    let (asq2, bsq2) = (fp_sq(&f[1].re), fp_sq(&f[1].im));
    let ab2 = fp_add(&asq2, &bsq2);
    
    //we share the inversion to compute zeta1 and zeta2
    let mut biginv = fp_mul(&ab1, &ab2);
    biginv = fp_inv(&biginv);
    let ab1inv = fp_mul(&biginv, &ab2);
    let ab2inv = fp_mul(&biginv, &ab1);
    
    //these are the reduced results of the pairing, norm 1 elements of fp2 so we only 
    //need the real part to compute the orders
    let zeta1 = fp_mul(&fp_sub(&asq1, &bsq1), &ab1inv);
    let mut zeta2 = fp_mul(&fp_sub(&asq2, &bsq2), &ab2inv);

    //first, we compute the order of zeta1, we get some set of "missing" torsion
    //which describes which ell_i P1 OR R1 is missing
    let order1: BigUint = lucascompute_order_product_tree(&zeta1);
    let rem_ord1 = &order_to_vec(&order1);

    let mut m = (*BIGONE).clone();
    //kill the order1 part of zeta 2, and check if it has any of the remaining order
    zeta2 = lucaspow(&zeta2, &order1);
    for ord in rem_ord1 {
        let check = lucaspow(&zeta2, &prodvecmiss(&rem_ord1, *ord));
        if check != *BIGONE {
            m *= ord.to_biguint().unwrap();
        }
    }

    //we assume that the remaining order is what is missing in P, not in Q1 and Q2 both
    //so we should set Q1 to Q1 + ( p + 1 / order1) Q2
    //and P to P1 + (p + 1 / m) P2
    //(this is the "Gauss step")
    println!("remaining order: {}", m);

    (p1x, p1y, q1x, q1y)
}