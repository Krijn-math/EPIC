//adapted elligator mapping from ctidh code

use num_bigint::{BigUint, ToBigUint};

use crate::{bigconstants::*, bigint::*, precomputed::*};

pub fn elligator(
    mont_A: &BigUint,
) -> (BigUint, BigUint, BigUint){
    let mut u: usize = 2;

    loop {
        u += 1;
        let u2 = fp_sq(&u.to_biguint().unwrap());
        let d = fp_sub(&u2, &BIGONE);
        if is_zero(&d) { continue }

        let mut m = fp_mul(&u2, mont_A);
        let mut t = fp_mul(&m, mont_A);
        let mut p = mont_A.clone();

        if is_zero(mont_A) {
            p = (*BIGONE).clone();
            m = (*BIGONE).clone();
            t = (*BIGONE).clone();
        }

        let D2 = fp_sq(&d);
        t = fp_add(&t, &D2);
        t = fp_mul(&t, &d);
        t = fp_mul(&t, &p);

        m = fp_neg(&m);

        return (p, m, d);
    }
}

pub fn rhs(x: &BigUint, mont_A: &BigUint) -> BigUint {
    let x2 = fp_sq(x);
    let tmp = fp_add(x, mont_A);
    let mut res = fp_mul(&x2, &tmp);
    res = fp_add(&res, x);
    res 
}

pub fn shitty_elligator(i: usize, mont_A: &BigUint) -> (BigUint, BigUint, BigUint, BigUint){
    //takes the precomputed random values in precompute.rs and returns you a point on the curve and the twist
    //uses precomputation so that we only need one exp to figure out which point is where and to get two affine points
    let (u, alpha) = precomputed_pair(i);
    let px = fp_mul(&alpha, mont_A);
    let mx = fp_neg(&fp_add(&px, mont_A));
    let rhspx = rhs(&px, mont_A);
    let tmpy = fp_pow(&rhspx, &SQRT);
    let tmpy2 = fp_sq(&tmpy);
    let py: BigUint;
    let my: BigUint;

    //if tmpy^2 == rhs then rhs had a root in Fp and so px is the point on E
    //otherwise, -tmpy^2 == rhs and so px is the point on the twist with py = i*tmpy
    //in both cases though, the y value of mx should bu i*u*yp
    
    if tmpy2 == rhspx {
        py = tmpy;
        my = fp_mul(&u, &py);

        //p on E, m on twist
        assert!(rhs(&px, mont_A) == fp_sq(&py));
        assert!(rhs(&mx, mont_A) == fp_neg(&fp_sq(&my)));
        (px, py, mx, my)
    } else {
        assert!(tmpy2 == fp_neg(&rhspx));
        py = tmpy;
        my = fp_neg(&fp_mul(&u, &py));           //py is i*something so multiplying with i*u gives the fp_neg

        //p on twist, m on E
        assert!(rhs(&px, mont_A) == fp_neg(&fp_sq(&py)));
        assert!(rhs(&mx, mont_A) == fp_sq(&my));
        (mx, my, px, py)
    }
}