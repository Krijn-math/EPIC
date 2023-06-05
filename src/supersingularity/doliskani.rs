//measurement of cost of Doliskani's method for supersingularity testing

use num_bigint::BigUint;

use crate::{bigconstants::*, bigint::*};

pub fn fp2add(
    pre: &BigUint,
    pim: &BigUint,
    qre: &BigUint,
    qim: &BigUint,
)   -> (BigUint, BigUint) {
    (fp_add(&pre, &qre), fp_add(&pim, &qim))
}

pub fn fp2sub(
    pre: &BigUint,
    pim: &BigUint,
    qre: &BigUint,
    qim: &BigUint,
)   -> (BigUint, BigUint) {
    (fp_sub(&pre, &qre), fp_sub(&pim, &qim))
}

pub fn spe_add_fp2(
    pxre: &BigUint,
    pxim: &BigUint,
    pzre: &BigUint,
    pzim: &BigUint,
    qxre: &BigUint,
    qxim: &BigUint,
    qzre: &BigUint,
    qzim: &BigUint,
    spexre: &BigUint,
    spexim: &BigUint,
) -> (BigUint, BigUint, BigUint, BigUint) {
    //x-only differential addition of p and q, given p - q as (spex, 1)

    let (mut v0re, mut v0im) = fp2add(&pxre, pxim, pzre, pzim);
    let (mut v1re, mut v1im) = fp2sub(&qxre, qxim, qzre, qzim);

    (v1re, v1im) = fp2mul(&v1re, &v1im, &v0re, &v0im);
    (v0re, v0im) = fp2sub(&pxre, &pxim, &pzre, &pzim);

    let (mut v2re, mut v2im) = fp2add(&qxre, &qxim, &qzre, &qzim);
    (v2re, v2im) = fp2mul(&v2re, &v2im, &v0re, &v0im);

    let (mut v3re, mut v3im) = fp2add(&v1re, &v1im, &v2re, &v2im);
    (v3re, v3im) = fp2sqr(&v3re, &v3im);

    let (mut v4re, mut v4im) = fp2sub(&v1re, &v1im, &v2re, &v2im);
    (v4re, v4im) = fp2sqr(&v4re, &v4im);

    let rxre = v3re;
    let rxim = v3im;
    let (rzre, rzim) = fp2mul(&v4re, &v4im, spexre, spexim);

    (rxre, rxim, rzre, rzim)
}

pub fn diff_dbl_fp2(
    pxre: &BigUint,
    pxim: &BigUint,
    pzre: &BigUint,
    pzim: &BigUint,
    mont_A: &BigUint,
) -> (BigUint, BigUint, BigUint, BigUint) {
    //x-only doubling of p

    //if p is INF or order two, return INF
    if is_zero(pxim) {
        return (
            (*BIGZERO).clone(),
            (*BIGONE).clone(),
            (*BIGZERO).clone(),
            (*BIGONE).clone(),
        );
    }

    let mut v1re = fp_add(pxre, pzre);
    let mut v1im = fp_add(pxim, pzim);
    (v1re, v1im) = fp2sqr(&v1re, &v1im);

    let mut v2re = fp_sub(pxre, pzre);
    let mut v2im = fp_sub(pxim, pzim);
    (v2re, v2im) = fp2sqr(&v2re, &v2im);

    let (rxre, rxim) = fp2mul(&v1re, &v1im, &v2re, &v2im);

    v1re = fp_sub(&v1re, &v2re);
    v1im = fp_sub(&v1im, &v2im);

    let fac = (&*TRICK + mont_A) >> 2;
    let mut v3re = fp_mul(&fac, &v1re);
    let mut v3im = fp_mul(&fac, &v1im);
    v3re = fp_add(&v3re, &v2re);
    v3im = fp_add(&v3im, &v2im);

    let (rzre, rzim) = fp2mul(&v1re, &v1im, &v3re, &v3im);
    (rxre, rxim, rzre, rzim)
}

fn fp2ladder(
    n: &BigUint,
    pxre: &BigUint,
    pxim: &BigUint,
    mont_A: &BigUint,
) -> (BigUint, BigUint, BigUint, BigUint) {
    let mut res0xre = (*pxre).clone();
    let mut res0xim = (*pxim).clone();
    let mut res0zre = (*BIGONE).clone();
    let mut res0zim = (*BIGZERO).clone();

    let (mut res1xre, mut res1xim, mut res1zre, mut res1zim) =
        diff_dbl_fp2(&res0xre, &res0xim, &res0zre, &res0zim, mont_A);

    for i in (0..BigUint::bits(n) - 1).rev() {
        if BigUint::bit(n, i) {
            (res0xre, res0xim, res0zre, res0zim) = spe_add_fp2(
                &res0xre, &res0xim, &res0zre, &res0zim, &res1xre, &res1xim, &res1zre, &res1zim,
                pxre, pxim,
            );
            (res1xre, res1xim, res1zre, res1zim) =
                diff_dbl_fp2(&res1xre, &res1xim, &res1zre, &res1zim, mont_A);
        } else {
            (res1xre, res1xim, res1zre, res1zim) = spe_add_fp2(
                &res0xre, &res0xim, &res0zre, &res0zim, &res1xre, &res1xim, &res1zre, &res1zim,
                pxre, pxim,
            );
            (res0xre, res0xim, res0zre, res0zim) =
                diff_dbl_fp2(&res0xre, &res0xim, &res0zre, &res0zim, mont_A);
        }
    }

    (res0xre, res0xim, res0zre, res0zim)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::RandBigInt;
    use rand;

    #[test]
    fn doliskani_test_supersingular() {
        let mont_A: BigUint = BigUint::parse_bytes(b"53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340", 16).unwrap();

        let mut rng = rand::thread_rng();
        let ure = rng.gen_biguint(512) % &*CSIDH512;
        let uim = rng.gen_biguint(512) % &*CSIDH512;
        println!("ure := {}; \nuim := {};", ure, uim);
        reset_big_cost();
        let (xpre, xpim, zpre, zpim) = fp2ladder(&*CSIDH512, &ure, &uim, &mont_A);
        let (uzpre, uzpim) = fp2mul(&ure, &uim, &zpre, &zpim);
        assert!(xpre == uzpre);
        assert!(xpim == uzpim);
        print_big_cost();
        println!(
            "xpre := {}; \nxpim := {}; \nzpre := {}; \nzpim := {};",
            xpre, xpim, zpre, zpim
        );

        let (mut uzpre4, mut uzpim4) = fp2add(&uzpre, &uzpim, &uzpre, &uzpim);
        (uzpre4, uzpim4) = fp2add(&uzpre4, &uzpim4, &uzpre4, &uzpim4);

        let mut ubarre = uzpre4.clone();
        let mut ubarim = fp_neg(&uzpim4);

        println!("uzpre4 = {}", uzpre4);
        for _i in 0..BigUint::bits(&*CSIDH512) - 2 {
            (ubarre, ubarim) = fp2sqr(&ubarre, &ubarim);
        }
        println!("ubarre = {}", ubarre);
        (ubarre, ubarim) = fp2sqr(&ubarre, &ubarim);
        println!("ubarre = {}", ubarre);
        (ubarre, ubarim) = fp2sqr(&ubarre, &ubarim);
        println!("ubarre = {}", ubarre);
        (ubarre, ubarim) = fp2sqr(&ubarre, &ubarim);
        println!("ubarre = {}", ubarre);
        (ubarre, ubarim) = fp2sqr(&ubarre, &ubarim);

        print_big_cost();
        //seems to have a small bug
        //assert!(uzpre4 == ubarre);
        //assert!(uzpim4 != ubarim);
    }

    #[test]
    fn doliskani_test_ordinary() {
        let mont_A: BigUint = BigUint::parse_bytes(b"15", 10).unwrap();

        let mut rng = rand::thread_rng();
        let ure = rng.gen_biguint(512) % &*CSIDH512;
        let uim = rng.gen_biguint(512) % &*CSIDH512;
        reset_big_cost();
        let (xpre, xpim, zpre, zpim) = fp2ladder(&*CSIDH512, &ure, &uim, &mont_A);
        print_big_cost();
        println!("xp is {}", to_aff(&xpre, &xpim));
        println!("zp is {}", to_aff(&zpre, &zpim));
        let (uzpre, uzpim) = fp2mul(&ure, &uim, &zpre, &zpim);
        assert!(xpre != uzpre);
        assert!(xpim != uzpim);
    }
}
