//contains everyting we need for x-only arithmetic on mont curves
use crate::bigconstants::*;
use crate::bigint::*;
use num_bigint::BigUint;

pub fn diff_add(
    px: &BigUint,
    pz: &BigUint,
    qx: &BigUint,
    qz: &BigUint,
    p_qx: &BigUint,
    p_qz: &BigUint,
    mont_A: &BigUint,
) -> (BigUint, BigUint) {
    //x-only differential addition of p and q, given p - q as (p_qx, p_qz)
    if is_zero(p_qx) & is_zero(p_qz) {
        return diff_dbl(px, pz, mont_A);
    }

    let mut v0 = fp_add(px, pz);
    let mut v1 = fp_sub(qx, qz);
    v1 = fp_mul(&v1, &v0);
    v0 = fp_sub(px, pz);

    let mut v2 = fp_add(qx, qz);
    v2 = fp_mul(&v2, &v0);

    let mut v3 = fp_add(&v1, &v2);
    v3 = fp_sq(&v3);

    let mut v4 = fp_sub(&v1, &v2);
    v4 = fp_sq(&v4);

    let rx = fp_mul(&v3, p_qz);
    let rz = fp_mul(&v4, p_qx);

    (rx, rz)
}

pub fn spe_add(
    px: &BigUint,
    pz: &BigUint,
    qx: &BigUint,
    qz: &BigUint,
    spex: &BigUint,
) -> (BigUint, BigUint) {
    //x-only differential addition of p and q, given p - q as (spex, 1)

    let mut v0 = fp_add(px, pz);
    let mut v1 = fp_sub(qx, qz);
    v1 = fp_mul(&v1, &v0);
    v0 = fp_sub(px, pz);

    let mut v2 = fp_add(qx, qz);
    v2 = fp_mul(&v2, &v0);

    let mut v3 = fp_add(&v1, &v2);
    v3 = fp_sq(&v3);

    let mut v4 = fp_sub(&v1, &v2);
    v4 = fp_sq(&v4);

    let rx = v3;
    let rz = fp_mul(&v4, spex);

    (rx, rz)
}

pub fn diff_dbl(px: &BigUint, pz: &BigUint, mont_A: &BigUint) -> (BigUint, BigUint) {
    //x-only doubling of p

    //if p is INF or order two, return INF
    if is_zero(px) {
        return ((*BIGZERO).clone(), (*BIGZERO).clone());
    }

    let mut v1 = fp_add(px, pz);
    v1 = fp_sq(&v1);

    let mut v2 = fp_sub(px, pz);
    v2 = fp_sq(&v2);

    let rx = fp_mul(&v1, &v2);

    v1 = fp_sub(&v1, &v2);

    let mut v3 = (&*TRICK + mont_A) >> 2;
    v3 = fp_mul(&v3, &v1);
    v3 = fp_add(&v3, &v2);

    let rz = fp_mul(&v1, &v3);
    (rx, rz)
}

pub fn ladder(n: &BigUint, px: &BigUint, pz: &BigUint, mont_A: &BigUint) -> (BigUint, BigUint) {
    //general montgomery ladder for p = (px : pz)
    if *pz == *BIGONE {
        return spe_ladder(n, px, mont_A);
    }

    let mut res0x = (*px).clone();
    let mut res0z = (*pz).clone();

    let mut res1x: BigUint;
    let mut res1z: BigUint;
    (res1x, res1z) = diff_dbl(&res0x, &res0z, mont_A);

    for i in (0..BigUint::bits(n) - 1).rev() {
        if BigUint::bit(n, i) {
            (res0x, res0z) = diff_add(&res0x, &res0z, &res1x, &res1z, px, pz, mont_A);
            (res1x, res1z) = diff_dbl(&res1x, &res1z, mont_A);
        } else {
            (res1x, res1z) = diff_add(&res0x, &res0z, &res1x, &res1z, px, pz, mont_A); //other way around res0 gets mutated before use
            (res0x, res0z) = diff_dbl(&res0x, &res0z, mont_A);
        }
    }

    (res0x, res0z)
}

pub fn spe_ladder(n: &BigUint, spex: &BigUint, mont_A: &BigUint) -> (BigUint, BigUint) {
    //specialized montgomery ladder for p = (spex : 1)

    let mut res0x = (*spex).clone();
    let mut res0z = (*BIGONE).clone();

    let mut res1x: BigUint;
    let mut res1z: BigUint;
    (res1x, res1z) = diff_dbl(&res0x, &res0z, mont_A);

    for i in (0..BigUint::bits(n) - 1).rev() {
        if BigUint::bit(n, i) {
            (res0x, res0z) = spe_add(&res0x, &res0z, &res1x, &res1z, spex);
            (res1x, res1z) = diff_dbl(&res1x, &res1z, mont_A);
        } else {
            (res1x, res1z) = spe_add(&res0x, &res0z, &res1x, &res1z, spex); //other way around res0 gets mutated before use
            (res0x, res0z) = diff_dbl(&res0x, &res0z, mont_A);
        }
    }

    (res0x, res0z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::RandBigInt;
    use rand;

    #[test]
    fn random_point() {
        let mont_A: BigUint = BigUint::parse_bytes(b"53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340", 16).unwrap();

        let mut rng = rand::thread_rng();
        let x = rng.gen_biguint(512) % &*CSIDH512;
        let pplusone = &*CSIDH512 + &(*BIGONE);
        let (_xx, zz) = ladder(&pplusone, &x, &*BIGONE, &mont_A);
        assert!(is_zero(&zz));
    }

    #[test]
    fn random_itself() {
        let mont_A: BigUint = BigUint::parse_bytes(b"53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340", 16).unwrap();

        let mut rng = rand::thread_rng();
        let x = rng.gen_biguint(512) % &*CSIDH512;
        let pplustwo = &*CSIDH512 + &(*BIGTWO);
        let (xx, zz) = ladder(&pplustwo, &x, &*BIGONE, &mont_A);
        println!("x: {}, xx: {}, zz: {}", x, xx, zz);
        assert!(x == to_aff(&xx, &zz));
    }
}