use crate::{bigconstants::*};
use num_bigint::{BigUint, ToBigUint, BigInt};

pub fn fp_mul(lhs: &BigUint, rhs: &BigUint) -> BigUint {
    if lhs == &*BIGONE {return rhs % &*CSIDH512; }
    if rhs == &*BIGONE {return lhs % &*CSIDH512; }
    if lhs == &*BIGZERO {return (*BIGZERO).clone(); }
    if rhs == &*BIGZERO {return (*BIGZERO).clone(); }

    if lhs == rhs {
        unsafe { BIGCOST[1] += 1 }
    } else {
        unsafe { BIGCOST[0] += 1 }
    }

    (lhs * rhs) % &*CSIDH512
}

pub fn fp_neg(val: &BigUint) -> BigUint {
    &*CSIDH512 - val
}

pub fn fp_add(lhs: &BigUint, rhs: &BigUint) -> BigUint {
    unsafe { BIGCOST[2] += 1 }

    (lhs + rhs) % &*CSIDH512
}

pub fn fp_sub(lhs: &BigUint, rhs: &BigUint) -> BigUint {
    fp_add(lhs, &fp_neg(rhs))
}

pub fn fp_pow(a: &BigUint, n: &BigUint) -> BigUint {
    let mut res = (*a).clone();

    for i in (0..BigUint::bits(n) - 1).rev() {
        if BigUint::bit(n, i) {
            res = fp_mul(&res, &res);
            res = fp_mul(&res, a);
        } else {
            res = fp_mul(&res, &res)
        }
    }

    res
}

pub fn fp_sq(a: &BigUint) -> BigUint {
    fp_mul(a, a)
}

pub fn fp_inv(a: &BigUint) -> BigUint {
    if *a == *BIGZERO {
        panic!("Division by zero")
    }

    fp_pow(a, &INV)
}

pub fn xgcd(m: &BigUint, n: &BigUint) -> (BigInt, BigInt, BigInt) {
    let mut v1: BigInt = (*m).clone().into();
    let mut v2: BigInt = (*BIGONE).clone().into();
    let mut v3: BigInt = (*BIGZERO).clone().into();

    let mut w1: BigInt = (*n).clone().into();
    let mut w2: BigInt = (*BIGZERO).clone().into();
    let mut w3: BigInt = (*BIGONE).clone().into();

    while w1.to_biguint().unwrap() != (*BIGZERO) {
        let x1 = &v1 - (&v1 / &w1)*&w1;
        let x2 = &v2 - (&v1 / &w1)*&w2;
        let x3 = &v3 - (&v1 / &w1)*&w3;

        v1 = w1; v2 = w2; v3 = w3;
        w1 = x1; w2 = x2; w3 = x3;
    }

    (v1, v2, v3)
}

pub fn non_constant_inversion(m : &BigUint) -> BigUint {
    let (d, _a, b) = xgcd(&*CSIDH512, &m);
    assert!(d.to_biguint().unwrap() == (*BIGONE));
    let mut c: BigInt = (*CSIDH512).clone().into();
    c = &b + &c;
    c.to_biguint().unwrap() % &*CSIDH512
}

pub fn inv_muls() -> u64 {
    BigUint::count_ones(&INV) - 1
}

pub fn inv_sqs() -> u64 {
    BigUint::bits(&INV) - 1
}

pub fn fp2mul(are: &BigUint, aim: &BigUint, bre: &BigUint, bim: &BigUint) -> (BigUint, BigUint) {
    let mut tt1 = fp_mul(are, bre);
    let mut tt2 = fp_mul(aim, bim);
    let t1 = fp_add(are, aim);
    let t2 = fp_add(bre, bim);

    let resre = fp_sub(&tt1, &tt2);
    tt1 = fp_add(&tt1, &tt2);
    tt2 = fp_mul(&t1, &t2);
    let resim = fp_sub(&tt2, &tt1);

    (resre, resim)
}

pub fn fp2sqr(are: &BigUint, aim: &BigUint) -> (BigUint, BigUint) {
    let t1 = fp_add(are, aim);
    let t2 = fp_sub(are, aim);
    let t3 = fp_add(are, are);
    let resre = fp_mul(&t1, &t2);
    let resim = fp_mul(&t3, aim);
    (resre, resim)
}

pub fn fp2pow(are: &BigUint, aim: &BigUint, n: &BigUint) -> (BigUint, BigUint) {
    let mut resre: BigUint = (*are).clone();
    let mut resim: BigUint = (*aim).clone();

    for i in (0..BigUint::bits(n) - 1).rev() {
        if BigUint::bit(n, i) {
            (resre, resim) = fp2sqr(&resre, &resim);
            (resre, resim) = fp2mul(&resre, &resim, are, aim);
        } else {
            (resre, resim) = fp2sqr(&resre, &resim);
        }
    }

    (resre, resim)
}

pub fn lucasseq(n: &BigUint, p: &BigUint) -> (BigUint, BigUint) {
   
    let two = BigUint::parse_bytes(b"2", 10).unwrap();
    let mut v0 = (two).clone();
    let mut v1 = (*p).clone();

    for i in (0..BigUint::bits(n)).rev() {
        if BigUint::bit(n, i) {
            v0 = fp_sub(&fp_mul(&v0, &v1), &p);
            v1 = fp_sub(&fp_sq(&v1), &two);
        } else {
            v1 = fp_sub(&fp_mul(&v0, &v1), &p);
            v0 = fp_sub(&fp_sq(&v0), &two);
        }    
    }

    (v0, v1)    //V_n(P, 1) and V_{n-1}(P-1)
}

pub fn lucaspow(are: &BigUint, n: &BigUint) -> BigUint {
    //computes the lucas sequence of (are, aim)
    //computes the real coordinate of (are + aim*i)^n

    let mut a: BigUint = (*are).clone();
    let p: BigUint = fp_add(&a, &a);

    let (vn, _vnmin) = lucasseq(&n, &p);    //we do not need v_{n-1}, just Re(zeta^m)

    if BigUint::bit(&vn, 0) {     //vn is odd, so can divide using + p and shift
        a = (&vn + &*CSIDH512) >> 1;
    } else {                                //vn is even, so can divide using shift
        a = &vn >> 1;
    }

    a
}

pub fn lucaspow_both(are: &BigUint, _aim: &BigUint, n: &BigUint) -> (BigUint, BigUint) {
    //computes the lucas sequence of (are, aim)
    //the inverse can be done with xGCD instead of fp_inv to be fast

    let mut a: BigUint = (*are).clone();
    let b: BigUint;
    let p: BigUint = fp_add(&a, &a);

    let (vk, _vkmin) = lucasseq(&n, &p);

    if BigUint::bit(&vk, 0) {     //vn2a is odd, so can divide using + p and shift
        a = (&vk + &*CSIDH512) >> 1;
    } else {            //vn2a is even, so can divide using shift
        a = &vk >> 1;
    }

    let mut b2 = fp_sq(&a);
    b2 = fp_sub( &*BIGONE, &b2);
    b = fp_pow(&b2, &*SQRT);

    assert!(a != *BIGONE);
    assert!(b2 == fp_sq(&b));
    assert!(is_one(&fp_add(&fp_sq(&a), &fp_sq(&b))));

    (a, b)
}

pub fn lucasonepow(are: &BigUint, n: &BigUint) -> bool {
    //computes the lucas sequence of (are, aim)
    //computes the real coordinate of (are + aim*i)^n

    let mut a: BigUint = (*are).clone();
    let p: BigUint = fp_add(&a, &a);

    let (vn2a, _vn2a1) = lucasseq(&n, &p);

    if BigUint::bit(&vn2a, 0) {     //vn2a is odd, so can divide using + p and shift
        a = (&vn2a + &*CSIDH512) >> 1;
    } else {            //vn2a is even, so can divide using shift
        a = &vn2a >> 1;
    }

    assert!(a == *BIGONE);

    true
}

pub fn bigF(x: u128) -> BigUint {
    x.to_biguint().unwrap()
}

pub fn is_inf(tx: &BigUint, ty: &BigUint, tz: &BigUint) -> bool {
    is_zero(tx) & !is_zero(ty) & is_zero(tz)
}

pub fn is_two(tx: &BigUint, ty: &BigUint, tz: &BigUint) -> bool {
    is_zero(tx) & is_zero(ty) & !is_zero(tz)
}

pub fn to_aff(tx: &BigUint, tz: &BigUint) -> BigUint {
    if is_zero(tx) {
        return (*BIGZERO).clone();
    }

    fp_mul(tx, &fp_inv(tz))
}

pub fn is_fp2_one(are: &BigUint, aim: &BigUint) -> bool {
    is_one(are) & is_zero(aim)
}

pub fn recover(tx: &BigUint, tz: &BigUint, mont_A: &BigUint) -> (BigUint, BigUint) {
    let tx2: BigUint = fp_sq(tx);
    let invz: BigUint = fp_inv(tz);
    let rhs1: BigUint = fp_mul(&fp_mul(tx, &tx2), &invz);
    let rhs2: BigUint = fp_mul(mont_A, &tx2);
    let rhs3: BigUint = fp_mul(tx, tz);

    let rhs: BigUint = fp_add(&rhs1, &fp_add(&rhs2, &rhs3));

    let ty: BigUint = fp_pow(&rhs, &SQRT);

    (fp_mul(tx, &invz), ty)
}

pub fn set_inf(mut _tx: &BigUint, mut _ty: &BigUint, mut _tz: & BigUint) {
    _tx = &*BIGZERO;
    _ty = &*BIGONE;
    _tz = &*BIGZERO;
    assert!(is_inf(_tx, _ty, _tz));
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::RandBigInt;
    use rand;

    #[test]
    fn random_inverse() {
        let mut rng = rand::thread_rng();
        let x = rng.gen_biguint(512) % &*CSIDH512;
        assert!(fp_mul(&fp_inv(&x), &x) == *BIGONE);
        assert!(fp_pow(&x, &*BIGTWO) == fp_mul(&x, &x))
    }

    #[test]
    fn random_negative() {
        let mut rng = rand::thread_rng();
        let x = rng.gen_biguint(512) % &*CSIDH512;
        assert!(fp_add(&x, &fp_mul(&x, &fp_neg(&*BIGONE))) == *BIGZERO);
    }

    #[test]
    fn random_distribute() {
        let mut rng = rand::thread_rng();
        let x = rng.gen_biguint(512) % &*CSIDH512;
        let y = rng.gen_biguint(512) % &*CSIDH512;
        let lhs = fp_pow(&fp_add(&x, &y), &*BIGTWO);
        let mut rhs = fp_pow(&x, &*BIGTWO);
        rhs = fp_add(&rhs, &fp_mul(&x, &y));
        rhs = fp_add(&rhs, &fp_mul(&x, &y));
        rhs = fp_add(&rhs, &fp_pow(&y, &*BIGTWO));

        assert!(lhs == rhs);
    }

    #[test]
    fn both_inversions() {
        let mut rng = rand::thread_rng();
        let x = rng.gen_biguint(512) % &*CSIDH512;
        let inv1 = fp_inv(&x); let inv2 = non_constant_inversion(&x);
        assert!(inv1 == inv2);
    }
}
