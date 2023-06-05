//given two points T+ and T- checks if they are a basis for p+1 torsion
//first computes the ReducedTatePairing of T+ and T-
//then checks if the result is a primitive root of unity

use num_bigint::{BigUint, ToBigUint};
use num_traits::One;

use crate::{bigconstants::*, bigint::*};
use crate::windowmiller::bigwindowmiller::windowmiller;

pub fn lucasPowOrderRec(are: &BigUint, low: usize, upp: usize, mut m: BigUint) -> BigUint {
    //variation of OrderRec in 2022/880 implemented specifically for roots of unity in fp2
    //only works for p + 1 div 4!!

    if upp - low == 0 {
        if is_one(are) {
            return m;
        }
        let k: BigUint = ELLS[upp].to_biguint().unwrap();
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
    let k: BigUint = ProdVec(&ELLS, low, mid);
    let lre: BigUint;
    lre = lucaspow(are, &k);
    m = lucasPowOrderRec(&lre, mid + 1, upp, m);

    let k: BigUint = ProdVec(&ELLS, mid + 1, upp);
    let rre: BigUint;
    rre = lucaspow(are, &k);

    lucasPowOrderRec(&rre, low, mid, m)
}

pub fn lucascompute_order_product_tree(are: &BigUint) -> BigUint {
    //computes the order of zeta = are + i*aim in F_p^2 as a root of unity of (p + 1)/
    //full order is (p + 1)/4, e.g. zeta^(p+1)/4 = 1
    //we always raise to the power 4 to remove this factor

    //uses product tree to be faster

    let a4re: BigUint;
    let mut res: BigUint = BigUint::one();
    (a4re) = lucaspow(are, &4.to_biguint().unwrap());

    res = lucasPowOrderRec(&a4re, 0, *NELLS - 1, res);

    res
}

pub fn lucaspairing_algorithm(
    tpx: &BigUint,
    tpy: &BigUint,
    tmx: &BigUint,
    tmy: &BigUint,
    mont_A: &BigUint,
) -> bool {

    //our degree r and power for reducing the tate pairing
    let r: BigUint = (*CSIDH512).clone() + &(*BIGONE);

    //we dont need this if we are clever
    //let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

    let f = windowmiller(&r, tpx, tpy, tmx, tmy, mont_A);
    print_big_cost();

    //we compute ^(p-1) by computing frob(f) and f^-1 and then multiply
    //frob(f) is practically free, f^-1 is one fp-inv so about 510 sq and 255 muls
    //however, in verification we can use the non constant xgcd
    let (asq, bsq) = (fp_sq(&f.re), fp_sq(&f.im));
    let abinv = non_constant_inversion(&fp_add(&asq, &bsq));

    //in lucaspairing, we dont need the imaginairy value of zeta, only the real part
    let zetare = fp_mul(&fp_sub(&asq, &bsq), &abinv);

    print_big_cost();
    reset_big_cost();
    let order: BigUint = lucascompute_order_product_tree(&zetare);
    print_big_cost();

    order == *SQRT
}

pub fn reduced_pairing(
    tpx: &BigUint,
    tpy: &BigUint,
    tmx: &BigUint,
    tmy: &BigUint,
    mont_A: &BigUint,
) -> BigUint {
    let r: BigUint = (*CSIDH512).clone() + &(*BIGONE);
    let f = windowmiller(&r, tpx, tpy, tmx, tmy, mont_A);
    let (asq, bsq) = (fp_sq(&f.re), fp_sq(&f.im));
    let abinv = non_constant_inversion(&fp_add(&asq, &bsq));
    let zetare = fp_mul(&fp_sub(&asq, &bsq), &abinv);
    zetare
}


#[cfg(test)]
mod tests {
    use crate::{tradmiller::bigmiller::miller};

    use super::*;

    #[test]
    fn verified_basis_CSIDH512() {
        let mont_A: BigUint = BigUint::parse_bytes(b"53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340", 16).unwrap();
        let tpx: BigUint = fp_neg(&BigUint::parse_bytes(b"17", 10).unwrap());
        let tmx: BigUint = BigUint::parse_bytes(b"16", 10).unwrap();
        let tpy: BigUint = BigUint::parse_bytes(b"B2EF33A6FACF4ED15A1AC03EA6C9473E08DAAC754D9CCC272C775D3B89A4D42DA1B8BB396430D3198FC61F45D5927F183424A5848FB23A374AB8C57B72E19F", 16).unwrap();
        let tmy: BigUint = BigUint::parse_bytes(b"28302633848ABDDFEC053E481685A2F4BF36C4405D9A5FB226160A1076AB579EC5C08E118E02E4DB1D5DF0F5C2E5DAEA84F3690086295C19A9A4EC3DE615A505", 16).unwrap();

        assert!(lucaspairing_algorithm(&tpx, &tpy, &tmx, &tmy, &mont_A));
    }

    #[test]
    fn lucaspow() {
        let mont_A: BigUint = BigUint::parse_bytes(b"53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340", 16).unwrap();
        let tpx: BigUint = fp_neg(&BigUint::parse_bytes(b"17", 10).unwrap());
        let tmx: BigUint = BigUint::parse_bytes(b"16", 10).unwrap();
        let tpy: BigUint = BigUint::parse_bytes(b"B2EF33A6FACF4ED15A1AC03EA6C9473E08DAAC754D9CCC272C775D3B89A4D42DA1B8BB396430D3198FC61F45D5927F183424A5848FB23A374AB8C57B72E19F", 16).unwrap();
        let tmy: BigUint = BigUint::parse_bytes(b"28302633848ABDDFEC053E481685A2F4BF36C4405D9A5FB226160A1076AB579EC5C08E118E02E4DB1D5DF0F5C2E5DAEA84F3690086295C19A9A4EC3DE615A505", 16).unwrap();

        let r: BigUint = (*CSIDH512).clone() + &(*BIGONE);

        let (fxre, fxim) = miller(&r, &tpx, &tpy, &tmx, &tmy, &mont_A);
        let (asq, bsq) = (fp_sq(&fxre), fp_sq(&fxim));
        let abinv = non_constant_inversion(&fp_add(&asq, &bsq));
        let (mut zetare, mut zetaim) = fp2sqr(&fxre, &fp_neg(&fxim));
        zetare = fp_mul(&zetare, &abinv);
        zetaim = fp_mul(&zetaim, &abinv);

        let fifteen: BigUint = 15.to_biguint().unwrap();

        assert!(lucaspow_both(&zetare, &zetaim, &fifteen) == fp2pow(&zetare, &zetaim, &fifteen));
    }
}
