//testing supersingularity by first computing pairing, then checking order is > 4sqrt(p)

use num_bigint::{BigUint, ToBigUint};
use num_traits::One;

use crate::{bigconstants::*, bigint::*, windowmiller::bigwindowmiller::*, elligator::shitty_elligator};

pub fn random_curve_twist_points(_mont_A: &BigUint) -> (BigUint, BigUint, BigUint, BigUint) {
    //for now hardcoded, should use elligator 2
    let tpx: BigUint = fp_neg(&BigUint::parse_bytes(b"17", 10).unwrap());
    let tmx: BigUint = BigUint::parse_bytes(b"16", 10).unwrap();
    let tpy: BigUint = BigUint::parse_bytes(b"B2EF33A6FACF4ED15A1AC03EA6C9473E08DAAC754D9CCC272C775D3B89A4D42DA1B8BB396430D3198FC61F45D5927F183424A5848FB23A374AB8C57B72E19F", 16).unwrap();
    let tmy: BigUint = BigUint::parse_bytes(b"28302633848ABDDFEC053E481685A2F4BF36C4405D9A5FB226160A1076AB579EC5C08E118E02E4DB1D5DF0F5C2E5DAEA84F3690086295C19A9A4EC3DE615A505", 16).unwrap();
    
    (tpx, tpy, tmx, tmy)
}

pub fn lucasPowOrderRecSQRT(
    are: &BigUint,
    low: usize,
    upp: usize,
    mut m: BigUint,
) -> BigUint {
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
    m = lucasPowOrderRecSQRT(&lre, mid + 1, upp, m);

    if BigUint::ge(&m, &*BOUND)
    {
        return m;
    }

    let k: BigUint = ProdVec(&ELLS, mid + 1, upp);
    let rre: BigUint;
    rre = lucaspow(are, &k);

    lucasPowOrderRecSQRT(&rre, low, mid, m)
}

pub fn has_large_order(zetare: &BigUint) -> bool{
    let a4re: BigUint;
    (a4re) = lucaspow(zetare, &4.to_biguint().unwrap());
    let ord = lucasPowOrderRecSQRT(&a4re, 0, *NELLS - 1, BigUint::one());
    BigUint::ge(&ord, &*BOUND)
}

pub fn is_supersingular(mont_A: &BigUint) -> bool{
    let (px, py, qx, qy) = shitty_elligator(5, &mont_A);      //can also be implemented based on Elligator-2 map
    let f = windowmiller(&*PPLUSONE, &px, &py, &qx, &qy, &mont_A);
    let (asq, bsq) = (fp_sq(&f.re), fp_sq(&f.im));
    let abinv = non_constant_inversion(&fp_add(&asq, &bsq));
    let zetare = fp_mul(&fp_sub(&asq, &bsq), &abinv);

    has_large_order(&zetare)
}