//given two points T+ and T- checks if they are a basis for p+1 torsion
//first computes the ReducedTatePairing of T+ and T-
//then checks if the result is a primitive root of unity

use num_bigint::{BigUint, ToBigUint};
use num_traits::One;

use crate::{bigconstants::*, bigint::*, bigmiller::*, inshallah};

pub fn compute_order_naive(are: &BigUint, aim: &BigUint) -> BigUint {
    //computes the order of zeta = are + i*aim in F_p^2 as a root of unity of (p + 1)/
    //full order is (p + 1)/4, e.g. zeta^(p+1)/4 = 1
    //we always raise to the power 4 to remove this factor

    let a4re: BigUint;
    let a4im: BigUint;
    let mut res: BigUint = BigUint::one();
    (a4re, a4im) = fp2pow(are, aim, &4.to_biguint().unwrap());

    for i in 0..*NELLS {
        let k: BigUint = ((&*CSIDH512 + &*BIGONE) >> 2) / ELLS[i].to_biguint().unwrap();
        let resre: BigUint;
        let resim: BigUint;
        (resre, resim) = fp2pow(&a4re, &a4im, &k);
        if !is_fp2_one(&resre, &resim) {
            res *= &ELLS[i].to_biguint().unwrap();
        }
    }

    res
}

pub fn PowOrderRec(
    are: &BigUint,
    aim: &BigUint,
    low: usize,
    upp: usize,
    mut m: BigUint,
) -> BigUint {
    //variation of OrderRec in 2022/880 implemented specifically for roots of unity in fp2
    //only works for p + 1 div 4!!

    if upp - low == 0 {
        if is_fp2_one(are, aim) {
            return m;
        }
        let zre: BigUint;
        let zim: BigUint;
        let k: BigUint = ELLS[upp].to_biguint().unwrap();
        (zre, zim) = fp2pow(are, aim, &k);
        if is_fp2_one(&zre, &zim) {
            return &k * &m;
        }
        panic!(
            "THIS CANNOT LEGALLY HAPPEN upp = {}, low = {}, and a is {} + i*{}",
            upp, low, &are, &aim
        );
    }

    let mid: usize = mid(upp, low);
    let k: BigUint = ProdVec(&ELLS, low, mid);
    let lre: BigUint;
    let lim: BigUint;
    (lre, lim) = fp2pow(are, aim, &k);
    m = PowOrderRec(&lre, &lim, mid + 1, upp, m);

    let k: BigUint = ProdVec(&ELLS, mid + 1, upp);
    let rre: BigUint;
    let rim: BigUint;
    (rre, rim) = fp2pow(are, aim, &k);

    PowOrderRec(&rre, &rim, low, mid, m)
}

pub fn compute_order_product_tree(are: &BigUint, aim: &BigUint) -> BigUint {
    //computes the order of zeta = are + i*aim in F_p^2 as a root of unity of (p + 1)/
    //full order is (p + 1)/4, e.g. zeta^(p+1)/4 = 1
    //we always raise to the power 4 to remove this factor

    //uses product tree to be faster

    let a4re: BigUint;
    let a4im: BigUint;
    let mut res: BigUint = BigUint::one();
    (a4re, a4im) = fp2pow(are, aim, &4.to_biguint().unwrap());

    res = PowOrderRec(&a4re, &a4im, 0, *NELLS - 1, res);

    res
}

pub fn pairing_algorithm(
    tpx: &BigUint,
    tpy: &BigUint,
    tmx: &BigUint,
    tmy: &BigUint,
    mont_A: &BigUint,
) -> bool {
    //our degree r and power for reducing the tate pairing
    let r: BigUint = (*CSIDH512).clone() + &(*BIGONE);
    let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

    //these should multiply to p^2 - 1 in general
    inshallah!(&r * &to_red == &*CSIDH512 * &*CSIDH512 - &*BIGONE);

    //miller: uses the fact that fp-values dont matter for zeta, as it gets raised to p - 1, hence we can
    //use anything we want for fz and ellz in the loop, don't even need it

    let (fxre, fxim) = miller(&r, tpx, tpy, tmx, tmy, mont_A);

    print_big_cost();
    //we compute ^(p-1) by computing frob(f) and f^-1 and then multiply
    //frob(f) is practically free, f^-1 is one fp-inv so about 510 sq and 255 muls
    //fp-inv can be improved with addition chains to about 500 sq and a few muls
    let (asq, bsq) = (fp_sq(&fxre), fp_sq(&fxim));
    let abinv = non_constant_inversion(&fp_add(&asq, &bsq));
    let (mut zetare, mut zetaim) = fp2sqr(&fxre, &fp_neg(&fxim));
    zetare = fp_mul(&zetare, &abinv);
    zetaim = fp_mul(&zetaim, &abinv);

    print_big_cost();

    let order: BigUint = compute_order_product_tree(&zetare, &zetaim);

    order == *SQRT
}

pub fn pairing_cost() {
    println!("PAIRING ALGORITHM:");

    let mut mul = 0;
    let mut sq = 0;

 
    //pairing per bit
    mul += 510*18;
    sq += 510*4;

    //pairing additional when bit == 1
    mul += 256*18;
    sq += 256*2;

    println!("  -  TATE PAIRING:   MULLS {}, SQUARES {}", mul, sq);

    //final power, per bit
    mul += 510*2;
    sq += 0;

    //final power, when bit == 1
    mul += 256*2;
    sq += 0;
    println!("  -  REDUCING F:     MULLS {}, SQUARES {}", mul, sq);

    //kill power of 4
    mul += 2*2;
    sq += 2*0;

    //check order using prod tree
    mul += 508*log2(*NELLS)*3 + 256*log2(*NELLS)*2;   

    sq += 0;

    println!("  -  FULL ALGORITHM: MULLS {}, SQUARES {}", mul, sq);

    println!("  -  COMPLEXITY: log(p)*pairing_cost + ( log(p) + log(n) * log(L) )*fp2pow")
}
