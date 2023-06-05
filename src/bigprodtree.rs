//product tree algorithm for verification of full order points

use num_bigint::{BigUint, ToBigUint};
use num_traits::One;

use crate::{bigconstants::*, bigcurve::*, bigint::to_aff};

pub fn OrderRec(
    qx: &BigUint,
    qz: &BigUint,
    low: usize,
    upp: usize,
    mut m: BigUint,
    mont_A: &BigUint,
) -> BigUint {
    //variation of OrderRec in 2022/880 implemented specifically for VERIFICATION
    //only works for p + 1 div 4!!

    if upp - low == 0 {
        if is_zero(qz) {
            return m;
        }
        let _zx: BigUint;
        let zz: BigUint;
        let k: BigUint = ELLS[upp].to_biguint().unwrap();
        (_zx, zz) = ladder(&k, qx, qz, mont_A);
        if is_zero(&zz) {
            return &k * &m;
        }
        panic!(
            "THIS SHOULD NOT HAPPEN. ORDINARY CURVE upp = {}, low = {}, and q is {}",
            upp,
            low,
            to_aff(qx, qz)
        );
    }

    let mid: usize = mid(upp, low);
    let k: BigUint = ProdVec(&ELLS, low, mid);
    let lx: BigUint;
    let lz: BigUint;
    (lx, lz) = ladder(&k, qx, qz, mont_A);
    m = OrderRec(&lx, &lz, mid + 1, upp, m, mont_A);

    let k: BigUint = ProdVec(&ELLS, mid + 1, upp);
    let rx: BigUint;
    let rz: BigUint;
    (rx, rz) = ladder(&k, qx, qz, mont_A);

    OrderRec(&rx, &rz, low, mid, m, mont_A)
}

pub fn product_tree_algorithm(x: &BigUint, mont_A: &BigUint) -> bool {
    let mut qx: BigUint = (*x).clone();
    let mut qz: BigUint = BigUint::one();

    for _i in 0..cofac() {
        (qx, qz) = diff_dbl(&qx, &qz, mont_A);
    }

    let mut m: BigUint = BigUint::one();

    //this computes the order of x
    m = OrderRec(&qx, &qz, 0, *NELLS - 1, m, mont_A);
    m == ((&*CSIDH512 + &*BIGONE) >> 2)
}


pub fn OrderRecSS(
    qx: &BigUint,
    qz: &BigUint,
    low: usize,
    upp: usize,
    mut m: BigUint,
    mont_A: &BigUint,
) -> BigUint {
    //variation of OrderRec in 2022/880 implemented specifically for VERIFICATION
    //only works for p + 1 div 4!!

    if upp - low == 0 {
        if is_zero(qz) {
            return m;
        }
        let _zx: BigUint;
        let zz: BigUint;
        let k: BigUint = ELLS[upp].to_biguint().unwrap();
        (_zx, zz) = ladder(&k, qx, qz, mont_A);
        if is_zero(&zz) {
            return &k * &m;
        }
        panic!(
            "THIS SHOULD NOT HAPPEN. ORDINARY CURVE upp = {}, low = {}, and q is {}",
            upp,
            low,
            to_aff(qx, qz)
        );
    }

    let mid: usize = mid(upp, low);
    let k: BigUint = ProdVec(&ELLS, low, mid);
    let lx: BigUint;
    let lz: BigUint;
    (lx, lz) = ladder(&k, qx, qz, mont_A);
    m = OrderRecSS(&lx, &lz, mid + 1, upp, m, mont_A);

    if BigUint::ge(&m, &*BOUND)
    {
        return m;
    }

    let k: BigUint = ProdVec(&ELLS, mid + 1, upp);
    let rx: BigUint;
    let rz: BigUint;
    (rx, rz) = ladder(&k, qx, qz, mont_A);

    OrderRecSS(&rx, &rz, low, mid, m, mont_A)
}

pub fn product_tree_algorithmSS(x: &BigUint, mont_A: &BigUint) -> bool {
    let mut qx: BigUint = (*x).clone();
    let mut qz: BigUint = BigUint::one();

    for _i in 0..cofac() {
        (qx, qz) = diff_dbl(&qx, &qz, mont_A);
    }

    let mut m: BigUint = BigUint::one();

    //this computes the order of x
    m = OrderRecSS(&qx, &qz, 0, *NELLS - 1, m, mont_A);
    m == ((&*CSIDH512 + &*BIGONE) >> 2)
}




pub fn product_tree_cost() {

    println!("PRODUCT TREE ALGORITHM:");

    let mut mul = 0;
    let mut sq = 0;

    //clear cofac
    mul += *COST_DBL_MUL*cofac(); 
    sq += *COST_DBL_SQR*cofac();

    //ladders in prod tree
    mul += *COST_LDD_MUL*ellbits()*log2(*NELLS);
    sq += *COST_LDD_SQR*ellbits()*log2(*NELLS);

    println!("  -  ASYMPTOTIC COST: MULLS {}, SQUARES {}", mul, sq);
    println!("  -  COMPLEXITY: log(l) * dbls + log(n) * log(L) * ladders")
}