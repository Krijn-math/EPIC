//naive algorithm for verification of full order points

use num_bigint::BigUint;
use num_bigint::ToBigUint;
use num_traits::One;

use crate::bigconstants::*;
use crate::bigcurve::*;
use crate::bigint::*;

pub fn naive_algorithm(px: &BigUint, mont_A: &BigUint) -> bool {
    //verifies that the point p is of full order by confirming every p+1 div 4 div ell is not INF
    //naievest way to check full torsion
    //cost per ell_i with bit length of (p+1)/ell as b_i is 3M + 2S + b_i*(6M + 4S)
    //by rough approximation b_i = log(p), in reality a few bits less, for CSIDH-512 for largest ell_i differs by 10 bits on 512 total

    let mut qx: BigUint = (*px).clone();
    let mut qz: BigUint = BigUint::one();

    for _i in 0..cofac() {
        (qx, qz) = diff_dbl(&qx, &qz, mont_A);
    }

    let p4 = to_aff(&qx, &qz);

    for i in 0..*NELLS {
        let k: BigUint = ((&*CSIDH512 + &*BIGONE) >> 2) / ELLS[i].to_biguint().unwrap();
        let x: BigUint;
        let z: BigUint;
        (x, z) = spe_ladder(&k, &p4, &mont_A);
        if is_zero(&x) | is_zero(&z) {
            return false;
        }
    }

    return true;

}

pub fn cost_naive() {
    //estimates cost given PRIME

    let mut mul = 0;
    let mut sq = 0;

    //clearing cofac
    mul += *COST_DBL_MUL*cofac();
    sq += *COST_DBL_SQR*cofac();

    //map to affine
    mul += 1;
    mul += inv_muls() as usize;
    sq += inv_sqs() as usize;

    for i in 0..ELLS.len() {
        mul += 3;
        sq += 2;

        let n: BigUint = ((&*CSIDH512 + &*BIGONE) >> 2) / ELLS[i].to_biguint().unwrap();
        let length = BigUint::bits(&n) as usize;
        mul += 6 * (length - 1);
        sq += 4 * (length - 1);
    }

    println!("NAIVE ALGORITHM:");
    println!("  -  PRECISE COST: MULLS {}, SQUARES {}", mul, sq);

    let mut mul = 0;
    let mut sq = 0;

    mul += *COST_DBL_MUL*(cofac() + 1) + *COST_SPE_MUL*ellbits()*(*NELLS - 1) - *COST_DBL_MUL*(*NELLS);
    sq += *COST_DBL_SQR*(cofac() + 1) + *COST_SPE_SQR*ellbits()*(*NELLS - 1) - *COST_DBL_SQR*(*NELLS);

    println!("  -  ASYMPTOTIC COST: MULLS {}, SQUARES {}", mul, sq);
    println!("  -  COMPLEXITY: log(l) * dbls + n * log(L) * ladders")
}

pub fn cost_naive_smallx() {
    //estimates cost given PRIME and assumes px is small enough to use adds

    let mut mul = 0;
    let mut sq = 0;

    for i in 0..ELLS.len() {
        mul += 3;
        sq += 2;

        let n: BigUint = ((&*CSIDH512 + &*BIGONE)) / ELLS[i].to_biguint().unwrap();
        let length = BigUint::bits(&n);
        mul += 5 * (length - 1);
        sq += 4 * (length - 1);
    }

    println!("NAIVE ALGORITHM WITH SMALL X:");
    println!("  -  PRECISE COST: MULLS {}, SQUARES {}", mul, sq);

    let mut mul = 0;
    let mut sq = 0;

    mul += (*COST_SPE_MUL - 1)*ellbits()*(*NELLS - 1) - *COST_DBL_MUL*(*NELLS);
    sq += *COST_SPE_SQR*ellbits()*(*NELLS - 1) - *COST_DBL_SQR*(*NELLS);

    println!("  -  ASYMPTOTIC COST: MULLS {}, SQUARES {}", mul, sq);
    println!("  -  COMPLEXITY: n * log(L) * spe_ladders")
}
