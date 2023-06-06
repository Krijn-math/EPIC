//here we use an additional piece of information for faster verification of a torsion basis given in a public key
//we assume a primitive root of unity zeta is given as a system parameter, and the public key consists of 
//the montgomery coefficient A, the Elligator seed u, and the value of the discrete logarithm Log(zeta, zeta')
//where zeta' is the value of the reduced Tate pairing given Tp and Tm

use num_bigint::BigUint;
use num_bigint::{BigInt};

use crate::{bigint::*, bigconstants::*};

pub fn check_if_primitive(zetare: &BigUint, log: &BigUint) -> bool {

    //the primitive root can be given by just the a-value, and we can precompute for the Lucas exponentiation too
    let genre: BigUint = BigUint::parse_bytes(b"51A8DB5255FBE2EC8BF0C8DDDFA3F39A42EEF5271E66D200F11FA0124E506031D72C25A773C8C3A9ED8BFA3C52C613584A177C5C8B2854C2DAE8F54E62686517", 16).unwrap();

    assert!(lucaspow(&genre, log) == *zetare);
    let (d, _a, _b) = xgcd(log, &PPLUSONE);
    let e = BigInt::to_biguint(&d).unwrap();
    
    if (e == *BIGONE) | (e == *BIGTWO) | (e == fp_add(&BIGTWO, &BIGTWO)) {
       return true
    }
    
    false
}




//the imaginairy part can be computed as sqrt(1 - a^2) if necessary
//2710764576908909416734924963482085518589387842138928894596636637582471416431352469057118340502133833206580988396908790817709662181679892592067160446321968*i