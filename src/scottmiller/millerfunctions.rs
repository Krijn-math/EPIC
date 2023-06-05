//took out of bigscottmiller to make it cleaner

use super::constructionhelperfunctions::*;
use super::evaluationhelperfunctions::*;
use crate::bigconstants::*;
use crate::bigint::*;
use num_bigint::BigUint;

pub fn double(t: Point, mont_A: &BigUint) -> (Point, Line) {
    //on input t shoud give the point 2t
    //and the line ell_{t, t}


    //check if t is inf or two
    if is_inf(&t.xz, &t.yz, &t.z2) | is_two(&t.xz, &t.yz, &t.z2) {
        //square f, set t to INF
        let x2 = t.x2.clone();
        let (xz, yz, z2) = ((*BIGZERO).clone(), (*BIGONE).clone(), (*BIGZERO).clone());
        let t_res: Point = Point { x2, xz, yz, z2 };
        let f_pre =  Line { lx: (*BIGZERO).clone(),  ly: (*BIGZERO).clone(),  l0: (*BIGONE).clone(), } ;
        return (t_res, f_pre);
    } else {
        //double T, get line coeffs
        let (t_res, step) = dbl_and_line(t, mont_A);
        return (t_res, step);
    }
}


pub fn add(
    t: Point,
    p: &AffPoint,
    mont_A: &BigUint,
) -> (Point, Line) {
    //on input t, p shoud give the point t + p
    //and the line ell_{t, p}

    if is_inf(&t.xz, &t.yz, &t.z2) {
        let f_pre =  Line { lx: (*BIGZERO).clone(),  ly: (*BIGZERO).clone(),  l0: (*BIGONE).clone(), } ;
        return (t, f_pre); 
    }

    //if t is two, we do nothing will ell or f
    if is_two(&t.xz, &t.yz, &t.z2) {
        //add P to T, get line coeffs, here tx2 is not actually computed
        return add_and_line(t, &p, mont_A);
    }

    //add P to T, get line coeffs, here tx2 is not actually computed
    let (res_t, step) = add_and_line(t, &p, mont_A);

    if is_inf(&res_t.xz, &res_t.yz, &res_t.z2) {              
        let f_pre =  Line { lx: (*BIGZERO).clone(),  ly: (*BIGZERO).clone(),  l0: (*BIGONE).clone(), } ;
        return (res_t, f_pre);
    }

    (res_t, step)
}
