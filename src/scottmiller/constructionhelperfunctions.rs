//all helper functions in the miller loop

use num_bigint::BigUint;

use crate::{bigint::*, bigconstants::{BIGONE, BIGZERO}};

pub struct Point {
    pub x2: BigUint,
    pub xz: BigUint,
    pub z2: BigUint,
    pub yz: BigUint,
}

pub struct AffPoint {
    pub x: BigUint,
    pub y: BigUint,
}

#[derive(Clone)]
pub struct Result {
    pub re: BigUint,
    pub im: BigUint,
}

#[derive(Clone)]
pub struct Line {
    pub lx: BigUint,
    pub ly: BigUint,
    pub l0: BigUint,
}

pub fn point_cast(px: &BigUint, py: &BigUint) -> Point {
    let xz = (*px).clone();
    let yz = (*py).clone();
    let z2 = (*BIGONE).clone();
    let x2 = fp_sq(&xz);
    Point { x2, xz, z2, yz }
}

pub fn neg_pro(p: &Point) -> Point{
    Point { x2: p.x2.clone(), xz: p.xz.clone(), z2: p.z2.clone(), yz: fp_neg(&p.yz)}
}

pub fn neg(p: &AffPoint) -> AffPoint{
    AffPoint { x: p.x.clone(), y: fp_neg(&p.y)}
}

pub fn dbl_and_line(
    t: Point,
    mont_A: &BigUint,
) -> (
    Point,
    Line,
) {
    let mut x2: BigUint = (t.x2).clone();
    let mut xz: BigUint = (t.xz).clone();
    let mut z2: BigUint = (t.z2).clone();
    let mut yz: BigUint = (t.yz).clone();

    let mut lx: BigUint;
    let mut ly: BigUint;
    let mut l0: BigUint;
    let v0: BigUint;

    let mut xx2: BigUint;
    let t0: BigUint;

    xx2 = fp_add(&yz, &yz);
    ly = fp_sq(&xx2);
    l0 = fp_sub(&x2, &z2);
    v0 = fp_sq(&l0);
    l0 = fp_mul(&xx2, &l0);
    lx = fp_mul(&xz, &l0);
    xx2 = fp_mul(&yz, &ly);
    lx = fp_add(&xx2, &lx);
    yz = fp_add(&x2, &z2);
    yz = fp_mul(&mont_A, &yz);
    xx2 = fp_add(&xz, &xz);
    yz = fp_add(&xx2, &yz);
    yz = fp_add(&xx2, &yz);
    yz = fp_mul(&xx2, &yz);

    xx2 = fp_sq(&v0);
    t0 = fp_sq(&l0);
    z2 = fp_sq(&ly);
    yz = fp_add(&v0, &yz);
    yz = fp_mul(&l0, &yz);

    ly = fp_mul(&xz, &ly);
    l0 = fp_mul(&x2, &l0);

    x2 = xx2;
    xz = t0;

    (Point {x2, xz, z2, yz}, Line{lx, ly, l0})
}


pub fn add_and_line(
    t: Point,
    p: &AffPoint,
    mont_A: &BigUint,
) -> (Point, Line) {
    let mut lambdax: BigUint;
    let mut lambdaz: BigUint;
    let mut lx: BigUint;
    let mut l0: BigUint;

    lambdax = fp_mul(&p.y, &t.z2);
    lambdax = fp_sub(&lambdax, &t.yz);

    lambdaz = fp_mul(&p.x, &t.z2); //additional trick is that px can be low so implemented as adds
    lambdaz = fp_sub(&lambdaz, &t.xz);

    if lambdaz == *BIGZERO {
        let inf = Point {x2: (*BIGZERO).clone(), xz: (*BIGZERO).clone(), z2: (*BIGZERO).clone(), yz: (*BIGONE).clone()}; 

        let mut lambdax = fp_mul(&mont_A, &t.xz);
        lambdax = fp_add(&lambdax, &lambdax);
        lambdax = fp_add(&lambdax, &t.x2);
        lambdax = fp_add(&lambdax, &t.x2);
        lambdax = fp_add(&lambdax, &t.x2);
        lambdax = fp_add(&lambdax, &t.z2);

        let mut lx = fp_mul(&t.z2, &lambdax);

        let mut tmp = fp_mul(&t.yz, &t.yz);
        tmp = fp_add(&tmp, &tmp);
        let mut l0 = fp_mul(&lambdax, &t.xz);
        l0 = fp_sub(&l0, &tmp);

        let mut ly = fp_mul(&t.yz, &t.z2);
        ly = fp_add(&ly, &ly);

        let ly = (*BIGONE).clone();
        let lx = fp_neg(&(BIGONE));
        let l0 = (*BIGZERO).clone();

        return (inf, Line{ lx, ly, l0 })
    }


    let lambdaz2 = fp_sq(&lambdaz); //tmp1 = pow(lambdaz, 2);
    let tylambdaz2 = fp_mul(&t.yz, &lambdaz2);

    //Get ell = line(t, p, r) = (ry - ty - lambda(rx - tx)) / (rx + tx + px - lambda*lambda + a)
    //note here that rx and px are related for specific points as sampled by us, namely rx = -px

    let lambdaxlambdaz = fp_mul(&lambdax, &lambdaz);
    l0 = fp_mul(&t.xz, &lambdaxlambdaz);
    l0 = fp_sub(&l0, &tylambdaz2);

    lx = fp_mul(&t.z2, &lambdaxlambdaz);

    //it feels like we can improve this by using the rep (XX : XZ : ZZ : YZ)
    //right now we just use the (XZ : ZZ : YZ) part and fix it in the end with v
    //perhaps all of this can also be combined well with ellxre

    let tzlambdaz2 = fp_mul(&t.z2, &lambdaz2); //tmp2 = tz * tmp1;
    let lambdax2 = fp_sq(&lambdax); //tmp3 = pow(lambdax, 2);
    let txlambdaz2 = fp_mul(&t.xz, &lambdaz2); //tmp4 = tx * tmp1; //used in step 1.b and step 3
    let lambdaz3 = fp_mul(&lambdaz2, &lambdaz); //tmp5 = lambdaz * tmp1;
    let tzlambdax2 = fp_mul(&t.z2, &lambdax2); //tmp6 = tz * tmp3;

    let mut newtx = fp_add(&p.x, mont_A);
    newtx = fp_mul(&newtx, &fp_neg(&tzlambdaz2));
    newtx = fp_sub(&newtx, &txlambdaz2);
    newtx = fp_add(&newtx, &tzlambdax2);

    let tmp = fp_mul(&lambdaz3, &t.yz);
    let mut newty = fp_sub(&txlambdaz2, &newtx);
    newty = fp_mul(&newty, &lambdax);
    newty = fp_sub(&newty, &tmp);

    let x2 = fp_sq(&newtx);
    let xz = fp_mul(&tzlambdaz2, &newtx);
    let z2 = fp_sq(&tzlambdaz2);
    let mut yz = fp_mul(&lambdaz, &newty);
    yz = fp_mul(&t.z2, &yz);

    let ly = tzlambdaz2;

    (Point {x2, xz, z2, yz}, Line{ lx, ly, l0})
}
