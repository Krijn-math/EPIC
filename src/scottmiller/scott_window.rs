//lots of intricate stuff, don't touch

use num_bigint::BigUint;

use crate::{bigconstants::*, bigint::*};

use super::bigscottmiller::*;
use super::evaluationhelperfunctions::*;
use super::millerfunctions::*;
use super::constructionhelperfunctions::*;

pub fn printf(f: &Result){
    println!("f is: {}", to_aff(&f.re, &f.im))
}

pub fn neg(p: &AffPoint) -> AffPoint{
    AffPoint { x: p.x.clone(), y: fp_neg(&p.y)}
}

pub fn precompute_const(
    px: &BigUint,
    py: &BigUint,
    mont_A: &BigUint,
) -> (
    Vec<AffPoint>,
    Vec<Line>,
) {
    //given p precomputes the points p1, p3, p5, p7 and p9 used in the window
    //also computes the required lines 

    let mut p: Vec<Point> = Default::default();


    let mut t = point_cast(&px, &py);
    let p_pre = AffPoint {  x: (*px).clone(), y: (*py).clone() } ;
    let mut step =  Line { lx: (*BIGZERO).clone(),  ly: (*BIGZERO).clone(),  l0: (*BIGONE).clone(), } ;
    let mut intstep: Line;

    let mut pres: Vec<AffPoint> = vec![ p_pre ];
    let mut f: Vec<&Line> = vec![ &step ];

    p.push(t);                  //contains in the end p1, p3, p5 etc in pro form
    f.push(&step);               //will contain the lines ell_{1, 1} needed in precompute_eval

    (t, step) = double(t, mont_A);         //2T
    (t, step) = add(t, &p_pre, mont_A);    //T becomes 3T, step here should be ell_(2,1)

    p.push(t); 
    f.push(&step);

    for i in 2..11 {
        (t, intstep) = add(t, &p_pre, mont_A);    
        (t, step) = add(t, &p_pre, mont_A);       
        step = combine(intstep, step);            
        p.push(t);
        f.push(&step);
    }

    let mut biginv = (*BIGONE).clone();

    for point in p {
        biginv = fp_mul(&biginv, &point.z2);
    }

    biginv = fp_inv(&biginv);              
    let mut prodA = (*BIGONE).clone();
    prodA = fp_mul(&prodA, &p[2].z2);
    prodA = fp_mul(&prodA, &p[3].z2);
    prodA = fp_mul(&prodA, &p[4].z2);
    prodA = fp_mul(&prodA, &p[5].z2);
    prodA = fp_mul(&prodA, &p[6].z2);

    let mut prodAA = (*BIGONE).clone();
    prodAA = fp_mul(&prodAA, &p[7].z2);
    prodAA = fp_mul(&prodAA, &p[8].z2);
    prodAA = fp_mul(&prodAA, &p[9].z2);

    let prodAB = fp_mul(&p[11].z2, &p[10].z2);

    let mut prodB = (*BIGONE).clone();
    prodB = fp_mul(&prodB, &p[7].z2);
    prodB = fp_mul(&prodB, &p[8].z2);
    prodB = fp_mul(&prodB, &p[9].z2);
    prodB = fp_mul(&prodB, &p[10].z2);
    prodB = fp_mul(&prodB, &p[11].z2);

    let mut prodBA = (*BIGONE).clone();
    prodBA = fp_mul(&prodBA, &p[2].z2);
    prodBA = fp_mul(&prodBA, &p[3].z2);
    prodBA = fp_mul(&prodBA, &p[4].z2);

    let prodBB = fp_mul(&p[6].z2, &p[5].z2);

    let prod35 = fp_mul(&p[3].z2, &p[2].z2);
    let prod1315 = fp_mul(&p[8].z2, &p[7].z2); //TODO swap around with above for saves

    let invA = fp_mul(&biginv, &prodA);
    let invAA = fp_mul(&invA, &prodAA); // has 19 21 left
    let invAB = fp_mul(&invA, &prodAB); // has 13 15 17 left

    let inv19 = fp_mul(&invAA, &p[11].z2);
    let inv21 = fp_mul(&invAA, &p[10].z2);
    let inv17 = fp_mul(&invAB, &prod1315);
    let inv1315 = fp_mul(&invAB, &p[9].z2);
    let inv13 = fp_mul(&inv1315, &p[8].z2);
    let inv15 = fp_mul(&inv1315, &p[7].z2);

    let invB = fp_mul(&biginv, &prodB);
    let invBA = fp_mul(&invB, &prodBA); // has 9 11 left
    let invBB = fp_mul(&invB, &prodBB); // has 3 5 7 left

    let inv9 = fp_mul(&invBA, &p[6].z2);
    let inv11 = fp_mul(&invBA, &p[5].z2);
    let inv7 = fp_mul(&invBB, &prod35);
    let inv35 = fp_mul(&invBB, &p[4].z2);
    let inv3 = fp_mul(&inv35, &p[3].z2);
    let inv5 = fp_mul(&inv35, &p[2].z2);

    
    let p1 = AffPoint { x: p[1].xz, y: p[1].yz };
    pres.push(p1);          //index 0 cannot be used so we add something useless
    pres.push(p1);
    pres.push(neg(&p1));

    let p3 = AffPoint {
        x: fp_mul(&p[2].xz, &inv3),
        y: fp_mul(&p[2].yz, &inv3),
    };

    pres.push(p3);
    pres.push(neg(&p3));

    let p5 = AffPoint {
        x: fp_mul(&p[3].xz, &inv5),
        y: fp_mul(&p[3].yz, &inv5),
    };
    
    pres.push(p5);
    pres.push(neg(&p5));

    let p7 = AffPoint {
        x: fp_mul(&p[4].xz, &inv7),
        y: fp_mul(&p[4].yz, &inv7),
    };

    pres.push(p7);
    pres.push(neg(&p7));

    let p9 = AffPoint {
        x: fp_mul(&p[5].xz, &inv9),
        y: fp_mul(&p[5].yz, &inv9),
    };

    pres.push(p9);
    pres.push(neg(&p9));

    let p11 = AffPoint {
        x: fp_mul(&p[6].xz, &inv11),
        y: fp_mul(&p[6].yz, &inv11),
    };

    pres.push(p11);
    pres.push(neg(&p11));

    let p13 = AffPoint {
        x: fp_mul(&p[7].xz, &inv13),
        y: fp_mul(&p[7].yz, &inv13),
    };

    pres.push(p13);
    pres.push(neg(&p13));

    let p15 = AffPoint {
        x: fp_mul(&p[8].xz, &inv15),
        y: fp_mul(&p[8].yz, &inv15),
    };

    pres.push(p15);
    pres.push(neg(&p15));

    let p17 = AffPoint {
        x: fp_mul(&p[9].xz, &inv17),
        y: fp_mul(&p[9].yz, &inv17),
    };

    pres.push(p17);
    pres.push(neg(&p17));

    let p19 = AffPoint {
        x: fp_mul(&p[10].xz, &inv19),
        y: fp_mul(&p[10].yz, &inv19),
    };

    pres.push(p19);
    pres.push(neg(&p19));

    let p21 = AffPoint {
        x: fp_mul(&p[11].xz, &inv21),
        y: fp_mul(&p[11].yz, &inv21),
    };
    
    pres.push(p21);
    pres.push(neg(&p21));

    (pres, f)
}

pub fn precompute_eval(
    fpre: Vec<Line>,
    rx: &BigUint,
    ry: &BigUint,
)   -> Vec<Result> {

    let mut f = Result{ re: (*BIGONE).clone(), im: (*BIGZERO).clone() };        //f1
    let res: Vec<Result> = vec![f] ;        //to make index 0 useless

    res.push(f);
    let mut minf = Result{ re: f.re, im: fp_neg(&f.im) };
    res.push(minf);

    let f2 = apply(fpre[1], &rx, &ry);
    
    for i in 2..fpre.len() {
        let eval = apply(fpre[i], rx, ry);
        f = absorb(f, eval);
        res.push(f);                             //should be f3, f5, etc.
        minf = Result{ re: f.re, im: fp_neg(&f.im) };
        res.push(minf);                          //should be f3min, f5min etc.
    }

    res
}