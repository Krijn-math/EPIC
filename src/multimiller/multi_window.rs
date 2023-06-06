//lots of intricate stuff, don't touch

use num_bigint::BigUint;

use crate::{bigconstants::*, bigint::*};

use super::bigmultimillerhelper::*;
use super::bigmultimiller::*;

pub fn printf(f: &Result){
    println!("f is: {}", to_aff(&f.re, &f.im))
}

pub fn precompute(
    p: &AffPoint,
    rx: Vec<&BigUint>,
    ry: Vec<&BigUint>,
    mont_A: &BigUint,
) -> (
    AffPoint,
    AffPoint,
    AffPoint,
    AffPoint,
    AffPoint,
    AffPoint,
    AffPoint,
    AffPoint,
    AffPoint,
    AffPoint,
    AffPoint,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
    Vec<Result>,
) {
    //given p precomputes the points p1, p3, p5, p7 and p9 used in the window
    //also computes the required miller values f1, f3, f5, f7, f9 for addition and subtraction


    let mut t = point_cast(&p.x, &p.y);
    let mut f: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let f1: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f3: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f5: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f7: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f9: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f11: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f13: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f15: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f17: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f19: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f21: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];

    let mut f1min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f3min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f5min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f7min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f9min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f11min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f13min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f15min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f17min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f19min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];
    let mut f21min: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.clone().len() ];

    let p1 = point_cast(&p.x, &p.y);

    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //3T

    let p3 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz.clone(),
    };
    for i in 0..f.len(){
        f3[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }


    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //2T

    let p5 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz.clone(),
    };
    for i in 0..f.len(){
        f5[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }

    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //2T

    let p7 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz.clone(),
    };
    for i in 0..f.len(){
        f7[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }

    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //2T

    let p9 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz.clone(),
    };
    for i in 0..f.len(){
        f9[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }

    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //2T

    let p11 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz.clone(),
    };
    for i in 0..f.len(){
        f11[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }

    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //2T

    let p13 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz.clone(),
    };
    for i in 0..f.len(){
        f13[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }

    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //2T

    let p15 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz.clone(),
    };
    for i in 0..f.len(){
        f15[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }

    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //2T

    let p17 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz.clone(),
    };
    for i in 0..f.len(){
        f17[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }

    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //2T

    let p19 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz.clone(),
    };
    for i in 0..f.len(){
        f19[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }

    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A);
    (t, f) = add(t, p, f, f1.clone(), rx.clone(), ry.clone(), mont_A); //2T

    let p21 = Point {
        x2: t.x2.clone(),
        xz: t.xz.clone(),
        z2: t.z2.clone(),
        yz: t.yz,
    };
    for i in 0..f.len(){
        f21[i] = Result { re: f[i].re.clone(), im: f[i].im.clone(),  };
    }

    let mut biginv = (p1.z2).clone(); //<-- this is one, does not need to be killed
    biginv = fp_mul(&biginv, &p3.z2);
    biginv = fp_mul(&biginv, &p5.z2);
    biginv = fp_mul(&biginv, &p7.z2);
    biginv = fp_mul(&biginv, &p9.z2);
    biginv = fp_mul(&biginv, &p11.z2);
    biginv = fp_mul(&biginv, &p13.z2);
    biginv = fp_mul(&biginv, &p15.z2);
    biginv = fp_mul(&biginv, &p17.z2);
    biginv = fp_mul(&biginv, &p19.z2);
    biginv = fp_mul(&biginv, &p21.z2);
    biginv = non_constant_inversion(&biginv);

    //prod tree A kills 3 5 7 9 11
    //prod tree AA kills 13 15 17,     AB kills 19 21

    //prod tree B kills 13 15 17 19 21
    //prod tree BA kill 3 5 7,          BB kills 9 11

    let mut prodA = (*BIGONE).clone();
    prodA = fp_mul(&prodA, &p3.z2);
    prodA = fp_mul(&prodA, &p5.z2);
    prodA = fp_mul(&prodA, &p7.z2);
    prodA = fp_mul(&prodA, &p9.z2);
    prodA = fp_mul(&prodA, &p11.z2);

    let mut prodAA = (*BIGONE).clone();
    prodAA = fp_mul(&prodAA, &p13.z2);
    prodAA = fp_mul(&prodAA, &p15.z2);
    prodAA = fp_mul(&prodAA, &p17.z2);

    let prodAB = fp_mul(&p19.z2, &p21.z2);

    let mut prodB = (*BIGONE).clone();
    prodB = fp_mul(&prodB, &p13.z2);
    prodB = fp_mul(&prodB, &p15.z2);
    prodB = fp_mul(&prodB, &p17.z2);
    prodB = fp_mul(&prodB, &p19.z2);
    prodB = fp_mul(&prodB, &p21.z2);

    let mut prodBA = (*BIGONE).clone();
    prodBA = fp_mul(&prodBA, &p3.z2);
    prodBA = fp_mul(&prodBA, &p5.z2);
    prodBA = fp_mul(&prodBA, &p7.z2);

    let prodBB = fp_mul(&p9.z2, &p11.z2);

    let prod35 = fp_mul(&p3.z2, &p5.z2);
    let prod1315 = fp_mul(&p13.z2, &p15.z2); //TODO swap around with above for saves

    let invA = fp_mul(&biginv, &prodA);
    let invAA = fp_mul(&invA, &prodAA); // has 19 21 left
    let invAB = fp_mul(&invA, &prodAB); // has 13 15 17 left

    let inv19 = fp_mul(&invAA, &p21.z2);
    let inv21 = fp_mul(&invAA, &p19.z2);
    let inv17 = fp_mul(&invAB, &prod1315);
    let inv1315 = fp_mul(&invAB, &p17.z2);
    let inv13 = fp_mul(&inv1315, &p15.z2);
    let inv15 = fp_mul(&inv1315, &p13.z2);

    let invB = fp_mul(&biginv, &prodB);
    let invBA = fp_mul(&invB, &prodBA); // has 9 11 left
    let invBB = fp_mul(&invB, &prodBB); // has 3 5 7 left

    let inv9 = fp_mul(&invBA, &p11.z2);
    let inv11 = fp_mul(&invBA, &p9.z2);
    let inv7 = fp_mul(&invBB, &prod35);
    let inv35 = fp_mul(&invBB, &p7.z2);
    let inv3 = fp_mul(&inv35, &p5.z2);
    let inv5 = fp_mul(&inv35, &p3.z2);

    let p1 = AffPoint { x: p1.xz, y: p1.yz };
    let p3 = AffPoint {
        x: fp_mul(&p3.xz, &inv3),
        y: fp_mul(&p3.yz, &inv3),
    };
    let p5 = AffPoint {
        x: fp_mul(&p5.xz, &inv5),
        y: fp_mul(&p5.yz, &inv5),
    };
    let p7 = AffPoint {
        x: fp_mul(&p7.xz, &inv7),
        y: fp_mul(&p7.yz, &inv7),
    };
    let p9 = AffPoint {
        x: fp_mul(&p9.xz, &inv9),
        y: fp_mul(&p9.yz, &inv9),
    };
    let p11 = AffPoint {
        x: fp_mul(&p11.xz, &inv11),
        y: fp_mul(&p11.yz, &inv11),
    };
    let p13 = AffPoint {
        x: fp_mul(&p13.xz, &inv13),
        y: fp_mul(&p13.yz, &inv13),
    };
    let p15 = AffPoint {
        x: fp_mul(&p15.xz, &inv15),
        y: fp_mul(&p15.yz, &inv15),
    };
    let p17 = AffPoint {
        x: fp_mul(&p17.xz, &inv17),
        y: fp_mul(&p17.yz, &inv17),
    };
    let p19 = AffPoint {
        x: fp_mul(&p19.xz, &inv19),
        y: fp_mul(&p19.yz, &inv19),
    };
    let p21 = AffPoint {
        x: fp_mul(&p21.xz, &inv21),
        y: fp_mul(&p21.yz, &inv21),
    };

    for i in 0..f.len(){
        f1min[i] = Result {
            re: f1[i].re.clone(),
            im: fp_neg(&f1[i].im),
        }; //this trick is great
        f3min[i] = Result {
            re: f3[i].re.clone(),
            im: fp_neg(&f3[i].im),
        }; //the trick being that we only need f up to f_p value
        f5min[i] = Result {
            re: f5[i].re.clone(),
            im: fp_neg(&f5[i].im),
        }; //and this is (a - ib) up to fp value
        f7min[i] = Result {
            re: f7[i].re.clone(),
            im: fp_neg(&f7[i].im),
        };
        f9min[i] = Result {
            re: f9[i].re.clone(),
            im: fp_neg(&f9[i].im),
        };
        f11min[i] = Result {
            re: f11[i].re.clone(),
            im: fp_neg(&f11[i].im),
        };
        f13min[i] = Result {
            re: f13[i].re.clone(),
            im: fp_neg(&f13[i].im),
        };
        f15min[i] = Result {
            re: f15[i].re.clone(),
            im: fp_neg(&f15[i].im),
        };
        f17min[i] = Result {
            re: f17[i].re.clone(),
            im: fp_neg(&f17[i].im),
        };
        f19min[i] = Result {
            re: f19[i].re.clone(),
            im: fp_neg(&f19[i].im),
        };
        f21min[i] = Result {
            re: f21[i].re.clone(),
            im: fp_neg(&f21[i].im),
        };
    }


    (
        p1, p3, p5, p7, p9, p11, p13, p15, p17, p19, p21, f1, f1min, f3, f3min, f5, f5min, f7,
        f7min, f9, f9min, f11, f11min, f13, f13min, f15, f15min, f17, f17min, f19, f19min, f21,
        f21min,
    )
}

pub fn manual_window(
    t_in: Point,
    f_in: Vec<Result>,
    rx: Vec<&BigUint>,
    ry: Vec<&BigUint>,
    mont_A: &BigUint,
    p1: &AffPoint,
    p3: &AffPoint,
    p5: &AffPoint,
    p7: &AffPoint,
    p9: &AffPoint,
    p11: &AffPoint,
    p13: &AffPoint,
    p15: &AffPoint,
    p17: &AffPoint,
    p19: &AffPoint,
    p21: &AffPoint,
    f1: Vec<Result>,
    f1min: Vec<Result>,
    f3: Vec<Result>,
    f3min: Vec<Result>,
    f5: Vec<Result>,
    f5min: Vec<Result>,
    f7: Vec<Result>,
    f7min: Vec<Result>,
    f9: Vec<Result>,
    f9min: Vec<Result>,
    f11: Vec<Result>,
    _f11min: Vec<Result>,
    f13: Vec<Result>,
    f13min: Vec<Result>,
    f15: Vec<Result>,
    f15min: Vec<Result>,
    f17: Vec<Result>,
    f17min: Vec<Result>,
    f19: Vec<Result>,
    f19min: Vec<Result>,
    f21: Vec<Result>,
    f21min: Vec<Result>,
) -> Vec<Result> {
    let mut t = t_in;
    let mut f = f_in;

    //1       <-- can be ignored
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p13, f, f13min.clone(), rx.clone(), ry.clone(), mont_A);           //SUB 13
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p19, f, f19min.clone(), rx.clone(), ry.clone(), mont_A);           //SUB 19
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p9, f, f9.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 9
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p3, f, f3min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 3
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = add(t, p1, f, f1.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 1
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p9, f, f9min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 9
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = add(t, p1, f, f1.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 1
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = add(t, p1, f, f1.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 1
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p15, f, f15min.clone(), rx.clone(), ry.clone(), mont_A);           //SUB 15
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p7, f, f7.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 7
    (t, f) = n_double(11, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p7, f, f7min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 7
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p11, f, f11.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 11
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p5, f, f5min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 5
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p13, f, f13.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 13
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p11, f, f11.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 11
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p7, f, f7min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 7
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p7, f, f7min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 7
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p3, f, f3.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 3
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p19, f, f19.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 19
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p21, f, f21min.clone(), rx.clone(), ry.clone(), mont_A);           //SUB 21
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p17, f, f17.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 17
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p19, f, f19min.clone(), rx.clone(), ry.clone(), mont_A);           //SUB 19
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = add(t, p1, f, f1.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 1
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p15, f, f15.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 15
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p17, f, f17min, rx.clone(), ry.clone(), mont_A);           //SUB 17
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p9, f, f9min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 9
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p17, f, f17.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 17
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p3, f, f3.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 3
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = add(t, p1, f, f1.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 1
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p19, f, f19.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 19
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p5, f, f5min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 5
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = subtract(t, p1, f, f1min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 1
    (t, f) = n_double(7, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p13, f, f13min.clone(), rx.clone(), ry.clone(), mont_A);           //SUB 13
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p13, f, f13.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 13
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p19, f, f19.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 19
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p9, f, f9.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 9
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p5, f, f5min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 5
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p5, f, f5.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 5
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p3, f, f3min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 3
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p19, f, f19min.clone(), rx.clone(), ry.clone(), mont_A);           //SUB 19
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = add(t, p1, f, f1.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 1
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p21, f, f21min.clone(), rx.clone(), ry.clone(), mont_A);           //SUB 21
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p5, f, f5min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 5
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p7, f, f7.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 7
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p7, f, f7min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 7
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p21, f, f21min.clone(), rx.clone(), ry.clone(), mont_A);           //SUB 21
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = add(t, p1, f, f1.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 1
    (t, f) = n_double(6, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p13, f, f13min, rx.clone(), ry.clone(), mont_A);           //SUB 13
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p21, f, f21.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 21
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = add(t, p1, f, f1.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 1
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p7, f, f7.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 7
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p5, f, f5.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 5
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p11, f, f11.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 11
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p7, f, f7.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 7
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p3, f, f3, rx.clone(), ry.clone(), mont_A);                     //ADD 3
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p13, f, f13.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 13
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = subtract(t, p1, f, f1min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 1
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = add(t, p1, f, f1, rx.clone(), ry.clone(), mont_A);                     //ADD 1
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p15, f, f15min, rx.clone(), ry.clone(), mont_A);           //SUB 15
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p19, f, f19min, rx.clone(), ry.clone(), mont_A);           //SUB 19
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p15, f, f15, rx.clone(), ry.clone(), mont_A);                   //ADD 15
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p19, f, f19, rx.clone(), ry.clone(), mont_A);                   //ADD 19
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = subtract(t, p1, f, f1min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 1
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p5, f, f5, rx.clone(), ry.clone(), mont_A);                     //ADD 5
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p7, f, f7min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 7
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p7, f, f7.clone(), rx.clone(), ry.clone(), mont_A);                     //ADD 7
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p3, f, f3min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 3
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p11, f, f11, rx.clone(), ry.clone(), mont_A);                   //ADD 11
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p21, f, f21min, rx.clone(), ry.clone(), mont_A);           //SUB 21
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p21, f, f21.clone(), rx.clone(), ry.clone(), mont_A);                   //ADD 21
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p13, f, f13, rx.clone(), ry.clone(), mont_A);                   //ADD 13
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p9, f, f9, rx.clone(), ry.clone(), mont_A);                     //ADD 9
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p9, f, f9min, rx.clone(), ry.clone(), mont_A);             //SUB 9
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p7, f, f7, rx.clone(), ry.clone(), mont_A);                     //ADD 7
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p7, f, f7min.clone(), rx.clone(), ry.clone(), mont_A);             //SUB 7
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p21, f, f21, rx.clone(), ry.clone(), mont_A);                   //ADD 21
    (t, f) = n_double(1, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p3, f, f3min, rx.clone(), ry.clone(), mont_A);             //SUB 3
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p7, f, f7min, rx.clone(), ry.clone(), mont_A);             //SUB 7
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(3, t, f, rx.clone(), ry.clone(), mont_A);    
    (t, f) = subtract(t, p5, f, f5min, rx.clone(), ry.clone(), mont_A);             //SUB 5  
    (t, f) = n_double(2, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = n_double(5, t, f, rx.clone(), ry.clone(), mont_A);                   
    (t, f) = add(t, p17, f, f17, rx.clone(), ry.clone(), mont_A);  
    (t, f) = n_double(4, t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A); 
    (t, f) = subtract(t, p1, f, f1min, rx.clone(), ry.clone(), mont_A);             //SUB 1
    (t, f) = double(t, f, rx.clone(), ry.clone(), mont_A);

    f
}
