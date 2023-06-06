//here, I am doing one pairing computation applied to multiple points R_1, R_2, ..., R_n
//this will be useful for faster finding of full torsion points
//this is essentially bigwindowmiller.rs but adapted to multiple evaluations rx, ry


use crate::bigconstants::*;
use crate::bigint::*;
use crate::multimiller::bigmultimillerhelper::*;
use crate::multimiller::multi_window::*;
use num_bigint::BigUint;

pub fn double(
    t: Point,
    f: Vec<Result>,
    rx: Vec<&BigUint>,
    ry: Vec<&BigUint>,
    mont_A: &BigUint,
) -> (Point, Vec<Result>) {
    let mut f_res: Vec<Result> = f;
    //check if t is inf or two
    if is_inf(&t.xz, &t.yz, &t.z2) | is_two(&t.xz, &t.yz, &t.z2) {
        //square f, set t to INF
        let x2 = t.x2;
        let (xz, yz, z2) = ((*BIGZERO).clone(), (*BIGONE).clone(), (*BIGZERO).clone());
        let t_res = Point { x2, xz, yz, z2 };

        for i in 0..rx.len() {
            let (re, im) = fp2sqr(&f_res[i].re, &f_res[i].im);
            f_res[i] = Result { re, im };
        }

        (t_res, f_res)
    } else {
        //double T, get line coeffs
        let (t_res, (lambdax, lambday, lambdaz)) = dbl_and_line(t, mont_A);

        for i in 0..rx.len() {
            //evaluate line in rx, ry
            let (ellxre, ellxim) = give_line(&lambdax, &lambday, &lambdaz, rx[i], ry[i]);

            //square f and multiply by evaluation
            f_res[i] = square_and_absorb(f_res[i].clone(), &ellxre, &ellxim);
        }

        (t_res, f_res)
    }
}

pub fn n_double(
    n: usize,
    t: Point,
    f: Vec<Result>,
    rx: Vec<&BigUint>,
    ry: Vec<&BigUint>,
    mont_A: &BigUint,
) -> (Point, Vec<Result>) {
    let mut t_res = t;
    let mut f_res: Vec<Result> = f;

    for _i in 0..n {
        (t_res, f_res) = double(t_res, f_res, rx.clone(), ry.clone(), mont_A);      
    }

    (t_res, f_res)
}

pub fn add(
    t: Point,
    p: &AffPoint,
    f: Vec<Result>,
    fadd: Vec<Result>,
    rx: Vec<&BigUint>,
    ry: Vec<&BigUint>,
    mont_A: &BigUint,
) -> (Point, Vec<Result>) {
    let mut f_res: Vec<Result> = f;
    if is_inf(&t.xz, &t.yz, &t.z2) {
        return (t, f_res); 
    }

    //if t is two, we do nothing will ell or f
    if is_two(&t.xz, &t.yz, &t.z2) {
        //add P to T, get line coeffs, here tx2 is not actually computed
        let (t_res, (_lx, _lambday, _l0)) = add_and_line(t, p, mont_A);
        return (t_res, f_res);
    }

    //add P to T, get line coeffs, here tx2 is not actually computed
    let (res_t, (lx, lambday, l0)) = add_and_line(t, p, mont_A);
    
    if !is_inf(&res_t.xz, &res_t.yz,&res_t.z2){
    for i in 0..rx.len() {
        //evaluate line in rx, ry
        let (ellxre, ellxim) = give_line_bit(&lx, &lambday, &l0, rx[i], ry[i]);

        //absorb evaluation
        let res_f = absorb(f_res[i].clone(), &ellxre, &ellxim);
        
        let (resre, resim) = fp2mul(&res_f.re, &res_f.im, &fadd[i].re, &fadd[i].im);
        
        f_res[i] = Result{ re: resre, im: resim };
    }}


    (res_t, f_res)
}

pub fn subtract(
    t: Point,
    p: &AffPoint,
    f: Vec<Result>,
    fsub: Vec<Result>,
    rx: Vec<&BigUint>,
    ry: Vec<&BigUint>,
    mont_A: &BigUint,
) -> (Point, Vec<Result>) {
    add(t, &neg(p), f, fsub, rx, ry, mont_A)
}

pub fn multimiller(
    _n: &BigUint,
    px: &BigUint,
    py: &BigUint,
    rx: Vec<&BigUint>,
    ry: Vec<&BigUint>,
    mont_A: &BigUint,
) -> Vec<Result> {
    //cleanmiller but implemented using a window and NAF, using the fact that we can easily add or substract

    let t: Point = point_cast(px, py);
    let mut f: Vec<Result> = vec![ Result { re: (*BIGONE).clone(), im: (*BIGZERO).clone() }; rx.len() ];

    let p: AffPoint = AffPoint {
        x: (*px).clone(),
        y: (*py).clone(),
    };

    //precomputing the window
    let (p1, p3, p5, p7, p9, p11, p13, p15, p17, p19, p21, f1, f1min, f3, f3min, f5, f5min, f7, f7min, f9, f9min, f11, f11min, f13, f13min, f15, f15min, f17, f17min, f19, f19min, f21, f21min) =
        precompute(&p, rx.clone(), ry.clone(), mont_A);

    f = manual_window(t, f, rx, ry, mont_A, &p1, &p3, &p5, &p7, &p9, &p11, &p13, &p15, &p17, &p19, &p21, f1, f1min, f3, f3min, f5, f5min, f7, f7min, f9, f9min, f11, f11min, f13, f13min, f15, f15min, f17, f17min, f19, f19min, f21, f21min);
    
    for i in 0..f.len(){
        (f[i].re, f[i].im) = fp2sqr(&f[i].re, &f[i].im);        //last squaring instead of a doubling of T and sq f
    }

    f
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ezero_test() {
        //this is a verified example with magma for E0
        let mont_A: BigUint = (*BIGZERO).clone();

        let tmx: BigUint = BigUint::parse_bytes(b"7", 10).unwrap();
        let tmy: BigUint = BigUint::parse_bytes(b"2D6D8560086EF2C1C9ECB70B54BE94C5D4089051E5CD25E1C812AAFB4D0B876019A9BDC1259CCED43FC83FE5B4D277080BAB880649618832D87643A7111497D4", 16).unwrap();

        let tpx: BigUint = fp_neg(&tmx); //EXTRASPECIAL FORM
        let tpy: BigUint = tmy.clone(); //E0 form

        let rx = vec![ &tmx, &tmx];
        let ry = vec![ &tmy, &tmy];

        let r: BigUint = (*PPLUSONE).clone();
        let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

        //UPDATE
        let f = multimiller(&r, &tpx, &tpy, rx, ry, &mont_A);
        let (zetare, zetaim) = fp2pow(&f[0].re, &f[0].im, &to_red);

        //computed with magma
        let magmare: BigUint = BigUint::parse_bytes(b"43E0EDF3667B03F03E05B72AC7D87032CE217C9837B6A28E0BEB265AE92D37A7B694DEA79B3D06711631B3CF92B2E5F8FE537F630CACC01D98562E2A8DDB2AA", 16).unwrap();
        let magmaim: BigUint = BigUint::parse_bytes(b"20701C645F34AE64D163CB21503A8FA0D3B4521E890EEA504ED152C53709DCCD234AFC78A7B3308E1A48E0339FC75E9151A168D12A9501B3024AE552B9430804", 16).unwrap();
        assert!(zetare == magmare);
        assert!(zetaim == magmaim);
    }

    #[test]
    fn ethree_test() {
        //this is a verified example with magma for E0
        let mont_A: BigUint = BigUint::parse_bytes(b"53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340", 16).unwrap();

        let tmx: BigUint = BigUint::parse_bytes(b"2", 10).unwrap();
        let tmy: BigUint = BigUint::parse_bytes(b"15871D66F2CE52EE533800A4716AD8DEB15022BE55E48E041BB09FBD31BE148E164FFCF7B10E76E4AAC5946FD3C072765C14FEEFDE36893CB78B23161BFB6D1A", 16).unwrap();

        let tpx: BigUint = fp_neg(&tmx); //EXTRASPECIAL FORM
        let tpy: BigUint = BigUint::parse_bytes(b"12C9163DF1EC9421CF06E661549E2C9069143312489A7C7AD4F41BF2C2F8FB90954C96691F3691649FCB4FD89C6E398462F72D441778945E39D8414C654AB56D", 16).unwrap();

        let rx = vec![ &tmx, &tmx];
        let ry = vec![ &tmy, &tmy];

        let r: BigUint = (*PPLUSONE).clone();
        let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

        let f = multimiller(&r, &tpx, &tpy, rx, ry, &mont_A);
        let (zetare, zetaim) = fp2pow(&f[1].re, &f[1].im, &to_red);

        // //computed with magma
        let magmare: BigUint = BigUint::parse_bytes(b"7230A53FA464DF7559C7716823F50D114A2A96069A89F0027C01837A5BF7912642D0897BE621FADC7134A98E787C135501F462C1AEFE5DD1DDA8ED8BBE42865", 16).unwrap();
        let magmaim: BigUint = BigUint::parse_bytes(b"2892CA051611D66B7A0DE4470C028535C04A57870FFA9C5FAC75644774EDE7AB2C9728E62D54CC77E33D7DC61088F4D6C873D167D0651F620C3904AA374C8C69", 16).unwrap();

        assert!(zetare == magmare);
        assert!(zetaim == magmaim);
    }

    #[test]
    fn one_more_pairing() {
        //A is five 3-isogs awayfrom 0
        let mont_A: BigUint = BigUint::parse_bytes(b"314B6A52B5BCE8757BCC6AD02E2B2EA18CC9A4EAE3A7B22D5B6FB6B95B9EDCFC3330735D75D5F3DC229E79030D539AA332A5A6B9DE43B8999926C4780CBE1BE6", 16).unwrap();

        let tpx: BigUint = BigUint::parse_bytes(b"4555F8EBE1991B11F49794D988AA064E9F362F9EBF562A13F2019288E8FAA0E2866F305EBA78776A84EB18E7DBCF1CF6EEEE527035BAE27879C48B6AC8F4A772", 16).unwrap();
        let tpy: BigUint = BigUint::parse_bytes(b"2148F30D572F2CD6C3D309CC76086590C87142FA1DE6BDFE9CDE38247E1893E91A155065F071887BF404FAD294DD0D21FC39940A582B368E6E80D3A59BA09DC6", 16).unwrap();

        let r1x: BigUint = BigUint::parse_bytes(b"11AB41BF68B65E9675B416AA21256F4E925BA04B0A3E376DB0476DCB258A0CDE09F768D6018C691531088CDD6228C1D20CB8FD8307FA37B6C73B0C3CDCA7F499", 16).unwrap();
        let r1y: BigUint = BigUint::parse_bytes(b"AABDC7B7DB3C0D57DAE369270E985A1D86B6EAADE80027A8B38E4324EEC505C787BCA8F9F003FF070CFF859C3A50EB8D0E2F93157BD4645A59DAA774693ED0A", 16).unwrap();

        let r2x: BigUint = BigUint::parse_bytes(b"2B5037BC34E2CBF99BC92B29D6E779B0B4C3AEC2E8ACE8389EA04A58C124EB2CD4CD1D3C7DB4D4C3C55B365BFA1F12603D46E5438835D7959EF9E60FFC66D351", 16).unwrap();
        let r2y: BigUint = BigUint::parse_bytes(b"4E610802507D287484C84137F0934EDDAF81E3B72CF6E26A81220DFE76C94BB22EC42BCA1BECD4C50B1B51AE435CD8440515053C42858DF755F1C0B6257BAB5D", 16).unwrap();

        let rx = vec![ &r1x, &r2x];
        let ry = vec![ &r1y, &r2y];

        let r: BigUint = (*PPLUSONE).clone();
        let to_red: BigUint = (*CSIDH512).clone() - &(*BIGONE);

        let f = multimiller(&r, &tpx, &tpy, rx, ry, &mont_A);
        let (zetare1, zetaim1) = fp2pow(&f[0].re, &f[0].im, &to_red);
        let (zetare2, zetaim2) = fp2pow(&f[1].re, &f[1].im, &to_red);

        //computed with magma
        let magmare1: BigUint = BigUint::parse_bytes(b"565EB32F823FEC56C6F94ED82F471E566B8148DA7227194BCDB8CE62C36CDD162BA1CA873DF4DAE3E1C9C9B5B189957D86BC31A7D7395759E111E71B6B93EE5D", 16).unwrap();
        let magmaim1: BigUint = BigUint::parse_bytes(b"375F69FD90A60E31F5F1DE4EDC8DA971DB8DE118A6135E4D401B6FCABBCF9C7820D8912BA9BEE8B27650C1DB409CF86591D84C09E41288FDEF144B9D325BE81D", 16).unwrap();

        let magmare2: BigUint = BigUint::parse_bytes(b"2A810B8D4CDEDD1A6DA09B4DCC772C55EDFAF399E43BE437DE76F130FD04720221D4A649F28F7C6643E8B05B5DBD8F6EF26E7948F1918BC0451FAE9959C4C64", 16).unwrap();
        let magmaim2: BigUint = BigUint::parse_bytes(b"4C38636851F9D54CA25509AEA0FC5045BCA727E9DCC8E049A332FA1D842881AC1471DE591C7C6AAE43C2D2A2F41ACC66F84EC694567E84E6E5135C3CB1C9E7A2", 16).unwrap();

        assert!(zetare1 == magmare1);
        assert!(zetaim1 == magmaim1);
        assert!(zetare2 == magmare2);
        assert!(zetaim2 == magmaim2);
    }

    #[test]
    fn one_more_pairing_shared_inv() {
        //A is five 3-isogs awayfrom 0
        let mont_A: BigUint = BigUint::parse_bytes(b"314B6A52B5BCE8757BCC6AD02E2B2EA18CC9A4EAE3A7B22D5B6FB6B95B9EDCFC3330735D75D5F3DC229E79030D539AA332A5A6B9DE43B8999926C4780CBE1BE6", 16).unwrap();

        let tpx: BigUint = BigUint::parse_bytes(b"4555F8EBE1991B11F49794D988AA064E9F362F9EBF562A13F2019288E8FAA0E2866F305EBA78776A84EB18E7DBCF1CF6EEEE527035BAE27879C48B6AC8F4A772", 16).unwrap();
        let tpy: BigUint = BigUint::parse_bytes(b"2148F30D572F2CD6C3D309CC76086590C87142FA1DE6BDFE9CDE38247E1893E91A155065F071887BF404FAD294DD0D21FC39940A582B368E6E80D3A59BA09DC6", 16).unwrap();

        let r1x: BigUint = BigUint::parse_bytes(b"11AB41BF68B65E9675B416AA21256F4E925BA04B0A3E376DB0476DCB258A0CDE09F768D6018C691531088CDD6228C1D20CB8FD8307FA37B6C73B0C3CDCA7F499", 16).unwrap();
        let r1y: BigUint = BigUint::parse_bytes(b"AABDC7B7DB3C0D57DAE369270E985A1D86B6EAADE80027A8B38E4324EEC505C787BCA8F9F003FF070CFF859C3A50EB8D0E2F93157BD4645A59DAA774693ED0A", 16).unwrap();

        let r2x: BigUint = BigUint::parse_bytes(b"2B5037BC34E2CBF99BC92B29D6E779B0B4C3AEC2E8ACE8389EA04A58C124EB2CD4CD1D3C7DB4D4C3C55B365BFA1F12603D46E5438835D7959EF9E60FFC66D351", 16).unwrap();
        let r2y: BigUint = BigUint::parse_bytes(b"4E610802507D287484C84137F0934EDDAF81E3B72CF6E26A81220DFE76C94BB22EC42BCA1BECD4C50B1B51AE435CD8440515053C42858DF755F1C0B6257BAB5D", 16).unwrap();

        let rx = vec![ &r1x, &r2x];
        let ry = vec![ &r1y, &r2y];

        let f = multimiller(&*PPLUSONE, &tpx, &tpy, rx, ry, &mont_A);

        let (asq1, bsq1) = (fp_sq(&f[0].re), fp_sq(&f[0].im));
        let ab1 = fp_add(&asq1, &bsq1);

        let (asq2, bsq2) = (fp_sq(&f[1].re), fp_sq(&f[1].im));
        let ab2 = fp_add(&asq2, &bsq2);
        
        //we share the inversion to compute zeta1 and zeta2
        let mut biginv = fp_mul(&ab1, &ab2);
        biginv = fp_inv(&biginv);
        let ab1inv = fp_mul(&biginv, &ab2);
        let ab2inv = fp_mul(&biginv, &ab1);

        assert!(ab1inv == fp_inv(&fp_add(&asq1, &bsq1)));
        assert!(ab2inv == fp_inv(&fp_add(&asq2, &bsq2)));
        
        let zeta1 = fp_mul(&fp_sub(&asq1, &bsq1), &ab1inv);
        let zeta2 = fp_mul(&fp_sub(&asq2, &bsq2), &ab2inv);

        //computed with magma
        let magma1: BigUint = BigUint::parse_bytes(b"565EB32F823FEC56C6F94ED82F471E566B8148DA7227194BCDB8CE62C36CDD162BA1CA873DF4DAE3E1C9C9B5B189957D86BC31A7D7395759E111E71B6B93EE5D", 16).unwrap();
        let magma2: BigUint = BigUint::parse_bytes(b"2A810B8D4CDEDD1A6DA09B4DCC772C55EDFAF399E43BE437DE76F130FD04720221D4A649F28F7C6643E8B05B5DBD8F6EF26E7948F1918BC0451FAE9959C4C64", 16).unwrap();

        assert!(zeta1 == magma1);
        assert!(zeta2 == magma2);
    }
}
