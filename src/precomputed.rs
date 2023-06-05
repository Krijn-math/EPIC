use num_bigint::BigUint;

pub fn precomputed_root(i: usize) -> BigUint {
    //I just took 2 to 11 here, but I could imagine something more "random" would be better for if you actually want to use the points to perform isogenies

    let uu0 = BigUint::parse_bytes(b"2", 16).unwrap();
    let uu1 = BigUint::parse_bytes(b"3", 16).unwrap();
    let uu2 = BigUint::parse_bytes(b"4", 16).unwrap();
    let uu3 = BigUint::parse_bytes(b"5", 16).unwrap();
    let uu4 = BigUint::parse_bytes(b"6", 16).unwrap();
    let uu5 = BigUint::parse_bytes(b"7", 16).unwrap();
    let uu6 = BigUint::parse_bytes(b"8", 16).unwrap();
    let uu7 = BigUint::parse_bytes(b"9", 16).unwrap();
    let uu8 = BigUint::parse_bytes(b"A", 16).unwrap();
    let uu9 = BigUint::parse_bytes(b"B", 16).unwrap();

    let uu = [uu0, uu1, uu2, uu3, uu4, uu5, uu6, uu7, uu8, uu9].to_vec();
    uu[i].clone()
}


pub fn precomputed_factor(i: usize) -> BigUint {
    let alpha0 = BigUint::parse_bytes(b"21E6DA2FD15A833FFED8E59B1F6A196E3C0F02BE4F42D96B73A95442310B9899E28E4241CD511C57C5CD10440A591A61EB7B5EA6C7E43811B3D5E857114242D4", 16).unwrap();
    let alpha1 = BigUint::parse_bytes(b"3F90D919A889B617FDD6AE82DAE6EFAEB09C2524D49D57A978DD7DFC1BF5BE2088CABC3B60F8152492E07E7F9367117799875178B6CBE921313113A3405C3D4D", 16).unwrap();
    let alpha2 = BigUint::parse_bytes(b"6C7C53CC37880A6662B611F06486B7C726966F2DCA6F848B0BB7740703584EB93B60D405C436C118DF5D00D9BAB6BAD2F18AC87C1940B36BD912E77D04073C4", 16).unwrap();
    let alpha3 = BigUint::parse_bytes(b"3717228DB4331547FE20751C130C69532198647540CCA14E9BF328EB8FB2D7FA10272BAAEDA3CE0EA16D3A6E90D0CADF1EA879CF04D2DB1CC43B998D7C0BAC98", 16).unwrap();
    let alpha4 = BigUint::parse_bytes(b"2E7E6D0E60F129074EE050D4BD5E4EC310899D5CC478EFA94E2A0D2793B821BD1972A40278A9C0785FBA24EF9929BDC81E600CC7788971775D081A1101B9F54", 16).unwrap();
    let alpha5 = BigUint::parse_bytes(b"1B8B9146DA198AA3FF103A8E098634A990CC323AA06650A74DF99475C7D96BFD081395D576D1E70750B69D374868656F8F543CE782696D8E621DCCC6BE05D64C", 16).unwrap();
    let alpha6 = BigUint::parse_bytes(b"238421876D9BC673CE07CBF7D7C47C2A574CAD8A6B6A9AA1549902F0029E6F1B1E1B20D7389E1DB1491FF8A8CDE370F8DE507B8A26BE531EBC66305B36A6F0AD", 16).unwrap();
    let alpha7 = BigUint::parse_bytes(b"4D8CDFE6FBB23F6F309026D2D7DC4D5F5C958979BBB5B7BF656023F103642046697237901F432408C7AB7B9BA4789F8CCAAA3550B60D4042182C837A5DE12C11", 16).unwrap();
    let alpha8 = BigUint::parse_bytes(b"44D4B33A4BF5D42C98CFA39FCB6ACEC88973CF4453592DCAA50274C473178B28EB01CC570D762245995A85DF79DBB1AF88C42CC6FAB052CE9BAA8A25230A68A6", 16).unwrap();
    let alpha9 = BigUint::parse_bytes(b"33B33FEF5276D4F4CB0AC48C8FE83394E863B0FBD279252A43956D4B57981584464C250ABF821805DA725267C2FB1508871C238B2408D5816573025193EB72B6", 16).unwrap();

    let alpha = [alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9].to_vec();
    alpha[i].clone()
}

pub fn precomputed_pair(i: usize) -> (BigUint, BigUint) {
    (precomputed_root(i), precomputed_factor(i))
}









