// Dormand–Prince method 4(5)
use crate::SmpVector;

pub struct DOP54 {
    // deltaT
    dt: f64,
    abs_tol: f64,
    // c: Butcher tableau
    // c1: T,
    c2: f64,
    c3: f64,
    c4: f64,
    c5: f64,
    c6: f64,
    c7: f64,
    // A: Butcher tableau
    a21: f64,
    a31: f64,
    a32: f64,
    a41: f64,
    a42: f64,
    a43: f64,
    a51: f64,
    a52: f64,
    a53: f64,
    a54: f64,
    a61: f64,
    a62: f64,
    a63: f64,
    a64: f64,
    a65: f64,
    a71: f64,
    a72: f64,
    a73: f64,
    a74: f64,
    a75: f64,
    a76: f64,
    // b1: Butcher tableau
    b11: f64,
    b12: f64,
    b13: f64,
    b14: f64,
    b15: f64,
    b16: f64,
    b17: f64,
    // b1: Butcher tableau
    b21: f64,
    b22: f64,
    b23: f64,
    b24: f64,
    b25: f64,
    b26: f64,
    b27: f64,
    // dt scale factor
    dt_scale_up_factor: f64,
    dt_scale_down_factor: f64,
    dt_max: f64,
    // solver order
    order: f64,
    // for debug
    num_dt_switch: i64,
    num_steps: i64,
}

impl DOP54 {
    #[allow(dead_code)]
    pub fn new(
        dt1: f64,
        dt_max1: f64,
        abs_tol1: f64,
        dt_scale_up_factor1: f64,
        dt_sale_down_factor1: f64,
    ) -> DOP54 {
        DOP54 {
            // c: Butcher tableau
            dt: dt1,
            abs_tol: abs_tol1,
            //c1: 0.0,
            c2: 1.0 / 5.0,
            c3: 3.0 / 10.0,
            c4: 4.0 / 5.0,
            c5: 8.0 / 9.0,
            c6: 1.0,
            c7: 1.0,
            // A: Butcher tableau
            a21: 1.0 / 5.0,
            a31: 3.0 / 40.0,
            a32: 9.0 / 40.0,
            a41: 44.0 / 45.0,
            a42: -56.0 / 15.0,
            a43: 32.0 / 9.0,
            a51: 19372.0 / 6561.0,
            a52: -25360.0 / 2187.0,
            a53: 64448.0 / 6561.0,
            a54: -212.0 / 729.0,
            a61: 9017.0 / 3168.0,
            a62: -355.0 / 33.0,
            a63: 46732.0 / 5247.0,
            a64: 49.0 / 176.0,
            a65: -5103.0 / 18656.0,
            a71: 35.0 / 384.0,
            a72: 0.0,
            a73: 500.0 / 1113.0,
            a74: 125.0 / 192.0,
            a75: -2187.0 / 6784.0,
            a76: 11.0 / 84.0,
            // b1: Butcher tableau
            b11: 35.0 / 384.0,
            b12: 0.0,
            b13: 500.0 / 1113.0,
            b14: 125.0 / 192.0,
            b15: -2187.0 / 6784.0,
            b16: 11.0 / 84.0,
            b17: 0.0,
            // b1: Butcher tableau
            b21: 5179.0 / 57600.0,
            b22: 0.0,
            b23: 7571.0 / 16695.0,
            b24: 393.0 / 640.0,
            b25: -92097.0 / 339200.0,
            b26: 187.0 / 2100.0,
            b27: 1.0 / 40.0,
            // dt scale factor
            dt_scale_down_factor: dt_sale_down_factor1,
            dt_scale_up_factor: dt_scale_up_factor1,
            dt_max: dt_max1,
            // solver order
            order: 5.0,
            // for debug
            num_dt_switch: 0,
            num_steps: 0,
        }
    }

    #[inline(always)]
    pub fn solve_5th_order(
        &self,
        func: impl Fn(&SmpVector) -> SmpVector,
        x: &SmpVector,
    ) -> SmpVector {
        // 1
        let k1 = self.dt * &func(x);
        // 2
        let t21 = self.a21 * &k1;
        let mut t22 = x + &t21;
        t22.time = x.time + self.c2 * self.dt;
        let k2 = self.dt * &func(&t22);
        // 3
        let t31 = &(self.a31 * &k1) + &(self.a32 * &k2);
        let mut t32 = x + &t31;
        t32.time = x.time + self.c3 * self.dt;
        let k3 = self.dt * &func(&t32);
        // 4
        let t41 = &(self.a41 * &k1) + &(self.a42 * &k2);
        let t42 = &t41 + &(self.a43 * &k3);
        let mut t45 = x + &t42;
        t45.time = x.time + self.c4 * self.dt;
        let k4 = self.dt * &func(&t45);
        // 5
        let t51 = &(self.a51 * &k1) + &(self.a52 * &k2);
        let t52 = &t51 + &(self.a53 * &k3);
        let t53 = &t52 + &(self.a54 * &k4);
        let mut t54 = x + &t53;
        t54.time = x.time + self.c5 * self.dt;
        let k5 = self.dt * &func(&t54);
        // 6
        let t61 = &(self.a61 * &k1) + &(self.a62 * &k2);
        let t62 = &t61 + &(self.a63 * &k3);
        let t63 = &t62 + &(self.a64 * &k4);
        let t64 = &t63 + &(self.a65 * &k5);
        let mut t65 = x + &t64;
        t65.time = x.time + self.c6 * self.dt;
        let k6 = self.dt * &func(&t65);
        // 7
        let t71 = &(self.a71 * &k1) + &(self.a72 * &k2);
        let t72 = &t71 + &(self.a73 * &k3);
        let t73 = &t72 + &(self.a74 * &k4);
        let t74 = &t73 + &(self.a75 * &k5);
        let t75 = &t74 + &(self.a76 * &k6);
        let mut t76 = x + &t75;
        t76.time = x.time + self.c7 * self.dt;
        let k7 = self.dt * &func(&t76);
        // 5th order sol.
        let s1 = &(self.b11 * &k1) + &(self.b12 * &k2);
        let s2 = &s1 + &(self.b13 * &k3);
        let s3 = &s2 + &(self.b14 * &k4);
        let s4 = &s3 + &(self.b15 * &k5);
        let s5 = &s4 + &(self.b16 * &k6);
        let s6 = &s5 + &(self.b17 * &k7);
        //
        SmpVector {
            time: self.dt + x.time,
            vec: (x + &s6).vec,
        }
    }

    #[allow(dead_code)]
    #[inline(always)]
    pub fn solve_4th_order(
        &self,
        func: impl Fn(&SmpVector) -> SmpVector,
        x: &SmpVector,
    ) -> SmpVector {
        // 1
        let k1 = self.dt * &func(x);
        // 2
        let t21 = self.a21 * &k1;
        let mut t22 = x + &t21;
        t22.time = x.time + self.c2 * self.dt;
        let k2 = self.dt * &func(&t22);
        // 3
        let t31 = &(self.a31 * &k1) + &(self.a32 * &k2);
        let mut t32 = x + &t31;
        t32.time = x.time + self.c3 * self.dt;
        let k3 = self.dt * &func(&t32);
        // 4
        let t41 = &(self.a41 * &k1) + &(self.a42 * &k2);
        let t42 = &t41 + &(self.a43 * &k3);
        let mut t43 = x + &t42;
        t43.time = x.time + self.c4 * self.dt;
        let k4 = self.dt * &func(&t43);
        // 5
        let t51 = &(self.a51 * &k1) + &(self.a52 * &k2);
        let t52 = &t51 + &(self.a53 * &k3);
        let t53 = &t52 + &(self.a54 * &k4);
        let mut t54 = x + &t53;
        t54.time = x.time + self.c5 * self.dt;
        let k5 = self.dt * &func(&t54);
        // 6
        let t61 = &(self.a61 * &k1) + &(self.a62 * &k2);
        let t62 = &t61 + &(self.a63 * &k3);
        let t63 = &t62 + &(self.a64 * &k4);
        let t64 = &t63 + &(self.a65 * &k5);
        let mut t65 = x + &t64;
        t65.time = x.time + self.c6 * self.dt;
        let k6 = self.dt * &func(&t65);
        // 7
        let t71 = &(self.a71 * &k1) + &(self.a72 * &k2);
        let t72 = &t71 + &(self.a73 * &k3);
        let t73 = &t72 + &(self.a74 * &k4);
        let t74 = &t73 + &(self.a75 * &k5);
        let t75 = &t74 + &(self.a76 * &k6);
        let mut t76 = x + &t75;
        t76.time = x.time + self.c7 * self.dt;
        let k7 = self.dt * &func(&t76);
        // 5th order sol.
        let s1 = &(self.b21 * &k1) + &(self.b22 * &k2);
        let s2 = &s1 + &(self.b23 * &k3);
        let s3 = &s2 + &(self.b24 * &k4);
        let s4 = &s3 + &(self.b25 * &k5);
        let s5 = &s4 + &(self.b26 * &k6);
        let s6 = &s5 + &(self.b27 * &k7);
        //
        SmpVector {
            time: self.dt + x.time,
            vec: (x + &s6).vec,
        }
    }

    #[allow(dead_code)]
    #[inline(always)]
    pub fn solve(&mut self, func: impl Fn(&SmpVector) -> SmpVector, x: &SmpVector) -> SmpVector {
        // 1
        let k1 = self.dt * &func(x);
        // 2
        let t21 = self.a21 * &k1;
        let mut t22 = x + &t21;
        t22.time = x.time + self.c2 * self.dt;
        let k2 = self.dt * &func(&t22);
        // 3
        let t31 = &(self.a31 * &k1) + &(self.a32 * &k2);
        let mut t32 = x + &t31;
        t32.time = x.time + self.c3 * self.dt;
        let k3 = self.dt * &func(&t32);
        // 4
        let t41 = &(self.a41 * &k1) + &(self.a42 * &k2);
        let t42 = &t41 + &(self.a43 * &k3);
        let mut t45 = x + &t42;
        t45.time = x.time + self.c4 * self.dt;
        let k4 = self.dt * &func(&t45);
        // 5
        let t51 = &(self.a51 * &k1) + &(self.a52 * &k2);
        let t52 = &t51 + &(self.a53 * &k3);
        let t53 = &t52 + &(self.a54 * &k4);
        let mut t54 = x + &t53;
        t54.time = x.time + self.c5 * self.dt;
        let k5 = self.dt * &func(&t54);
        // 6
        let t61 = &(self.a61 * &k1) + &(self.a62 * &k2);
        let t62 = &t61 + &(self.a63 * &k3);
        let t63 = &t62 + &(self.a64 * &k4);
        let t64 = &t63 + &(self.a65 * &k5);
        let mut t65 = x + &t64;
        t65.time = x.time + self.c6 * self.dt;
        let k6 = self.dt * &func(&t65);
        // 7
        let t71 = &(self.a71 * &k1) + &(self.a72 * &k2);
        let t72 = &t71 + &(self.a73 * &k3);
        let t73 = &t72 + &(self.a74 * &k4);
        let t74 = &t73 + &(self.a75 * &k5);
        let t75 = &t74 + &(self.a76 * &k6);
        let mut t76 = x + &t75;
        t76.time = x.time + self.c7 * self.dt;
        let k7 = self.dt * &func(&t76);
        // 4th order sol.
        let s41 = &(self.b21 * &k1) + &(self.b22 * &k2);
        let s42 = &s41 + &(self.b23 * &k3);
        let s43 = &s42 + &(self.b24 * &k4);
        let s44 = &s43 + &(self.b25 * &k5);
        let s45 = &s44 + &(self.b26 * &k6);
        let s46 = &s45 + &(self.b27 * &k7);
        //
        let x_order_4 = SmpVector {
            time: self.dt + x.time,
            vec: (x + &s46).vec,
        };
        // 5th order sol.
        let s51 = &(self.b11 * &k1) + &(self.b12 * &k2);
        let s52 = &s51 + &(self.b13 * &k3);
        let s53 = &s52 + &(self.b14 * &k4);
        let s54 = &s53 + &(self.b15 * &k5);
        let s55 = &s54 + &(self.b16 * &k6);
        let s56 = &s55 + &(self.b17 * &k7);
        //
        let x_order_5 = SmpVector {
            time: self.dt + x.time,
            vec: (x + &s56).vec,
        };
        // adaptive dt
        let err_45 = (&x_order_5 - &x_order_4).norm();
        //println!("{:.14}", self.num_dt_switch);
        if err_45 > self.abs_tol {
            self.dt = self.dt_scale_down_factor
                * self.dt
                * (self.abs_tol / err_45).powf(1.0 / self.order);
            self.num_dt_switch += 1;
            self.num_steps += 1;
            x.clone()
        } else {
            self.dt *= self.dt_scale_up_factor;
            if self.dt > self.dt_max {
                self.dt = self.dt_max;
            }
            self.num_steps += 1;
            x_order_5
        }
    }

    #[allow(dead_code)]
    #[inline(always)]
    pub fn get_num_of_dt_switch(&self) -> i64 {
        self.num_dt_switch
    }

    #[allow(dead_code)]
    #[inline(always)]
    pub fn get_num_of_steps(&self) -> i64 {
        self.num_steps
    }

    #[allow(dead_code)]
    #[inline(always)]
    pub fn solve_to_end_time(
        &mut self,
        end_time: f64,
        func: impl Fn(&SmpVector) -> SmpVector + 'static + Copy,
        x: &SmpVector,
    ) -> SmpVector {
        let mut _y_new = x.clone();
        let mut y_old = x.clone();
        loop {
            _y_new = self.solve(func, &y_old);
            if _y_new.time > end_time {
                let to_end = end_time - y_old.time;
                self.set_delta_t(to_end);
                _y_new = self.solve_5th_order(func, &y_old);
                println!("{:.14} {:.14}", _y_new.time, _y_new.vec[0]);
                break;
            }
            y_old = _y_new.clone();
            println!("{:.14} {:.14}", _y_new.time, _y_new.vec[0]);
        }
        _y_new
    }

    pub fn set_delta_t(&mut self, dt: f64) {
        self.dt = dt;
    }

    #[allow(dead_code)]
    pub fn get_delta_t(&self) -> f64 {
        self.dt
    }
}

#[cfg(test)]
pub mod dop54_tests {
    use crate::num_core::dop54::DOP54;
    use crate::num_core::smp_vector::SmpVector;

    #[inline(always)]
    pub fn func1(x: &SmpVector) -> SmpVector {
        let ret = x.clone();
        ret.clone()
    }

    #[test]
    fn adaptive_time_steps1() {
        extern crate approx;

        // settings
        let init_dt = 1.0e-6;
        let dt_max = 0.2;
        let abs_tol = 1.0e-12;
        let dt_scale_up_factor = 1.005;
        let dt_scale_down_factor = 0.9;
        // ODE solver
        let mut sol = DOP54::new(
            init_dt,
            dt_max,
            abs_tol,
            dt_scale_up_factor,
            dt_scale_down_factor,
        );

        let mut y1 = SmpVector::set_one_fill(2);
        loop {
            let y2 = sol.solve(func1, &y1);
            y1 = y2;
            if y1.time > 10.0 {
                break;
            }
        }

        assert_eq!(sol.num_dt_switch, 98);
        assert_eq!(sol.num_steps, 4661);
        approx::assert_abs_diff_eq!(y1.time, 10.00061711082885, epsilon = 1.0e-14);
        approx::assert_abs_diff_eq!(y1.vec[0], 22040.062760356937, epsilon = 1.0e-14);
    }

    #[test]
    fn forth_order() {
        extern crate approx;

        // settings
        let init_dt = 1.0e-3;
        let dt_max = 0.2;
        let abs_tol = 1.0e-12;
        let dt_scale_up_factor = 1.001;
        let dt_scale_down_factor = 0.2;
        let s = DOP54::new(
            init_dt,
            dt_max,
            abs_tol,
            dt_scale_up_factor,
            dt_scale_down_factor,
        );

        let mut y1 = SmpVector::set_one_fill(1);
        for _i in 0..1000 {
            let y2 = s.solve_4th_order(func1, &y1);
            y1 = y2;
        }

        approx::assert_abs_diff_eq!(y1.vec[0], std::f64::consts::E, epsilon = 1.0e-11);
    }

    #[test]
    fn fifth_order() {
        extern crate approx;

        // settings
        let init_dt = 1.0e-3;
        let dt_max = 0.1;
        let abs_tol = 1.0e-12;
        let dt_scale_up_factor = 1.001;
        let dt_scale_down_factor = 0.2;
        let s = DOP54::new(
            init_dt,
            dt_max,
            abs_tol,
            dt_scale_up_factor,
            dt_scale_down_factor,
        );

        let mut y1 = SmpVector::set_one_fill(1);
        for _i in 0..1000 {
            let y2 = s.solve_5th_order(func1, &y1);
            y1 = y2;
        }

        approx::assert_abs_diff_eq!(y1.vec[0], std::f64::consts::E, epsilon = 1.0e-11);
    }

    #[test]
    fn calc_ode_napier() {
        // ODE solver
        let mut sol = DOP54::new(1.0e-1, 0.2, 1.0e-14, 1.003, 1.0);
        // initial vector
        let x = SmpVector::set_one_fill(100);
        let y1 = sol.solve_to_end_time(1.0, func1, &x);
        //
        approx::assert_abs_diff_eq!(y1.vec[99], std::f64::consts::E, epsilon = 1.0e-13);
    }
}
