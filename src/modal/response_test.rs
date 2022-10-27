#[cfg(test)]
pub mod test {
    extern crate nalgebra as na;
    use super::super::modes_test;
    use super::super::response::*;
    use na::DVector;
    use std::ops::Div;

    pub struct Response {
        pub modal: modes_test::test::Modal,
        pub free_transient: Transient,
    }

    pub struct Transient {
        pub init_cond: DVector<f64>,
        pub response: Vec<Snapshot>,
    }

    pub struct Snapshot {
        pub time: f64,
        pub values: DVector<f64>,
    }

    pub struct Setup {
        pub _2d: Response,
    }

    impl Setup {
        pub fn new() -> Self {
            let modal = modes_test::test::Setup::new();

            let response_t5s = Snapshot {
                time: 5.0,
                values: DVector::from_vec(vec![-0.0131029, -0.2396485]).div(100.0),
            };
            let response_t10s = Snapshot {
                time: 10.0,
                values: DVector::from_vec(vec![-0.7700099, 0.24229929]).div(100.0),
            };
            let hand_calc_sol = Transient {
                init_cond: DVector::from_vec(vec![1.0, 0.0]).div(100.0),
                response: vec![response_t5s, response_t10s],
            };

            let _2d: Response = Response {
                modal: modal._2d,
                free_transient: hand_calc_sol,
            };

            Self { _2d }
        }
    }

    #[test]
    fn et() {
        assert!(true)
    }
}
