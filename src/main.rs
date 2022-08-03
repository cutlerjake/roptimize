mod constraint;
// mod environment;
mod affine_expr;
mod model;
mod simplex;
mod solver;
mod tableau;
mod var;

use var::{Environment, VarType, Variable, VariableCollection, VariableDefinition};

use crate::{
    affine_expr::AffineExpression,
    constraint::{Comp, Constraint},
    model::{Model, OptDir},
    simplex::Simplex,
    tableau::Tableau,
    //solver::Solver,
};

use ndarray::{array, Array2};

fn main() {
    let mut coll = VariableCollection::new();

    let mut env = Environment::new();

    let vd1 = VariableDefinition::new(VarType::Float)
        .with_lb(0)
        .with_name(String::from("a"));

    let vd2 = VariableDefinition::new(VarType::Float)
        .with_lb(0.0)
        .with_name(String::from("b"));


    let var1 = Variable::new(&mut env, vd1);
    let var2 = Variable::new(&mut env, vd2);

    let exp1 = &var1 + &var2;
    let exp2 = 2* &var1 + &var2;


    let obj_exp = 40* &var1 - 30* &var2;

    let cons1 = Constraint::new(exp1.clone(), Comp::Le, 12.0_f64);
    let cons2 = Constraint::new(exp2, Comp::Le, 16.0_f64);
    let cons3 = Constraint::new(exp1.clone(), Comp::Ge, 4.0_f64);
    let mut model = Model::new(env.clone());

    model.set_obj_fn(OptDir::Max, obj_exp);
    model.add_constraint(cons1);
    model.add_constraint(cons2);
    model.add_constraint(cons3);

    println!("model:\n{}", model);

    println!("standardized model:\n{}", model.as_standard_form(true).mdl);
    // let std = model.as_standard_form(true);

    // println!("{}", std.mdl);

    // let arr:Array2<f64> = array![
    //     [2.0, 8.0, 4.0, 2.0],
    //     [2.0, 5.0, 1.0, 5.0],
    //     [4.0, 10.0, -1.0, 1.0]
    // ];

    let mut smplx = Simplex::new();
    smplx.solve(&model, true);
    println!("Model solved!");

    // let mut tbl = Tableau::new(arr);
    
    // println!("{:?}", tbl.tbl());

    // tbl.rref();
    // println!("");

    // println!("{:?}", tbl.tbl());

}
