#[allow(unused_variables, dead_code)]
#[cfg_attr(feature = "cargo-clippy", allow(dead_code, unused_variables))]
mod constraint;
// mod environment;
mod affine_expr;
mod model;
mod print_table;
mod simplex;
mod solver;
mod tableau;
mod var;

use var::{Environment, VarType, Variable, VariableCollection, VariableDefinition};

use crate::{
    constraint::{Comp, Constraint},
    model::{Model, OptDir},
    simplex::Simplex,
    tableau::{SparseTableau, Tableau, TableauIx},
    //solver::Solver,
};

fn main() {
    let _coll = VariableCollection::new();

    let mut env = Environment::new();

    let vd1 = VariableDefinition::new(VarType::Float)
        //.with_lb(0)
        .with_name(String::from("a"));

    let vd2 = VariableDefinition::new(VarType::Float)
        //.with_lb(4.0)
        .with_name(String::from("b"));

    let var1 = Variable::new(&mut env, vd1);
    let var2 = Variable::new(&mut env, vd2);

    let exp1 = &var1 + &var2;
    let exp2 = 2 * &var1 + &var2;

    let obj_exp = 40 * &var1 + 30 * &var2;

    let cons1 = Constraint::new(exp1.clone(), Comp::Le, 12.0_f64);
    let cons2 = Constraint::new(exp2, Comp::Le, 16.0_f64);
    let cons3 = Constraint::new(exp1, Comp::Ge, 4.0_f64);
    let mut model = Model::new(env.clone());

    model.set_obj_fn(OptDir::Max, obj_exp);
    model.add_constraint(cons1);
    model.add_constraint(cons2);
    model.add_constraint(cons3);

    println!("model:\n{}", model);

    let std_mdl = model.as_standard_form(true);

    let dt = std_mdl.mdl.as_tableau();

    let mut dense = std_mdl.mdl.as_tableau();
    let mut sp = SparseTableau::from(dense.clone());

    println!("standardized model:\n{}", std_mdl.mdl);
    println!("{}", dense);

    // for _ in 0..2 {
    //     println!("______________________________________________________");
    //     println!("Dense model:\n{}\n", dense);
    //     println!("Sparse:");
    //     println!(
    //         "\tConstraints:\n{}",
    //         format!("\t{}", sp.constraints.to_dense()).replace("\n", "\n\t")
    //     );
    //     println!(
    //         "\tBasis ({:?}):\n{}\n",
    //         sp.basis.shape(),
    //         format!("\t{}", sp.basis.to_dense()).replace("\n", "\n\t")
    //     );
    //     println!("\tRelative costs:\n\t{:?}", sp.r_costs().to_dense());
    //     println!(
    //         "\tSP to Dense: \n{}",
    //         format!("\t{}", Tableau::from(&sp)).replace("\n", "\n\t")
    //     );

    //     let (j, _) = sp
    //         .obj_fn
    //         .iter()
    //         .max_by(|(i1, val1), (i2, val2)| val1.partial_cmp(val2).unwrap())
    //         .unwrap();
    //     let i = 2;
    //     let ix = TableauIx::new(i, j);

    //     println!("Pivot ind: {:?}", ix);

    //     sp.pivot(ix);
    //     dense.pivot(&ix);
    // }

    // println!("Tableau: \n{}\n", sp.constraints.to_dense());
    // println!(
    //     "\tOriginal Basis:\n{}\n",
    //     format!("\t{}", sp.basis.to_dense()).replace("\n", "\n\t")
    // );
    // println!("Relative costs:\n\t{:?}", sp.r_costs().to_dense());

    // println!("SP to Dense: \n{}", Tableau::from(&sp));

    // let (j, _) = sp
    //     .obj_fn
    //     .iter()
    //     .max_by(|(i1, val1), (i2, val2)| val1.partial_cmp(val2).unwrap())
    //     .unwrap();
    // let i = 2;
    // let ix = TableauIx::new(i, j);

    // println!("Pivot ind: {:?}", ix);

    // sp.pivot(ix);

    // println!("Tableau: \n{}\n", sp.constraints.to_dense());
    // println!(
    //     "\tPivot 1 Basis:\n{}\n",
    //     format!("\t{}", sp.basis.to_dense()).replace("\n", "\n\t")
    // );
    // println!(
    //     "Reletive costs:\n\t{:?}",
    //     sp.nb_r_costs().collect::<Vec<(usize, f64)>>()
    // );

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
