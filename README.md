# roptimize
This is a very simple implementation of the SIMPLEX algorithm I implemented to improve my understanding. The library can handle modeling and solving continuous linear programs.

# Usage

```rust
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

    let mut smplx = Simplex::new();
    smplx.solve(&model, true);
    println!("Model solved!");
```

