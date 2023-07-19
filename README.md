# roptimize
This is a very simple implementation of the SIMPLEX algorithm I implemented to improve my understanding. The library can handle modeling and solving continuous linear programs.

# Usage

```rust
    // Create a new environment
    let mut env = Environment::new();

    // Create a new variable definition with a lower bound of 0 and a name of "a"
    let vd1 = VariableDefinition::new(VarType::Float)
        .with_lb(0)
        .with_name(String::from("a"));

    // Create a new variable definition with a lower bound of 4 and a name of "b"
    let vd2 = VariableDefinition::new(VarType::Float)
        .with_lb(4.0)
        .with_name(String::from("b"));

    // Create new variables with the variable definitions vd1 and vd2
    let var1 = Variable::new(&mut env, vd1);
    let var2 = Variable::new(&mut env, vd2);

    // Create an expression with the variables var1 and var2
    let exp1 = &var1 + &var2;
    let exp2 = 2 * &var1 + &var2;

    // Create an expression with the variables var1 and var2
    let obj_exp = 40 * &var1 + 30 * &var2;

    // Create constraints
    let cons1 = Constraint::new(exp1.clone(), Comp::Le, 12.0_f64);
    let cons2 = Constraint::new(exp2, Comp::Le, 16.0_f64);
    let cons3 = Constraint::new(exp1, Comp::Ge, 4.0_f64);

    // Create a new model
    let mut model = Model::new(env);

    // Set objective function and add constraints
    model.set_obj_fn(OptDir::Max, obj_exp);
    model.add_constraint(cons1);
    model.add_constraint(cons2);
    model.add_constraint(cons3);

    // Create simplex solver and solve the model
    let mut smplx = Simplex::new();
    smplx.solve(&model, true);
    println!("Model solved!");
```

