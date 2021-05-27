

int_fact = 1000000000.0
emp_grad_dval = 0.0001

assert int_fact * emp_grad_dval > 1.0

def integerize_single(val):
    return int(val*int_fact)

def integerize(iso):
    for mass, prob in zip(iso.masses, iso.probs):
        yield (integerize_single(mass), integerize_single(prob))


flows_loud = True
