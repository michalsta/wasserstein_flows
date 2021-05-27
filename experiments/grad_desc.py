import math
from collections import deque
from flows_cl import AbWSDistCalc
import numpy as np


def fix(point):
    return np.array([max(x, 0.00000001) for x in point])

def gradient_descent(start_point, val_f, grad_f):
    point = np.array(start_point)
    annealing_fact = 1.0
    changes_dir = grad_f(point)

    minima = deque()
    best_seen = val_f(point)
    best_point = point
    minima.append(point)

    for _ in range(100000):
        annealing_fact *= 0.999
        point -= changes_dir * annealing_fact
        point = fix(point)
        current_val = val_f(point)
        if current_val < best_seen:
            best_seen = current_val
            best_point = point
        minima.append(current_val)
        if len(minima) >= 100:
            minima.popleft()
            if max(minima) - 0.01 < best_seen:
                return point
        print(point)
        changes_dir = changes_dir * 0.9 + grad_f(point)

if __name__ == "__main__":

    EXP = [
            ((0.0, 0.0), 2.0),
            ((1.0, 2.0), 4.0),
          ]

    THE1 = [((0.0, 0.0), 1.0)]
    THE2 = [((1.0, 2.0), 1.0)]

    def euclidean_distance(point1, point2):
        return math.sqrt(sum((x_i - y_i)*(x_i - y_i) for x_i, y_i in zip(point1, point2)))

    wasserstein_distancer = AbWSDistCalc(EXP, [THE1, THE2], 5.0, 5.0, euclidean_distance)

    print("------------------------")

    gradient_descent([3245.0, 3563.0], wasserstein_distancer.value_at, wasserstein_distancer.gradient_at)
#    print(scipy.optimize.minimize(wasserstein_distancer.value_at, [12.0, 15.0], jac=wasserstein_distancer.gradient_at, method='SLSQP'))

