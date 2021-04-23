def fibonacci(max):
    values = [0, 1]
    iter = 0
    while values[0] + values[1] < max:
        values.append(values[iter] + values[iter+=1])
    return values


def factorial(value):
    if value == 0:
        return 1
    else:
        return value * factorial(value - 1)
