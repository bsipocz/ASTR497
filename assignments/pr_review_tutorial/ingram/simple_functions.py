def fibonacci(max):
    values = [0, 1]
    while values[0] + values[1] < max:
        values.append(values[0] + values[1])
    return values


def factorial(value):
    if value == 0:
        return 1
    else:
        return value * factorial(value - 1)
