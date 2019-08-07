a = 1
b = 1/400
c = 0.5
x = 665.1127

x.asinh <- asinh(a + b * x) + c

(sinh(x.asinh - c) - a) * 1/b


