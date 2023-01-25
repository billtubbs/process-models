% test bounded_angle.m

assert(bounded_angle(0) == 0)
assert(bounded_angle(pi) == pi)
assert(bounded_angle(-pi) == -pi)

assert(rad2deg(bounded_angle(1.25*pi)) == -135)
assert(rad2deg(bounded_angle(-1.25*pi)) == 135)

assert(rad2deg(bounded_angle(1.75*pi)) == -45)
assert(rad2deg(bounded_angle(-1.75*pi)) == 45)

assert(abs(rad2deg(bounded_angle(3.25*pi)) - (-135)) < 1e-13)
assert(abs(rad2deg(bounded_angle(-3.25*pi)) - 135) < 1e-13)

assert(abs(rad2deg(bounded_angle(3.75*pi)) - (-45)) < 1e-13)
assert(abs(rad2deg(bounded_angle(-3.75*pi)) - 45) < 1e-13)

assert(abs(rad2deg(bounded_angle(5.25*pi)) - (-135)) < 1e-13)
assert(abs(rad2deg(bounded_angle(-5.25*pi)) - 135) < 1e-13)