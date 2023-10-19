# rtqf
The form x^2 + y^2 + 10*z^2 (x, y, z integers) is known as [Ramanujan's ternary quadratic form](https://en.wikipedia.org/wiki/Ramanujan%27s_ternary_quadratic_form) after appearing in a 1916 paper of Ramanujan. There is a known formula producing all even natural numbers not representable in this form, but there are only 18 known odd natural numbers not representable:
- 3, 7, 21, 31, 33, 43, 67, 79, 87, 133, 217, 219, 223, 253, 307, 391, 679, 2719

In 1997 K. Ono and K. Soundararajan proved that assuming the [generalized Riemann hypothesis](https://en.wikipedia.org/wiki/Generalized_Riemann_hypothesis) (GRH) that this list is complete, i.e. every odd integer above 2719 is representable in this form, at the time this was verified for all odd integers between 2719 and 2e+10, I have currently verified this up to 2.5e+10. rtqf.c searches for a counterexample/a non-representable odd integer above 2e+10), if one is found this would disprove GRH.
