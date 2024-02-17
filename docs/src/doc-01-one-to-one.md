# Randomization, targeting _equal_  `1:1` allocation

The following settings are assumed:

- There are two _treatment arms_ investigated in a trial:
    + ``E`` is an _experimental_ treatment arm.
    + ``C`` is a _control_ treatment arm.
- ``n`` is a total _sample size_.
- ``N_1(j)`` and ``N_2(j)`` are _treatment numbers_, i.e., _sample sizes_ on treatments after the ``j^\text{th}`` allocation step (``N_1(j)+N_2(j) = j, \: j = 1, \ldots, n``).
- ``\delta_j`` is a _treatment indicator_:

```math
\delta_j = \left\{\begin{array}{rl}
1, &\text{if treatment }E \\
0, &\text{if treatment }C
\end{array}\right. .
```

Under these assumptions, a _restricted randomization_ procedure is defined as

```math
\begin{aligned}
\phi_1 &= \Pr(\delta_1 = 1) = 0.5; \\
\phi_j &= \Pr(\delta_j = 1|\delta_1, \ldots, \delta_{j-1}), \: j = 2, \ldots, n.
\end{aligned}
```