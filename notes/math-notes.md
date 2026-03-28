# Hormetic-and-hydra-effects-in-a-discrete-predator-prey-model
Appendix: Supplementary Information for Review Only
\appendix{Proof of the Flip bifurcation}
\label{TheoremA}
\noindent
\begin{proof}
Take $u=x-x_{3}$, $v=y-y_{3}$, $\gamma = q_2 -{q_2}^{*}$ and transform the unique positive fixed point $E_{3}=(x_{3}, y_{3})$ to the origin, $\gamma$ is a new dependent variable. The system (\ref{map_out}) is transformed into the following form
\begin{equation}
\left(\begin{array}{l}
u \\
v \\
\gamma
\end{array}\right) \rightarrow\left(\begin{array}{l}
a_{100} u + a_{010} v + a_{200} u^2 + a_{110} u v + a_{020} v^2 + a_{300} u^3 + a_{210} u^2 v + a_{120} u v^2 \\ + a_{030} v^3 + \mathcal{O}\left(|u, v, \gamma|^4\right) \\
b_{100} u + b_{010} v + b_{001} \gamma + b_{200} u^2 + b_{110} u v + b_{101} u \gamma + b_{011} v \gamma + b_{300} u^3 \\ + b_{210} u^2 v  +
b_{201} u^2 \gamma + b_{111} u v \gamma + \mathcal{O}\left(|u, v, \gamma|^4\right) \\
\gamma
\end{array}\right),
\label{Flip_map}
\end{equation}
where $a_{100}=\frac{K_1- r_1 x_3}{K_1}$, $a_{010}=-\beta x_3$, $a_{200}=\frac{r_1 (r_1 x_3-2 K_1)}{2 {K_1}^2}$, $a_{110}=-\frac{\beta (K_1-r_1 x_3)}{K_1}$, $a_{020}=\frac{\beta^2 x_3}{2}$, $a_{300}=\frac{{r_1}^2 (3 K_1-r_1 x_3)}{6 {K_1}^3}$, $a_{210}=\frac{\beta r_1 (2 K_1-r_1 x_3)}{2 {K_1}^2}$, $a_{120}=\frac{{\beta}^2 (K_1-r_1 x_3)}{2 K_1}$, $a_{030}=\frac{{\beta}^3 x_3}{6}$, $b_{100}=\beta \eta y_3$, $b_{010}=\frac{K_2-r_2 y_3}{K_2}$, $b_{001}=\frac{y_3}{1-\gamma}$, $b_{200}=\frac{y_3 \beta^2 \eta^2}{2}$, $b_{110}=\frac{\beta \eta (K_2-r_2 y_3)}{K_2}$, $b_{101}=-\frac{y_3 \beta \eta}{1-\gamma}$, $b_{020}=\frac{r_2 (r_2 y_3 - K_2)}{2 {K_2}^2}$, $b_{011}=-\frac{(r_2 y_3 - K_2)}{K_2 (1-\gamma)}$, $b_{300}=\frac{y_3 {\beta}^3 {\eta}^3}{6}$,
$b_{210}=\frac{{\beta}^2 {\eta}^2 (K_2-r_2 y_3)}{2 K_2}$, $b_{201}=-\frac{y_3 {\beta}^2 {\eta}^2}{2 (1-\gamma)}$, $b_{120}=\frac{r_2 \beta \eta (r_2 y_3 - 2 K_2)}{2 {K_2}^2}$, $b_{111}=\frac{\beta \eta(r_2 y_3 - K_2)}{K_2 (1-\gamma)}$, $b_{030}=\frac{{r_2}^2 (3 K_2 - r_2 y_3)}{6 {K_2}^3}$,
$b_{021}=\frac{r_2 (2 K_2 - r_2 y_3)}{2 (1-\gamma) {K_2}^2}$.

We construct an invertible matrix
\begin{equation*}
T_{M}=\left[\begin{array}{ccc}
1 & 1 & 1 \\
X_1 &X_2 &X_3 \\
0 & 0 &X_4
\end{array}\right]
\end{equation*}
with $X_1=\frac{b_{010}-a_{100}+\sqrt{\Delta}}{2 a_{010}}$, $X_2=\frac{b_{010}-a_{100}-\sqrt{\Delta}}{2 a_{010}}$,  $X_3=\frac{1-a_{100}}{a_{010}}$, $X_4=\frac{1+a_{100} b_{010}- a_{010} b_{100}-a_{100}+b_{010}}{a_{010} b_{001}}$, where $\Delta=4 a_{010} b_{100}+{a_{100}}^2+{b_{010}}^2-2 a_{100} b_{010}$.

Use the transformation
\begin{equation*}
\left(\begin{array}{l}
u \\
v \\
\gamma
\end{array}\right)=T_{M}\left(\begin{array}{l}
X \\
Y \\
Z \\
\end{array}\right),
\end{equation*}
then the system (\ref{Flip_map}) can be written as
\begin{equation*}
\left(\begin{array}{c}
X \\
Y \\
Z
\end{array}\right) \rightarrow\left(\begin{array}{ccc}
A_{100} & 0 & 0 \\
0 & -1 & 0  \\
0 & 0 &  1
\end{array}\right)\left(\begin{array}{c}
X \\
Y \\
Z
\end{array}\right)+\left(\begin{array}{c}
l(X, Y, Z) \\
m(X, Y, Z) \\
0
\end{array}\right),
\end{equation*}
where
\begin{equation*}
\begin{aligned}
l(X, Y, Z)=&A_{200}X^{2}+A_{110}X Y+A_{101}X Z+A_{020} Y^{2}+A_{011} Y Z+A_{002} Z^{2}+A_{300}X^{3}\\ &+A_{210}X^{2} Y
+A_{201}X^{2} Z+A_{120}X Y^{2}+A_{111}X Y Z+A_{102}X Z^{2}+A_{030} Y^{3}\\ &+A_{021} Y^{2} Z+A_{012} Y Z^{2}+A_{003} Z^{3}
\end{aligned}
\end{equation*}
and
\begin{equation*}
\begin{aligned}
m(X, Y, Z)=& B_{200}X^{2}+B_{110}X Y+B_{101}X Z+B_{020} Y^{2}+B_{011} Y Z+B_{002} Z^{2}+B_{300}X^{3}\\ & +B_{210}X^{2} Y
+B_{201} X^{2} Z+B_{120}X Y^{2}+B_{111}X Y Z+B_{102}X Z^{2}+B_{030} Y^{3}\\ & +B_{021} Y^{2} Z+B_{012} Y Z^{2}+B_{003} Z^{3},
\end{aligned}
\end{equation*}
where the coefficients of $l(X, Y, Z)$ and $m(X, Y, Z)$ are listed as follows:\\
$A_{100}=\frac{1}{2 a_{010}}(a_{100}+b_{010}+\sqrt{4 a_{010} b_{100}+{a_{100}}^2+{b_{010}}^2-2 a_{100} b_{010}})$, \\ $A_{200}=\frac{1}{X_2 -X_1}({X_{1}}^2X_{2} a_{020} - {X_{1}}^2 b_{020} +X_{1}X_{2} a_{110} -X_{1} b_{110} +X_{2} a_{200} - b_{200})$,\\
$A_{110}=\frac{2X_{1} {X_{2}}^2 a_{020} +X_{1}X_{2} a_{110} - 2X_{1}X_{2} b_{020} + {X_{2}}^2 a_{110} -X_{1} b_{110} + 2X_{2} a_{200} -X_{2} b_{110} - 2 b_{200}}{X_2 -X_1}$,\\
$A_{101}=\frac{1}{X_2 -X_1}(2X_{1}X_{2}X_{3} a_{020} +X_{1}X_{2} a_{110} - 2X_{1}X_{3} b_{020} -X_{1}X_{4} b_{011} +X_{2}X_{3} a_{110} -X_{1} b_{110} + 2X_{2} a_{200} -X_{3} b_{110} -X_{4} b_{101} - 2 b_{200})$,\\
$A_{020}=\frac{{X_{2}}^3 a_{020} + {X_{2}}^2 a_{110} - {X_{2}}^2 b_{020} +X_{2} a_{200} -X_{2} b_{110} - b_{200}}{X_2 -X_1}$,\\
$A_{011}=\frac{1}{X_2 -X_1}(2 {X_{2}}^2X_{3} a_{020} + {X_{2}}^2 a_{110} +X_{2}X_{3} a_{110} - 2X_{2}X_{3} b_{020} -X_{2}X_{4} b_{011} + 2X_{2} a_{200} -X_{2} b_{110} -X_{3} b_{110} -X_{4} b_{101} - 2 b_{200})$,\\
$A_{002}=\frac{X_{2} {X_{3}}^2  a_{020} + X_{2}X_{3} a_{110} - {X_{3}}^2 b_{020} -X_{3}X_{4} b_{011} +X_{2} a_{200} -X_{3} b_{110} -X_{4} b_{101} - b_{200}}{X_2 -X_1}$,\\
$A_{300}=\frac{{X_{1}}^3 X_{2} a_{030} - {X_{1}}^3 b_{030} + {X_{1}}^2 X_{2} a_{120} - {X_{1}}^2 b_{120} + X_{1} X_{2} a_{210} - X_{1} b_{210} + X_{2} a_{300} - b_{300}}{X_2 -X_1}$,\\
$A_{210}=\frac{1}{X_2 -X_1}(3 {X_{1}}^2 {X_{2}}^2 a_{030} + {X_{1}}^2 X_{2} a_{120} - 3 {X_{1}}^2 X_{2} b_{030} + 2 X_{1} {X_{2}}^2 a_{120} - {X_{1}}^2 b_{120} + 2 X_{1} X_{2} a_{210} - 2 X_{1} X_{2} b_{120} + {X_{2}}^2 a_{210} - 2 {X_{1}} b_{210} + 3 X_{2} a_{300} - X_{2} b_{210} - 3 b_{300})$,\\
$A_{201}=\frac{1}{X_2 -X_1} ({3 X_{1}}^2 X_{2} X_{3} a_{030} + {X_{1}}^2 X_{2} a_{120} - 3 {X_{1}}^2 X_{3} b_{030} - X_{1}^2 X_{4} b_{021} + 2 X_{1} X_{2} X_{3} a_{120} - {X_{1}}^2 b_{120} + 2 X_{1} X_{2} a_{210} - 2 X_{1} X_{3} b_{120} - X_{1} X_{4} b_{111} + X_{2} X_{3} a_{210} - 2 X_{1} b_{210} + 3 X_{2} a_{300} - X_{3} b_{210} - X_{4} b_{201} - 3 b_{300})$,\\
$A_{120}=\frac{1}{X_2 -X_1}(3 X_{1} {X_{2}}^3 a_{030} + 2 X_{1} {X_{2}}^2 a_{120} - 3 X_{1} {X_{2}}^2 b_{030} + X_{2}^3 a_{120} + X_{1} X_{2} a_{210} - 2 X_{1} X_{2} b_{120} + 2 X_{2}^2 a_{210} - {X_{2}}^2 b_{120} - X_{1} b_{210}+ 3 X_{2} a_{300} - 2 X_{2} b_{210} - 3 b_{300})$,\\
$A_{111}=\frac{1}{X_2 -X_1} (6 X_{1} {X_{2}}^2 X_{3} a_{030} + 2 X_{1} {X_{2}}^2 a_{120} + 2 X_{1} X_{2} X_{3} a_{120} - 6 X_{1} X_{2} X_{3} b_{030} - 2 X_{1} X_{2} X_{4} b_{021} + 2 {X_{2}}^2 X_{3} a_{120} + 2 X_{1} X_{2} a_{210} - 2 X_{1} X_{2} b_{120} - 2 X_{1} X_{3} b_{120} - X_{1} X_{4} b_{111} + 2 {X_{2}}^2 a_{210} + 2 X_{2} X_{3} a_{210} - 2 X_{2} X_{3} b_{120} - X_{2} X_{4} b_{111} - 2 X_{1} b_{210} + 6 X_{2} a_{300} - 2 X_{2} b_{210} - 2 X_{3} b_{210} - 2 X_{4} b_{201} - 6 b_{300})$,\\
$A_{102}=\frac{1}{X_2 -X_1} (3 X_{1} X_{2} {X_{3}}^2 a_{030} + 2 X_{1} X_{2} X_{3} a_{120} - 3 X_{1} {X_{3}}^2 b_{030} - 2 X_{1} X_{3} X_{4} b_{021} + X_{2} {X_{3}}^2 a_{120} + X_{1} X_{2} a_{210} - 2 X_{1} X_{3} b_{120} - X_{1} X_{4} b_{111} + 2 X_{2} X_{3} a_{210} - {X_{3}}^2 b_{120} - X_{3} X_{4} b_{111} - X_{1} b_{210} + 3 X_{2} a_{300} - 2 X_{3} b_{210} - 2 X_{4} b_{201} - 3 b_{300})$,\\
$A_{030}=\frac{{X_{2}}^4 a_{030} + {X_{2}}^3 a_{120} - {X_{2}}^3 b_{030} + {X_{2}}^2 a_{210} - {X_{2}}^2 b_{120} + X_{2} a_{300} - X_{2} b_{210} - b_{300}}{X_2 -X_1}$,\\
$A_{021}=\frac{1}{X_2 -X_1} (3 {X_{2}}^3 X_{3} a_{030} + {X_{2}}^3 a_{120} + 2 {X_{2}}^2 X_{3} a_{120} - 3 {X_{2}}^2 X_{3} b_{030} - {X_{2}}^2 X_{4} b_{021} + 2 {X_{2}}^2 a_{210} - {X_{2}}^2 b_{120} + X_{2} X_{3} a_{210}- 2 X_{2} X_{3} b_{120} - X_{2} X_{4} b_{111} + 3 X_{2} a_{300} - 2 X_{2} b_{210} - X_{3} b_{210} - X_{4} b_{201} - 3 b_{300})$,\\
$A_{012}=\frac{1}{X_2 -X_1} (3 {X_{2}}^2 {X_{3}}^2 a_{030} + 2 {X_{2}}^2 X_{3} a_{120} + X_{2} {X_{3}}^2 a_{120} - 3 X_{2} {X_{3}}^2 b_{030} - 2 X_{2} X_{3} X_{4} b_{021} + {X_{2}}^2 a_{210} + 2 X_{2} X_{3} a_{210} - 2 X_{2} X_{3} b_{120} - X_{2} X_{4} b_{111} - {X_{3}}^2 b_{120} - X_{3} X_{4} b_{111} + 3 X_{2} a_{300} - X_{2} b_{210} - 2 X_{3} b_{210} - 2 X_{4} b_{201} - 3 b_{300})$,\\
$A_{003}=\frac{1}{X_2 -X_1}(X_{2} {X_{3}}^3 a_{030} + X_{2} {X_{3}}^2 a_{120} - {X_{3}}^3 b_{030} - {X_{3}}^2 X_{4} b_{021} + X_{2} X_{3} a_{210} - {X_{3}}^2 b_{120} - X_{3} X_{4} b_{111} + X_{2} a_{300} - X_{3} b_{210} - X_{4} b_{201} - b_{300})$,\\																											
$B_{200}=\frac{{X_{1}}^3 a_{020} + {X_{1}}^2 a_{110} - {X_{1}}^2 b_{020} +X_{1}  a_{200} -X_{1} b_{110} - b_{200}}{X_1 -X_2}$,\\
$B_{110}=\frac{2 {X_{1}}^2X_{2} a_{020} + {X_{1}}^2 a_{110} +X_{1}X_{2} a_{110} - 2X_{1}X_{2} b_{020} + 2X_{1} a_{200} -X_{1} b_{110} -X_{2} b_{110}- 2 b_{200}}{X_1 -X_2}$,\\
$B_{101}=\frac{1}{X_1 -X_2}(2 {X_{1}}^2X_{3} a_{020} + {X_{1}}^2 a_{110} +X_{1}X_{3} a_{110} - 2X_{1}X_{3} b_{020} -X_{1}X_{4} b_{011} + 2X_{1} a_{200} -X_{1} b_{110}-X_{3} b_{110} -X_{4} b_{101} - 2 b_{200})$,\\
$B_{020}=\frac{X_{1} {X_{2}}^2 a_{020} +X_{1}X_{2} a_{110} - {X_{2}}^2 b_{020} +X_{1} a_{200} -X_{2} b_{110} - b_{200}}{X_1 -X_2}$,\\
$B_{011}=\frac{1}{X_1 -X_2}(2X_{1}X_{2}X_{3} a_{020} +X_{1}X_{2} a_{110} +X_{1}X_{3} a_{110} - 2X_{2}X_{3} b_{020} -X_{2}X_{4} b_{011} + 2X_{1} a_{200} -X_{2} b_{110} -X_{3} b_{110} -X_{4} b_{101} - 2 b_{200})$,\\
$B_{002}=\frac{X_{1} {X_{3}}^2 a_{020} +X_{1}X_{3} a_{110} - {X_{3}}^2 b_{020} -X_{3}X_{4} b_{011} +X_{1} a_{200} -X_{3} b_{110} -X_{4} b_{101} - b_{200}}{X_1 -X_2}$,\\
$B_{300}=\frac{X_{1}^4 a_{030} + {X_{1}}^3 a_{120} - {X_{1}}^3 b_{030} + {X_{1}}^2 a_{210} - {X_{1}}^2 b_{120} +X_{1} a_{300} -X_{1} b_{210} - b_{300}}{X_1 -X_2}$,\\
$B_{210}=\frac{1}{X_1 -X_2}(3 {X_{1}}^3X_{2} a_{030} + {X_{1}}^3 a_{120} + 2 {X_{1}}^2X_{2} a_{120} - 3 {X_{1}}^2X_{2} b_{030} + 2 {X_{1}}^2 a_{210} - {X_{1}}^2 b_{120} +X_{1}X_{2} a_{210} - 2X_{1}X_{2} b_{120} + 3X_{1} a_{300} - 2X_{1} b_{210} -X_{2} b_{210} - 3 b_{300})$,\\
$B_{201}=\frac{1}{X_1 -X_2}((3X_{3} a_{030} +  a_{120}){X_{1}}^3 + (2X_{3} a_{120} - 3X_{3} b_{030} -X_{4} b_{021} + 2 a_{210} - b_{120}) {X_{1}}^2 + (a_{210} - 2 b_{120})X_{1}X_{3} - (X_{4} b_{111} + 3 a_{300} - 2 b_{210})X_{1} -X_{3} b_{210} -X_{4} b_{201} - 3 b_{300})$,\\
$B_{120}=\frac{1}{X_1 -X_2}(3 {X_{1}}^2 {X_{2}}^2 a_{030} + 2 {X_{1}}^2X_{2} a_{120} +X_{1} {X_{2}}^2 a_{120} - 3X_{1} {X_{2}}^2 b_{030} + {X_{1}}^2 a_{210} + 2X_{1}X_{2} a_{210} - 2X_{1}X_{2} b_{120} - {X_{2}}^2 b_{120} + 3X_{1} a_{300} -X_{1} b_{210} - 2X_{2} b_{210} - 3 b_{300})$,
$B_{111}=\frac{1}{X_1 -X_2} ((X_{2}X_{3} a_{030} + 2 {X_{1}}^2X_{2} a_{120} + 2X_{3} a_{120} + 2  a_{210}){X_{1}}^2 + (2 a_{120} - 6 b_{030})X_{1}X_{2}X_{3} + (- 2X_{4} b_{021} + 2 a_{210} - 2 b_{120})X_{1}X_{2} + 2X_{1}X_{3} a_{210} - 2X_{1}X_{3} b_{120} -X_{1}X_{4} b_{111}
- 2X_{2}X_{3} b_{120} -X_{2}X_{4} b_{111} + 6X_{1} a_{300} - 2X_{1} b_{210} - 2X_{2} b_{210} - 2X_{3} b_{210} - 2X_{4} b_{201} - 6 b_{300})$,\\
$B_{102}=\frac{1}{X_1 -X_2} (3 {X_{1}}^2 {X_{3}}^2 a_{030} + 2 {X_{1}}^2 X_{3} a_{120} + X_{1} {X_{3}}^2 a_{120} - 3 X_{1} {X_{3}}^2 b_{030} - 2 X_{1} X_{3} X_{4} b_{021} + {X_{1}}^2 a_{210} + 2 X_{1} X_{3} a_{210} - 2 X_{1} X_{3} b_{120} - X_{1} X_{4} b_{111} - {X_{3}}^2 b_{120} - X_{3} X_{4} b_{111} + 3 X_{1} a_{300} - X_{1} b_{210} - 2 X_{3} b_{210} - 2 X_{4} b_{201} - 3 b_{300})$,\\
$B_{030}=\frac{X_{1} {X_{2}}^3 a_{030} +X_{1} {X_{2}}^2 a_{120} - {X_{2}}^3 b_{030} +X_{1}X_{2} a_{210} - {X_{2}}^2 b_{120} +X_{1} a_{300} -X_{2} b_{210}- b_{300}}{X_1 -X_2}$,\\
$B_{021}=\frac{1}{X_1 -X_2} (3X_{1} {X_{2}}^2X_{3} a_{030} +X_{1} {X_{2}}^2 a_{120} + 2X_{1}X_{2}X_{3} a_{120} - 3 {X_{2}}^2X_{3} b_{030} - {X_{2}}^2X_{4} b_{021} + 2X_{1}X_{2} a_{210} +X_{1}X_{3} a_{210} - {X_{2}}^2 b_{120} - 2X_{2}X_{3} b_{120} -X_{2}X_{4} b_{111} + 3X_{1} a_{300} - 2X_{2} b_{210} -X_{3} b_{210} -X_{4} b_{201} - 3 b_{300})$,\\
$B_{012}=\frac{1}{X_1 -X_2} (3X_{1} X_{2} {X_{3}}^2 a_{030} + 2X_{1} X_{2} X_{3} a_{120} +X_{1} {X_{3}}^2 a_{120} - 3 X_{2} {X_{3}}^2 b_{030} - 2 X_{2} X_{3} X_{4} b_{021} +X_{1} X_{2} a_{210} + 2 X_{1} X_{3} a_{210} - 2X_{2} X_{3} b_{120} -X_{2} X_{4} b_{111} -{X_{3}}^2 b_{120} -X_{3} X_{4} b_{111} + 3 X_{1} a_{300} -X_{2} b_{210} - 2 X_{3} b_{210} - 2X_{4} b_{201} - 3 b_{300})$,\\
$B_{003}=\frac{1}{X_1 -X_2}(X_{1} {X_{3}}^3 a_{030} +X_{1} {X_{3}}^2 a_{120} -{X_{3}}^3 b_{030} -{X_{3}}^2 X_{4} b_{021} +X_{1} X_{3} a_{210} -{X_{3}}^2 b_{120} -X_{3} X_{4} b_{111} +X_{1} a_{300} -X_{3} b_{210} -X_{4} b_{201} - b_{300})$.

By the center manifold theory \cite{guckenheimer2013nonlinear,layek2015introduction}, the center manifold can be represented as follows
\begin{equation*}
W^c(0) =\left\{(X, Y, Z) \in R^3 \mid X=H(Y, Z), H(0,0)=0, D H(0,0)=0\right\}
\end{equation*}
for $Y$ and $Z$ sufficiently small and we have
\begin{equation}
H(Y, Z)=h_1 Y^2+h_2 Y Z+h_3 Z^2+\mathcal{O}\left(|Y, Z|^3\right).
\label{heq3}
\end{equation}

One yields
\begin{equation}
\mathcal{N}(H(Y, Z))\!=H\Big(\!-Y \!+ m(H(Y, Z),Y,Z),Z\Big)\!-A_{100} H(Y, Z)\!-l(H(Y, Z),Y,Z)\!=0.
\label{heq4}
\end{equation}

Substituting (\ref{heq3}) into (\ref{heq4}) and comparing the coefficients of (\ref{heq4}), we obtain
$$
h_1= \frac{A_{020}}{1-A_{100}}, h_2=\frac{A_{011}}{-A_{100} -1 }, h_3=\frac{A_{002}}{1-A_{100}}.
$$
Thus, the system restricted to the center manifold is given by
\begin{equation}
\begin{aligned}
Y \rightarrow \tilde{n}(Y, Z) = & -Y + C_{011} Y Z + C_{020} Y^2 + C_{002} Z^2+ C_{030} Y^3+ C_{003} Z^3  \\ &+C_{021} Y^2 Z+C_{012} Y Z^2 + \mathcal{O}\left(|Y, Z|^4\right),
\end{aligned}
\label{flip_cen}
\end{equation}
where $C_{011}=B_{011}$, $C_{020}=B_{020}$, $C_{002}=B_{002}$, $C_{030}=B_{110} h_1 + B_{030}$, $C_{003}=B_{101} h_3 + B_{003}$, $C_{021}=B_{110} h_2+ B_{101} h_1 + B_{021}$, $C_{012}=B_{110} h_3+ B_{101} h_2 + B_{012}$.

We can see that $\tilde{n}(0, 0)=0$, $\frac{\partial \tilde{n}}{\partial Y}(0,0)=-1$,
$\frac{\partial^2 \tilde{n}}{\partial Z \partial Y}(0,0)=C_{011}  \neq 0$, and $\frac{1}{2}\left(\frac{\partial^2 \tilde{n}}{\partial Y^2}(0,0)\right)^2 + \frac{1}{3} \frac{\partial^3 \tilde{n}}{\partial Y^3}(0,0)= 2 C_{020}^2 + 2 C_{030} \neq 0$. According to \citet{kuznetsov1998elements}, the fixed point $E_{3}=(x_3, y_{3})$ undergoes a Flip bifurcation at $q_2={q_2}^*$ when the conditions of Theorem~\ref{flip_theorem} are satisfied. The proof is completed.
\end{proof}

\appendix{Proof of the Neimark-Sacker bifurcation}
\label{TheoremB}
\begin{proof}

From the results in Section~\ref{sec:existence}, the fixed point $E_{3} = (x_{3}, y_{3})$ of system~(\ref{map_out}) is non-hyperbolic when the conditions of Theorem~\ref{ns_theorem} are satisfied, with the eigenvalues of the Jacobian matrix $J_{3}$ given by $\lambda = \overline{\lambda}$ and $|\lambda| = |\overline{\lambda}| = 1$. The characteristic equation of the Jacobian is:
$$
\lambda^2 - m(q_2) \lambda +n(q_2)=0,
$$
where $m(q_2)=a_{100} + b_{010}=\frac{2 A - (r_1 B + r_2 C)}{A}$ and $n(q_2)=a_{100} b_{010}-a_{010} b_{100}=\frac{A + B C- (r_1 B + r_2 C)}{A}$.

Thus, the eigenvalues $\lambda, \bar{\lambda}$ are a pair of conjugate complexes with module 1 and satisfies
$$
\lambda, \bar{\lambda}=\frac{m(q_2)}{2} \pm \frac{\sqrt{4 n(q_2)-m^2(q_2)}}{2} i,
$$
and consequently we have the following transversality condition when $q_2$ is sufficiently closing to ${q_2}^*$,
\begin{equation*}
|\lambda, \bar{\lambda}|=\sqrt{n(q_2)}, \left.\frac{d|\lambda, \bar{\lambda}|}{d q_2}\right|_{q_2={q_2}^*} = \frac{B^* {r_{1}}^2 + (C^*-r_1) r_1 K_2 \beta + r_1 r_2}{A (1-{q_2}^*)} \neq 0 .
\end{equation*}

In addition, ${\lambda}^j, {\bar{\lambda}}^j \neq 1, j = 1,2,3,4$ must hold, which implies that $m^2({q_2}^*) - 4 n({q_2}^*) < 0$ with $n({q_2}^*)=1$ and $m({q_2}^*) \neq 0, -1$, i.e.
$$
\left\{\begin{array}{l}
	B^* C^* < 4 A,\\
	2 A - B^* C^* \neq 0,\\
	A - B^* C^* \neq 0.
\end{array}\right.
$$

We now compute the non-degeneracy condition for the NS bifurcation. Building on the previous analysis, we choose the natural enemy's killing rate $q_2$ as the bifurcation parameter and derive the conditions for the existence of NS bifurcation using bifurcation theory \cite{tang2002chaos,guckenheimer2013nonlinear}. Let $u = x - x_3$ and $v = y - y_3$, and transform the fixed point $E_3 = (x_3, y_3)$ to the origin. Expanding the right-hand side of system~(\ref{map_out}) around the origin, we obtain:
\begin{equation}
	\left(\begin{array}{l}
		u \\
		v \\
	\end{array}\right) \rightarrow\left(\begin{array}{l}
		a_{100} u + a_{010} v + a_{200} u^2 + a_{110} u v + a_{020} v^2 + a_{300} u^3 + a_{210} u^2 v \\+ a_{120} u v^2 + a_{030} v^3 +  \mathcal{O}\left(|u, v|^4\right) \\
		b_{100} u + b_{010} v + b_{200} u^2 + b_{110} u v + b_{300} u^3 + b_{210} u^2 v + \mathcal{O}\left(|u, v|^4\right) \\
	\end{array}\right),
	\label{Ns_Map}
\end{equation}
where the coefficients for each item are same as those in system (\ref{Flip_map}) by replacing $q^{*}$ by ${q_2}^*$ and $\gamma$ by $0$. Make the following transformation:
\begin{equation}
\left(\begin{array}{l}
u \\
v \\
\end{array}\right)=T_{M}\left(\begin{array}{l}
X \\
Y \\
\end{array}\right)
\end{equation}
with
\begin{equation}
	T_{M}=\left[\begin{array}{cc}
		1 & 0 \\
		\frac{b_{010}-a_{100}}{2 a_{010}} & \frac{\sqrt{4 n(q_2)-m^2(q_2)}}{2 a_{010}}\\
	\end{array}\right],
\end{equation}
and the system (\ref{Ns_Map}) can be described as
\begin{equation*}
\left(\begin{array}{c}
X \\
Y \\
\end{array}\right) \rightarrow\left(\begin{array}{cc}
\frac{a_{100}+b_{010}}{2} & \frac{\sqrt{4 n(q_2)-m^2(q_2)}}{2}  \\
-\frac{\sqrt{4 n(q_2)-m^2(q_2)}}{2} & \frac{a_{100}+b_{010}}{2}  \\
\end{array}\right)\left(\begin{array}{c}
X \\
Y \\
\end{array}\right)+\left(\begin{array}{c}
\hat{f}(X, Y) \\
\hat{g}(X, Y)\\
\end{array}\right),
\end{equation*}
where
\begin{equation*}
\begin{aligned}
\hat{f}(X, Y)\!=&\!A_{200}X^2 \!+ A_{110}X Y \!+A_{020} Y^2 \!+A_{300}X^3 \!+A_{210}X^2 Y \!+A_{120}X Y^2  \!+A_{030}Y^3\!+O\left(|x, y|^4\right),  \\
\hat{g}(X, Y)\!=&\!B_{200}X^2 \!+ B_{110}X Y \!+B_{020} Y^2 \!+B_{300}X^3 \!+B_{210}X^2 Y \!+B_{120}X Y^2 \!+ B_{030}Y^3\!+O\left(|x, y|^4\right),
\end{aligned}
\end{equation*}
and $A_{200}=\frac{4 a_{010}^2 a_{200} + 2 a_{110} (b_{010} - a_{100}) a_{010} + a_{020} (a_{100}-b_{010})^2}{4 a_{010}^2}$, \\
$A_{110}=\frac{(a_{020} (b_{010} - a_{100}) + a_{010} a_{110}) \sqrt{4 n(q_2)-m^2(q_2)}}{2 {a_{010}}^2}$, \\
$A_{020}=\frac{a_{020} (4 n(q_2)-m^2(q_2))}{4 {a_{010}}^2}$, \\
$A_{300}=\frac{8 a_{300} {a_{010}}^3 - 4 (a_{100}-b_{010}) a_{210} {a_{010}}^2 + 2 {(a_{100}-b_{010})}^2 a_{120} a_{010} - {(a_{100}-b_{010})}^3 a_{030}}{8 {a_{010}}^3}$, \\
$A_{210}=\frac{(3 a_{030} {a_{100}-b_{010}}^2 -a_{120} a_{010} (a_{100}-b_{010}) + a_{210} {a_{010}}^2) \sqrt{4 n(q_2)-m^2(q_2)}}{8 {a_{010}}^3}$, \\
$A_{120}=\frac{3 a_{030} (b_{010}-a_{100} + 2 a_{120} a_{010})(4 n(q_2)-m^2(q_2))}{8 {a_{010}}^3}$, \\
$A_{030}=\frac{a_{030} \left(4 n(q_2)-m^2(q_2)\right)^{3/2}}{8 {a_{010}}^3}$, \\
$B_{200}=\frac{8 b_{200} {a_{010}}^3 + 4 (a_{200} - b_{110}) (a_{100}-b_{010}) {a_{010}}^2 - 2 {(a_{100}-b_{010})}^2 (a_{110} - b_{020}) a_{010} + {(a_{100}-b_{010})}^3 a_{020}}{4 a_{010}^2 \sqrt{4 n(q_2)-m^2(q_2)}}$, \\
$B_{110}=\frac{2 b_{110} {a_{010}}^2 + (a_{110} - 2 b_{020}) (a_{100}-b_{010}) a_{010} - a_{020} {(a_{100}-b_{010})}^2}{2 {a_{010}}^2}$,\\
$B_{020}=-\frac{(a_{020} (b_{010} - a_{100})- 2 b_{020} a_{010}) \sqrt{4 n(q_2)-m^2(q_2)}}{4 {a_{010}}^2}$, \\
$B_{300}=\frac{16 b_{300} {a_{010}}^4 - 8 (a_{100}-b_{010}) (b_{210} - a_{300}) {a_{010}}^3 - 4 {(a_{100}-b_{010})}^2 (a_{210} - b_{120}) {a_{010}}^2 + 2 {(a_{100}-b_{010})}^3 a_{120} a_{010} - {(a_{100}-b_{010})}^4  a_{030}}{8 {a_{010}}^3 \sqrt{4 n(q_2)-m^2(q_2)}}$,\\
$B_{210}=\frac{8 b_{210} {a_{010}}^3 + 4 (a_{210} - 2 b_{120}) (a_{100}-b_{010}) {a_{010}}^2 - 4 {(a_{100}-b_{010})}^2 a_{120} a_{010} + 3 {(a_{100}-b_{010})}^3 a_{030}}{8 {a_{010}}^3}$, \\
$B_{120}=\frac{(4 b_{120} {a_{010}}^2 + 2 a_{120} a_{010} (a_{100}-b_{010}) - 3 a_{030} (a_{100}-b_{010})^2) \sqrt{4 n(q_2)-m^2(q_2)}}{8 {a_{010}}^3}$, \\
$B_{030}=\frac{a_{030} (b_{010}-a_{100}) (4 n(q_2)-m^2(q_2))}{8 {a_{010}}^3}$.

Thus, the additional nondegeneracy condition for NS bifurcation is given by
$$
\hat{a}=-\operatorname{Re}\left(\frac{(1-2 \lambda_1) {\lambda_2}^2}{1-\lambda_1} \varrho_{11} \varrho_{20}\right)-\frac{1}{2}\left|\varrho_{11}\right|^2-\left|\varrho_{02}\right|^2+\operatorname{Re}\left(\lambda_2 \varrho_{21}\right),
$$
where
$$
\begin{aligned}
& \varrho_{20}=\left.\frac{1}{8}\left[\hat{f}_{XX}-\hat{f}_{Y Y}+2 \hat{g}_{X Y}+i\left(\hat{g}_{XX}-\hat{g}_{Y Y}-2 \hat{f}_{X Y}\right)\right]\right|_{(0,0)}, \\
& \varrho_{11}=\left.\frac{1}{4}\left[\hat{f}_{XX}+\hat{f}_{Y Y}+i\left(\hat{g}_{XX}+\hat{g}_{Y Y}\right)\right]\right|_{(0,0)}, \\
& \varrho_{02}=\left.\frac{1}{8}\left[\hat{f}_{XX}-\hat{f}_{Y Y}-2 \hat{g}_{X Y}+i\left(\hat{g}_{XX}-\hat{g}_{Y Y}+2 \hat{f}_{X Y}\right)\right]\right|_{(0,0)}, \\
& \varrho_{21}=\left.\frac{1}{16}\left[\hat{f}_{XXX}-\hat{f}_{X Y Y}+\hat{g}_{XX Y}+\hat{g}_{Y Y Y}+i\left(\hat{g}_{XXX}-\hat{g}_{X Y Y}-\hat{f}_{XX Y}-\hat{f}_{Y Y Y}\right)\right]\right|_{(0,0)} .
\end{aligned}
$$

According to the references \cite{liu2007complex, liu2018bifurcation}, when $\hat{a}<0$, a stable invariant cycle bifurcates from the fixed point $E_{3}=(x_{3}, y_{3})$ via the NS bifurcation at $q_2={q_2}^*$; when $\hat{a} > 0$, an unstable invariant cycle bifurcates from $E_3$ via an NS bifurcation at $q_2 = {q_2}^*$. This completes the proof.
\end{proof}
