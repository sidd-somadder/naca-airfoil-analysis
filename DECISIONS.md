# Engineering Decisions Log

A record of technical decisions made during development, with the reasoning behind each. Entries are ordered primarily by significance.

---

## 1. Boundary Condition

**Issue:** A constant-strength vortex panel method requires a boundary condition to close the system for the unknown vortex strengths. Two options apply: enforcing zero normal velocity at the surface, which states impermeability directly, or enforcing zero tangential velocity at an interior point, which Katz & Plotkin (Section 11.2.3) recommends for this element type because the self-induced normal velocity at a panel's own collocation point is identically zero (Eq. 10.43).

**Decision:** The Neumann zero-normal-velocity condition, implemented with the K and L geometric-integral influence matrices.

**Reasoning:** The internal tangential condition was implemented first, following the recommendation. It returned circulation roughly 5.9× below thin airfoil theory. Four tests helped isolate the cause:

1. Panel-count sweep was run from N = 50 to N = 800. The ratio to TAT stayed at ~5.9. It was decided not a discretization error, which would otherwise change under refinement.

2. Thickness check using NACA 0001 test airfoil, and still ~4.7× off at 1% thickness. At 1% thickness the VPM should converge to thin airfoil theory by construction, so a persistent 4.7× discrepancy ruled out geometry fidelity.

3. The tangential-condition implementation built its influence matrix from a per-panel induced-velocity routine. That routine was checked with a single vortex panel against the analytical point-vortex limit $|V| = \frac{\Gamma}{2\pi r}$ for increasing $r$. The result of this test was that the ratio of theoretical $|V|$ to the induced velocity magnitude tended to 1. Hence, the induced velocity algorithm did not explain the discrepancy.

4. The solver was rebuilt temporarily using the classic sharp-TE scheme with no added closure panel. Results were unchanged, ruling out the trailing-edge treatment.

With discretization, geometry, the velocity routine, and trailing-edge treatment eliminated, the formulation itself remained. It was decided to attempt the Neumann boundary condition since there was no clear way forward to address the failure of the internal-tangential condition.

The Neumann condition applies a zero normal velocity condition on each panel's collocation point. Physically, this is intuitive since there is no fluid flow through a solid boundary.

Implementing the Neumann condition produced circulation in agreement with XFOIL inviscid results, with $C_L$ matching to within 0.04 across the sweep. Note that Katz flags the internal-tangential scheme as requiring modifications beyond what is presented, which is consistent with what was observed here. The conditioning of the resulting system is examined in decisions log #2.

---

## 2. Pressure Oscillation Diagnosis

**Issue:** The constant-strength pressure distribution has a sawtooth panel-to-panel oscillation that worsens with mesh refinement rather than improving. Since the lift coefficient computed via Kutta Joukowski agreed with XFOIL results, it was unclear if this was an implementation error or a property specific to the constant-strength formulation.  

**Decision:** A singular value decomposition (SVD) was run on the K geometric influence matrix used in solving for vortex strength distribution.

**Reasoning:** A smoothing algorithm was considered first and rejected, since it would remove information about the actual behaviour of the solved vortex strength.

Running SVD on the K matrix before implementing the Kutta Condition yielded a singular matrix, meaning there was no unique solution. Once the Kutta condition was implemented, the matrix was no longer singular, making it solvable.

With the Kutta condition applied, cond(K) grows in proportion to $N^2$. As N increases, the system increasingly fails to constrain the family of solutions where the vortex strength alternates sign from panel to panel. These alternating values cancel when summed along the surface, so lift remains accurate, but no cancellation occurs at each individual panel, which is what produces the sawtooth pattern.

This confirms the oscillation is a property of the constant-strength formulation rather than an implementation error, and that lift results remain usable while local pressure values do not.

---

## 3. Open Trailing Edge

**Issue:** The standard NACA thickness equation produced a nonzero thickness at $x/c = 1$, so the generated airfoils have blunt trailing edges (TE) by construction.  

**Decision:** An explicit closure panel was added connecting the last upper and lower points, and is included in the K influence matrix like all other panels. Original NACA thickness coefficients kept rather than modified.

**Reasoning:** Force closing the TE by setting the upper and lower endpoints y-values to 0 distorts geometry specifically where the Kutta condition is applied, and so where the unique solution for circulation is distorted too.

The NACA coefficients are the published definition of the section, and rather than modify them to suit the method, the method was built around the blunt TE as defined. Since the closure panel handles it directly, there was no reason to alter the published geometry.

Leaving the gap unpaneled leaves no boundary at the TE to enforce a condition on. This means the surface never closes, hence not a solid boundary, and the Kutta condition has no complete geometry to act on.

Keeping the closure panel is similar to how Drela's XFOIL Inviscid solver handles open TE, though the closure panel here carries vortex strength rather than source strength.

Since the closure panel is a mathematical artifact used to create a closed loop of points, the Kutta condition applies to the upper and lower adjacent panels and the equation for this closure panel is replaced by it.

---

## 4. Two-way $C_L$ Calculation 

**Issue:** The pressure distribution had no independent check for accuracy even if the oscillatory pattern was a consequence of ill-conditioning as discussed in #2.

**Decision:** $C_L$ was computed in two ways from the same solved $\gamma$ vector, (1) integrating pressure across the upper/lower surface panels, and (2) summing the local vortex strength across all panels, and then applying the Kutta-Joukowski theorem.

**Reasoning:** The two algorithms utilize the same circulation solution in different ways to derive the quantity. Hence, their agreement would verify the pressure distribution provided that $C_L$ computed from Kutta-Joukowski agrees with XFOIL Inviscid results. 

In practice, the two quantities agreed to within 0.01 across the angle sweep, greatest at $\alpha = 20 \degree$, and tended to zero as $\alpha \rightarrow 0 \degree$.

Note that due to panel masking as discussed in decision entry #6, the two algorithms use slightly different panel sets. Hence, a small residual is expected rather than being indicative of implementation error.

---

## 5. Half Cosine Spacing

**Issue:** The NACA airfoil thickness function has a $\sqrt{x}$ term, and so an infinite slope at the leading edge.

**Decision:** Half cosine spacing was used over uniform spacing, clustering panels at LE and TE.

**Reasoning:** 
Since the airfoil shape changes the greatest at the leading edge and trailing edge, discretizing such that the geometry is more accurate and precise at those locations is preferred. Hence, half-cosine spacing clusters panels at these two locations while uniform does not.

A measured consequence is that the first tangential angle for the NACA 0012 (199 panels) case is ~85° with half-cosine spacing while ~60° with uniform spacing. This is a result of the former placing the next panel node at about 0.025% of chord length while uniform places it at about 1.0%.

Notably, uniform spacing yields a slightly lower condition number (6.1e+3 vs 7.0e+3), so this is a geometric discretization issue rather than a conditioning problem as seen in #2.

---

## 6. TE-Adjacent Panel Masking

**Issue:** Panels indexed 0, M-2, and M-1 produced extremely high $C_P$ values far outside physical range. Since $C_n$ and $C_a$ are sums of $C_P$ weighted by panel length, a few values far outside physical range corrupt the whole integral.

**Decision:** Exclude those three from the pressure integration algorithm, masked with NaN in the length-M array rather than deletion.

**Reasoning:** Masking maintains the length-M array, so it stays index-aligned with the midpoints, panel length, and angle arrays. Otherwise, deletion would require reindexing everywhere to prevent indexing bugs.

These panels were masked since they are the shortest panels produced by cosine spacing and sit adjacent to the open trailing edge.

A consequence of this decision is that the two different $C_L$ algorithms use different panel sets, which is addressed in decisions log #4.

---

## 7. Vortex Panel over Source/Doublet

**Issue:** Katz & Plotkin present source, doublet, and vortex panel formulations. One elementary flow element had to be selected.

**Decision:** Vortex panels were chosen.

**Reasoning:** This project was initially motivated by implementing the vortex panel method discussed in classwork, and XFOIL's inviscid solver itself is a vortex panel method. The comparison is a direct one rather than across panel method types.

Moreover, source/sink flows are defined as radial flow with no tangential element. By definition they can't carry circulation, hence a source panel method, while able to define the flow potential, can't produce a lifting flow. 

Katz discusses how a spatially varying doublet panel method can produce a lifting flow. However, even though this doublet formulation is mathematically equivalent to a vortex formulation, it is more straightforward to use the vortex panel method to solve for circulation.

---

## 8. Reverse-Selig Panel Order

**Issue:** The direction in which coordinates are traversed determines the local unit tangent direction of the panels, which sets the sign convention within the influence matrix construction.

**Decision:** The in-repo NACA generator writes coordinates in reverse-Selig order. 

TE $\rightarrow$ Lower Surface $\rightarrow$ LE $\rightarrow$ Upper Surface $\rightarrow$ TE. 

**Reasoning:** Rather than adapt the discretization to a different convention, the ordering assumed by Katz was used. Otherwise, plausible-looking but mirrored results would be produced.

An early version using standard Selig ordering produced a $C_P$ distribution with the suction peak on the geometric lower surface at positive angle of attack, which is physically backwards. Coloring the upper and lower surfaces separately on the $C_P$ plot made the inversion immediately visible.

---

## 9. Constant Strength over Discrete Vortex Method

**Issue:** Per Katz & Plotkin, there are multiple vortex panel method formulations: discrete, constant-strength, and linear strength options. 

**Decision:** Constant-strength formulation was chosen for the panel method.

**Reasoning:** 
Discrete Vortex Method places vortex points along the mean camber line, which uses the same zero-thickness assumption of thin airfoil theory. Hence, it would be a numerical approximation of the analytic solution rather than a physical model at a different fidelity level. The aim of this project requires the latter.

Constant-strength places panels on the actual upper and lower surfaces of the airfoil, which takes the thickness into consideration, and hence allows us to distinguish airfoils by actual geometry rather than analytic definition.

Linear-strength formulation was not rejected, but deferred as the natural next step for this project pipeline: the conditioning behaviour found in #2 motivates this choice.

---

## 10. Scope Reduction (July)

**Issue:** The project was originally scoped to handle .dat importing, thin airfoil spline approximation for non-NACA geometry, and 5-digit generation. A GUI was planned as well. Timeline constraints left less time than planned.

**Decision:** All four planned features were cut to tighten the scope of the repository. NACA 4-digit geometry supported only with in-house generator.

**Reasoning:** 
The project's main aim is the multi-method comparison, not input generality. Time spent on deepening the comparison is more valuable than time spent on input robustness.

5-digit generation and GUI remain clean future additions. The in-house generator already handles thickness distribution and panel ordering, and the 5-digit NACA codes utilize an analytic camber line function similar to the 4-digit NACA codes. The GUI can be implemented as a wrap around the existing algorithms/methods.

