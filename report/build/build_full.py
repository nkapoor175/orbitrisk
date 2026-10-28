"""Builds the full OrbitRisk IEEE report inside the course's IEEE Word template.

Run from the repo root (python-docx and pandas needed):
    python report/build/build_full.py report/OrbitRisk_IEEE_Report.docx
All numbers in the text and tables are read from data/processed/ or quoted from the verified results.
"""
import math
import re
import sys

from report_lib import *  # noqa: F401,F403  (template document `d`, data frames, helpers)

OUT = Path(sys.argv[1])

# ---------------- title block -----------------
for r in list(P[0]._p.findall(qn("w:r"))):
    P[0]._p.remove(r)
P[1].text = ("OrbitRisk: Learning Two-Line-Element Position Uncertainty from Catalogue History "
             "for Space-Debris Collision Probability Estimation")
a = P[2]
for r in list(a.runs):
    r._r.getparent().remove(r._r)
lines = ["Navika Kapoor (25BCE0835), Parthvi Mishra (25BCE0827), Vishi Joshi (25BBS0082)",
         "SCOPE, Vellore Institute of Technology, Vellore, India",
         "navika.kapoor2025@vitstudent.ac.in, parthvi.mishra2025@vitstudent.ac.in, "
         "vishi.joshi2025@vitstudent.ac.in"]
for i, line in enumerate(lines):
    for part in re.split(r"(\[[^\]]*\])", line):
        if part:
            run = a.add_run(part)
            if part.startswith("["):
                run.font.highlight_color = WD_COLOR_INDEX.YELLOW
    if i < len(lines) - 1:
        run.add_break()

ab = P[3]
for r in list(ab.runs):
    r._r.getparent().remove(r._r)
r = ab.add_run("Abstract—")
r.bold = True
r.italic = True
ab.add_run(
    "Whether a satellite should move out of the way of a piece of debris is decided from the probability of "
    "collision (Pc), and Pc needs to know how wrong the predicted positions might be. The free public catalogue "
    "of two-line element sets (TLEs) gives positions but no uncertainty at all. We estimate that uncertainty "
    "from the catalogue’s own history. For 300 objects in the 700–900 km shell we propagated each TLE with SGP4 "
    "to the epoch of a later TLE of the same object, giving 145,467 residuals split along the three orbit "
    "axes. The error scale grows roughly linearly with TLE age, and the along-track scale is "
    "about three times the other two. A Student-t distribution fits the standardised residuals much better "
    "than a Gaussian (AIC lower by 19,245 to 364,512; 7.1, 1.0 and 4.5 degrees of freedom), and its 95% "
    "intervals cover 95.0–96.9% of held-out residuals. A simulation shows that the Wilcoxon test we first used "
    "to compare the models rejects a true Gaussian 10–30% of the time instead of 5%. Screening one week gave "
    "8,172 satellite–debris conjunctions. With the learned Gaussian scales the exact Pc has a median of "
    "2×10⁻¹⁸ and a maximum of 4.3×10⁻⁴, and a Poisson model gives a 0.22% chance of at least one collision. "
    "Importance sampling, which matches the exact integral to within 0.07%, lets us compute Pc under the "
    "Student-t model too. The total risk hardly moves, but the list of risky encounters changes: 38 exceed "
    "10⁻⁵ instead of 18, and a Gaussian alert rule misses 23 of them even though its ROC area is 0.995.")

body_el = d.element.body
started = False
for el in list(body_el):
    if el is ab._p:
        started = True
        continue
    if started and el.tag != qn("w:sectPr"):
        body_el.remove(el)

kw = d.add_paragraph(style="IndexTerms")
r = kw.add_run("Keywords—")
r.bold = True
r.italic = True
kw.add_run("space debris, collision probability, two-line elements, SGP4, Student-t distribution, maximum likelihood, "
           "Monte Carlo simulation, importance sampling, Bayesian credible interval, hypothesis testing")

# ---------------- I. Introduction -----------------
h1("Introduction")
body("Low Earth orbit is crowded. Next to the satellites that still work there is far more junk: retired "
     "spacecraft, upper stages left behind by launches, and pieces created when objects broke apart. Kessler "
     "and Cour-Palais [1] warned that once this junk is dense enough, crashes make new pieces faster than the "
     "thin upper atmosphere can pull them down, and the amount of debris then keeps rising by itself. "
     "Preventing a single crash therefore helps twice: it saves the satellite, and it keeps the altitude "
     "usable for everyone else.")
body("When two objects are predicted to pass near each other, an event called a conjunction, the number an "
     "operator watches is Pc. It answers one question: at the moment the two are nearest (TCA), how likely "
     "is it that the real gap between their centres is smaller than their two sizes added together? Neither "
     "position is known exactly, so answering it means integrating the error distribution over that tiny "
     "region. A value of about 10⁻⁴ or more [2] usually starts a discussion about moving the satellite.")
body("The catch is the error distribution. Satellite operators who track their own spacecraft get a "
     "covariance with every orbit estimate. Everyone else depends on the public catalogue that the US "
     "military tracking network releases through Space-Track [5]. Its two-line element sets are averaged "
     "orbital elements meant to be fed into the SGP4 propagator [3], [4], and not one of them says how "
     "accurate it is. So a student, a researcher or a small operator "
     "computing Pc from public data has to guess a covariance, usually a Gaussian with hand-picked standard "
     "deviations, and the answer is only as good as the guess.")
body("Our idea is to stop guessing and learn the error from the catalogue itself. An object’s TLE is "
     "updated every few hours to days. If we take an older TLE, propagate it to the time of a newer one and "
     "compare, the difference is exactly the error someone using the older TLE would have made. Doing this "
     "145,467 times gives enough data to model the error properly. The work has three parts. First, we "
     "measure how fast the error grows as a TLE gets older, separately along the three orbit axes defined in "
     "Section IV. Second, we ask whether the error is normally distributed, with likelihood fits, AIC, "
     "goodness-of-fit distances, a check on held-out data and a simulation of how often our test makes each "
     "kind of mistake; it is not, because its tails are far too heavy. Third, we search the catalogue for "
     "close satellite–debris passes and compute Pc for every one, exactly for the normal model and by Monte "
     "Carlo and importance sampling, with Bayesian intervals, for the heavy-tailed one, and finally combine "
     "the probabilities in a Poisson model.")
body("Section II states the problem and Section III reviews earlier work. Section IV covers the data and "
     "the methods, Section V the results, Section VI what they mean and where they fall short, and Section VII "
     "concludes.")

# ---------------- II. Problem statement -----------------
h1("Problem Statement")
body("Take a satellite and a debris object with predicted positions r₁ and r₂ at TCA and relative velocity "
     "v. Looking along v, the encounter plane is perpendicular to it, and in that plane the predicted miss "
     "vector μ has two components. The true miss vector is μ plus the position errors of both objects. A hit "
     "happens when this true vector ends up within a distance R of the origin, where R adds up the sizes we "
     "give the two objects; (13) states the probability of that.")
body("To evaluate (13) we need the error distribution of each object at the age its TLE has reached at TCA. "
     "Our problem is therefore: using nothing but public TLEs, estimate that distribution and use it to get "
     "collision probabilities we can trust. In statistical terms this splits into three questions. (a) How "
     "does the TLE position error grow with TLE age in each direction? (b) Is it Gaussian, as almost "
     "everyone assumes, or does it have heavier tails? (c) How much do Pc values, and the alerts based on "
     "them, change when we replace the Gaussian with the learned model?")
body("We restrict ourselves to the 700–900 km shell, which is crowded and full of debris. The error model "
     "uses 120 days of TLE history for 300 objects, the conjunction search covers seven days, and everything "
     "is written in Python and R.")

# ---------------- III. Literature review -----------------
h1("Literature Review")
h2("TLE accuracy and error estimation")
body("The SGP4 propagator was described in Spacetrack Report No. 3 [3]. In 2006 Vallado et al. [4] "
     "released a corrected and consistent version of the code, and that version is what our scripts call. "
     "Since TLEs have no covariance, their accuracy has been "
     "estimated indirectly. Flohrer, Krag and Klinkrad [6] did it for whole snapshots of the US catalogue by "
     "fitting orbits to pseudo-observations made from the TLEs, and grouped the objects into accuracy classes "
     "that ESA uses when it assesses collision risk. Osweiler [7] estimated covariances for single satellites "
     "from two weeks of TLE-to-TLE differences. Geul, Mooij and "
     "Noomen [8] made the differencing idea robust with weighted least squares, estimated both the starting "
     "uncertainty and how fast it grows, and checked the result against GPS orbits of GOCE.")
body("So TLE differencing is an accepted way to measure TLE error, and it is known that the error grows "
     "with propagation time. What these studies report, though, is a standard deviation, which quietly "
     "treats the error as Gaussian. None of them takes the learned distribution on to a collision "
     "probability.")
h2("Collision probability methods")
body("The Pc used in operations goes back to Foster and Estes [9]. A conjunction lasts only a moment, so "
     "the relative motion can be treated as a straight line, and the volume integral shrinks to an area "
     "integral over a disk in the plane at right angles to the relative velocity. Akella and Alfriend [10] "
     "derived the collision "
     "probability for Gaussian uncertainty in this setting, Alfano [11] turned the double integral into a "
     "single one evaluated with Simpson’s rule, and Chan [12] wrote a full treatment with series solutions. "
     "Every one of these methods assumes Gaussian errors. The covariance is then the only input that counts, "
     "and heavy tails are simply not part of the model.")
h2("Statistical tools")
body("Lange, Little and Taylor [13] showed that fitting a t distribution, and letting the data choose its "
     "degrees of freedom, makes estimates resistant to outliers; pushing that parameter to infinity gives "
     "back the normal curve. To choose between models we use Akaike’s criterion [14], which trades fit "
     "against the number of parameters, and the Wilcoxon signed-rank test [15] when two models are compared "
     "residual by residual. Because a collision almost never happens, a naive simulation wastes nearly all "
     "of its draws; importance sampling [16] puts the draws where hits occur and corrects with weights. When "
     "an estimate is a share of hits, Brown, Cai and DasGupta [17] found the Jeffreys interval reliable even "
     "with only a handful of hits.")
table("Related work and the gap each leaves",
      ["Ref.", "Contribution", "Gap for this work"],
      [["[3], [4]", "SGP4 model and reference code", "Positions only, no uncertainty"],
       ["[6]", "Catalogue-wide TLE error estimates", "Gaussian summary, no tails"],
       ["[7]", "TLE covariance from differences", "Few objects, no age law"],
       ["[8]", "Robust differencing, error growth", "No collision probability"],
       ["[9]–[12]", "2D Pc formulation and integration", "Gaussian errors assumed"],
       ["[13]", "Student-t robust modelling", "Not applied to TLE error"],
       ["[16], [17]", "Rare-event sampling, intervals", "Combined here for Pc"]],
      [0.55, 1.4, 1.35])
body("Put simply, the pieces exist separately: a propagator, the differencing idea, the Pc integral and the "
     "statistical tools. What is missing is an error model that depends on TLE age, allows heavy tails, is "
     "learned from the public catalogue, and is carried all the way to Pc with its uncertainty. That is what "
     "we built.")

# ---------------- IV. Methodology -----------------
h1("Methodology")
h2("Data and preprocessing")
body("We downloaded TLE histories from Space-Track [5] for 300 objects (120 satellites and 180 debris "
     "objects) with mean altitudes between 700 and 900 km, covering 8 June to 6 October 2026. Space-Track "
     "limits how fast it can be queried, so the download asks for 25 objects at a time with a pause after "
     "each request. Rows with a mean altitude outside 300–1200 km and objects with fewer than five TLEs were "
     "removed. Manoeuvres needed care: when a satellite fires its thrusters, the next TLE is far from the "
     "propagated old one, but that is not a prediction error. A manoeuvre shows up as a jump in semi-major "
     "axis, so any jump above 0.5 km between consecutive TLEs was flagged and no residual pair was allowed to "
     "cross one. Finally, 129 pairs with an error component above 100 km were dropped as corrupted TLEs. "
     "Table II summarises the data, and Fig. 1 shows the whole pipeline.")
figure(REPO / "report" / "figures" / "report_pipeline.png", "Processing pipeline of OrbitRisk.", width=3.3)
table("Dataset summary",
      ["Item", "Value"],
      [["Source", "Space-Track.org TLE history [5]"],
       ["Orbital shell (mean altitude)", "700–900 km"],
       ["Objects in error model", "300 (120 satellites, 180 debris)"],
       ["TLEs after cleaning", f"{len(tle):,}"],
       ["TLE epochs", "8 Jun – 6 Oct 2026 (120 days)"],
       ["Manoeuvre rule", "semi-major-axis jump > 0.5 km"],
       ["Residuals (TLE pairs)", f"{len(res):,}, ages 0.1–7 days"],
       ["Residual cap", "pairs with any component > 100 km dropped"],
       ["Objects in conjunctions", f"{conj.norad_id_1.nunique()} satellites, {conj.norad_id_2.nunique():,} debris"],
       ["Screening window", "6–13 Oct 2026 (7 days)"],
       ["Conjunctions kept", f"{len(conj):,} (miss < 5 km, rel. speed > 0.5 km/s)"]],
      [1.5, 1.8])

h2("Residuals in the orbital frame")
body("For each TLE we picked up to two newer TLEs of the same object at random, between 0.1 and 7 days "
     "later. Two per TLE keeps the data set a manageable size while still spreading the ages over the whole "
     "range. We run SGP4 on the older TLE up to the newer TLE’s epoch t₂ and compare with the newer TLE’s own "
     "position, the closest thing to the truth we have; the difference is Δr in (1), and t is the TLE age. "
     "An error up or down, an error forward or back along the path, and an error sideways out of the orbit "
     "plane have different causes, so Δr is split along the radial (R), along-track (T) and cross-track (N) "
     "axes of (2), which gives the three numbers in (3).")
eq("Δr = r_new(t₂) − SGP4(TLE₁, t₂),   t = t₂ − t₁", 1)
eq("unit vectors: R̂ ∥ r,  N̂ ∥ r × v,  T̂ = N̂ × R̂", 2)
eq("(δ_R, δ_T, δ_N) = (R̂·Δr, T̂·Δr, N̂·Δr)", 3)

h2("Age-dependent error scale")
body("Older TLEs should be less accurate. We model the spread with the power law (4), which can describe "
     "linear growth (b = 1) as well as slower or faster growth. The ages are cut into 20 bins with the same "
     "number of residuals, and in each bin the scale is measured with the median absolute deviation (5). For "
     "Gaussian data this equals the standard deviation, but a handful of huge residuals cannot drag it up. A "
     "straight-line fit of log scale against log age (6) then gives a_j, the scale at one day, and the "
     "growth exponent b_j.")
eq("σ_j(t) = a_j · t^{b_j},   j ∈ {R, T, N}", 4)
eq("σ̂_j = 1.4826 × MAD(δ_j)   (MAD: median absolute deviation)", 5)
eq("ln σ̂_j(t) ≈ ln a_j + b_j · ln t", 6)

h2("Distribution models and maximum likelihood")
body("Dividing each residual by the scale for its age (7) removes the age effect, so all standardised "
     "residuals z of one direction can be pooled to study their shape. Two shapes were tried. The normal "
     "model (8) needs only a centre μ and a spread σ. The t model (9) adds a third number, ν: when ν is small "
     "the tails are fat, and when it is large the curve can hardly be told apart from the normal. In effect "
     "the data decide how much tail they need. Both models were fitted by maximising the log-likelihood (10): "
     "the Gaussian has a closed-form solution, and the Student-t was fitted numerically with BFGS in R. Our "
     "first fit forced ν > 2, so that the variance exists, and the along-track estimate got stuck exactly at "
     "that bound. We removed the restriction, and ν then settled at about 1. The models are compared with AIC "
     "(10); k counts the fitted parameters, two for the normal and three for t, and the smaller AIC wins.")
eq("z = δ_j / σ_j(t)", 7)
eq("f_N(z) = (1/(σ√(2π))) exp(−(z−μ)²/(2σ²))", 8)
eq("f_t(z) ∝ (1/s) · (1 + (z − m)² / (ν s²))^{−(ν+1)/2}", 9)
eq("θ̂ = argmax_θ Σ_i log f(z_i; θ),   AIC = 2k − 2 log L(θ̂)", 10)

h2("Hypothesis tests, error rates and calibration")
body("We also used three checks. The Kolmogorov–Smirnov distance (11) measures how far the cumulative "
     "curve of the data, F_n, strays from the fitted curve F at its worst point. F was fitted to these same "
     "residuals, so we use D_n only to describe the fit, not to get a p-value. Then, for each residual we "
     "took the difference in "
     "log density between the two models, d_i = log f_t(z_i) − log f_N(z_i), and applied the two-sided "
     "Wilcoxon signed-rank test [15] of H₀: d is centred at zero.")
eq("D_n = largest |F_n(z) − F(z)| over all z", 11)
body("A test is only as good as its error rates, so we measured them. Rejecting H₀ when it is in fact "
     "true is a Type I error, with rate α; keeping H₀ when H₁ is in fact true is a Type II error, with rate "
     "β, and 1 − β is the power. To estimate both we generated fake data sets of 50, 100, 200, 500 and 1000 "
     "points, once from the fitted normal, where H₀ holds, and once from the fitted t, where H₁ holds, "
     "refitted both models every time, ran the test at the 5% level, and counted rejections over 1000 "
     "repetitions. The third check was calibration on data the models had not seen: "
     "the residuals were split at random 80/20, the models were fitted on the 80%, and we counted how many of "
     "the other 20% fell inside each model’s central 95% interval.")
body("To see how much the Gaussian assumption changes actual decisions, we treated the Student-t Pc as the "
     "reference and called a conjunction risky if its Student-t Pc is at least 10⁻⁵. We first tried 10⁻⁴, "
     "but only one conjunction passed, which makes any comparison meaningless, so we moved to 10⁻⁵. The "
     "Gaussian Pc is the score. Sweeping a threshold over it traces the receiver operating characteristic "
     "(ROC) curve, and we report the area under it (AUC) together with the false-alarm and missed-alert "
     "rates at the same 10⁻⁵ threshold.")

h2("Conjunction screening and the encounter plane")
body("For screening we took a snapshot of the catalogue: the latest TLE of every satellite and debris "
     "object in the shell that had been updated in the last 14 days. These were propagated for the seven days "
     "after the newest epoch. Comparing every satellite with every debris object at every step would be far "
     "too slow, so at each 30 s step we built k-d trees [18] of the satellite and debris positions and only "
     "pulled out pairs closer than a search radius. The radius is the 5 km limit plus a 5 km margin plus the "
     "distance two objects can close in half a step, so that no encounter can slip between steps. Each "
     "candidate was then refined at 1 s resolution to find the TCA and miss distance. We checked the method "
     "against a slow brute-force search, and it found all 14 test encounters. We kept encounters closer than "
     "5 km with a relative speed above 0.5 km/s; slower pairs are objects moving almost together, for which "
     "the short-encounter model does not apply. In every pair, object 1 is the satellite and object 2 the debris.")
body("In the encounter frame, z points along the relative velocity, x along the part of the satellite’s "
     "orbit normal that lies in the encounter plane, and y = z × x. For each object we take the R, T and N "
     "scales at its own TLE age, clamped to the fitted range of 0.1–7 days, as a diagonal covariance in its "
     "orbital frame. This is rotated into inertial axes by the frame matrix M_k. The two objects’ errors do "
     "not depend on each other, so their covariances can be summed, and multiplying by P, whose two columns "
     "are the plane’s x and y axes, gives the 2×2 covariance C in (12).")
eq("C = Pᵀ(Σ₁ + Σ₂)P,   Σ_k = M_k diag(σ_R², σ_T², σ_N²) M_kᵀ", 12)

h2("Collision probability")
body("Equation (13) is the short-encounter formula of [9]–[12] for normal errors: a 2D normal density, "
     "centred on the predicted miss vector μ with covariance C, integrated over a disk of radius R. We "
     "assumed 2.0 m for every satellite and 0.5 m for every piece of debris, so R = 2.5 m; Section V tests "
     "how much this choice matters. For the normal model "
     "we evaluate (13) by numerical quadrature in polar coordinates over the circle, with Gauss–Legendre "
     "points in radius and evenly spaced points in angle. The code passes eight tests against cases with a "
     "known answer and against Monte Carlo. The Student-t model has no closed form, so there Pc has to be "
     "estimated by sampling.")
eq("Pc = ∬ over the disk ‖x‖ ≤ R of  N₂(x; μ, C) dx", 13)

h2("Monte Carlo, importance sampling and Bayesian intervals")
body("Plain Monte Carlo draws n miss vectors from the error model and counts the k that land inside the "
     "circle (14). It is unbiased, but its relative standard error is about 1/√(n·Pc): to estimate a Pc of "
     "10⁻⁶ to within 10% it would need around 10⁸ samples. Importance sampling [16] draws from a proposal "
     "density g that is concentrated on the hard-body circle and weights each hit by f/g (15). The estimate "
     "stays unbiased and its variance drops by orders of magnitude. For the Student-t model each error is "
     "drawn directly as location plus scale times a t variate with ν degrees of freedom, which still works "
     "when ν ≤ 2 and the variance is infinite.")
eq("P̂c = k/n,   Var(P̂c) = Pc(1 − Pc)/n", 14)
eq("P̂c,IS = mean over draws x_i ~ g of [hit(x_i) · f(x_i) / g(x_i)]", 15)
body("For Monte Carlo we give a Bayesian credible interval. Starting from the Jeffreys prior Beta(½, ½) "
     "and observing k hits in n draws, Bayes’ rule gives the Beta posterior in (16). Cutting 2.5% off each "
     "end of it leaves a 95% credible interval that behaves sensibly even with no hits at all [17]. We report the "
     "posterior mean (k + ½)/(n + 1). One consequence caught us out at first: with no hits it equals "
     "½/(n + 1), not zero, so these values must not be added up across conjunctions.")
eq("Pc | k ~ Beta(k + ½, n − k + ½)", 16)
body("Importance-sampling weights are not 0/1 hits, so the Beta posterior does not apply to them. Each "
     "weight is, however, bounded by a known constant B (the disk area times the highest value the density "
     "can reach, divided by the Jacobian), and Hoeffding’s inequality then gives the 95% confidence interval "
     "(17). It holds for any sample size and any weight distribution, but it is wide. Both Monte Carlo and "
     "importance sampling used n = 10⁵ samples per conjunction and model. With 10⁴ importance samples the "
     "Student-t estimates were still noisy, which is why we went up to 10⁵. For the three riskiest "
     "conjunctions we also ran Monte Carlo with n = 10⁶. Finally, we recomputed Pc for R = 1, 2.5 and 5 m to "
     "see how much the assumed object size matters.")
eq("P̂c ± B √( ln(2/α) / (2n) ),   α = 0.05", 17)

h2("Fleet-level Poisson model")
body("Each conjunction is a rare yes/no event with probability Pc_i, and different conjunctions are nearly "
     "independent. The number of collisions N in the window is then a sum of many rare independent "
     "indicators, which is well approximated by a Poisson variable whose mean is the sum of the "
     "probabilities (18). From it we get the expected number of collisions and the chance of at least one.")
eq("N ~ Poisson(Λ),  Λ = Σ_i Pc_i,  P(N ≥ 1) = 1 − e^{−Λ}", 18)

# ---------------- V. Results -----------------
h1("Experimental Results")
h2("Error growth with TLE age")
body("Fig. 2 plots the median absolute along-track residual in each of the 20 age bins on log–log axes. "
     "The points fall close to a straight line, from about 0.03 km for TLEs a few hours old to about 0.9 km "
     "for TLEs six to seven days old, which is what the power law (4) predicts. Table III lists the fitted "
     "parameters. The exponent b is 1.00 radially and cross-track and 1.10 along-track, so the error grows "
     "almost exactly in proportion to TLE age. At one day the along-track scale is 0.126 km, about three "
     "times the radial (0.040 km) and cross-track (0.048 km) scales. At seven days the fitted scales reach "
     "about 1.07 km along-track, 0.28 km radially and 0.33 km cross-track.")
figure(REPO / "figures" / "final_error_growth_along.png",
       "Median absolute along-track residual per TLE-age bin (log–log axes).")
rows = []
for comp in ["radial", "along", "cross"]:
    g = params[(params.component == comp) & (params.model == "gaussian")].iloc[0]
    s = params[(params.component == comp) & (params.model == "student_t")].iloc[0]
    rows.append([comp.capitalize(), f"{fit.loc[comp, 'a']:.4f}", f"{g.b:.3f}", f"{g.a:.4f}", f"{s.a:.4f}", f"{s.nu:.2f}"])
table("Fitted error model (km, t in days)",
      ["Dir.", "a (MAD)", "b", "a_N", "a_t", "ν"], rows, [0.6, 0.6, 0.5, 0.55, 0.55, 0.5],
      "a_N, a_t: scale at t = 1 day of the fitted Gaussian and Student-t models.")

h2("Gaussian versus Student-t")
body("The standardised residuals are clearly not Gaussian. In the normal Q–Q plot of the along-track "
     "residuals (Fig. 3) the middle of the data follows the reference line, but both ends bend away sharply. "
     "The largest residuals sit hundreds of robust scale units from zero, where a Gaussian allows only a few. "
     "Table IV puts numbers on this. The Student-t model has the lower AIC in all three directions (by "
     "46,519 radially, 364,512 along-track and 19,245 cross-track) and the smaller Kolmogorov–Smirnov "
     "distance, most of all along-track (0.031 against 0.281). The fitted degrees of freedom are 7.1 "
     "radially, 4.5 cross-track and 1.0 along-track (Table III). The last value is practically a Cauchy "
     "distribution, which has no finite variance. The Wilcoxon test rejects ‘no difference’ in every "
     "direction, in favour of the Student-t.")
figure(REPO / "figures" / "qq_along.png", "Normal Q–Q plot of standardised along-track residuals.")
rows = []
for comp in ["radial", "along", "cross"]:
    c = cmp_.loc[comp]
    rows.append([comp.capitalize(), f"{c.AIC_gaussian:,.0f}", f"{c.AIC_student:,.0f}",
                 f"{c.KS_gaussian_statistic:.3f}", f"{c.KS_student_statistic:.3f}"])
table("Model comparison (n = 145,467 per direction)",
      ["Dir.", "AIC Gauss.", "AIC Stud.-t", "KS D Gauss.", "KS D Stud.-t"], rows, [0.6, 0.75, 0.75, 0.6, 0.6],
      "Wilcoxon signed-rank p < 10⁻¹⁰⁰ in all directions, favouring the Student-t model.")
body("The tails also spoil the Gaussian fit itself. Its standard deviation is pulled up by the extreme "
     "residuals, so the Gaussian along-track scale at one day comes out at 0.932 km, more than ten times the "
     "Student-t scale of 0.088 km (Table III). A Gaussian forced onto this data is too wide for ordinary "
     "errors and too thin in the tails at the same time. The radial residuals also have a systematic offset: "
     "the fitted location is 0.94 scaled units.")

h2("Held-out calibration")
body("Table V shows the check on the 29,093 held-out residuals. For the t model, 96.9% (R), 95.8% (T) and "
     "95.0% (N) of these unseen residuals fall inside its 95% interval, never more than two points from the "
     "goal. The Gaussian intervals contain 98.0% and 98.2% radially and along-track, too many because the "
     "inflated standard deviation makes them too wide, and 94.7% cross-track.")
table("Coverage of 95% intervals on the 20% held-out residuals (n = 29,093)",
      ["Direction", "Gaussian", "Student-t"],
      [[k.capitalize(), pct(v[0]), pct(v[1])] for k, v in calib.items()], [1.1, 1.1, 1.1])

h2("How often the model test is wrong")
body("Table VI and Fig. 4 give the simulated rejection rates of the Wilcoxon test at the 5% level. With "
     "truly Gaussian data it should reject 5% of the time. It actually rejects about 10–15% of the time for n "
     "between 50 and 200 and 24–30% at n = 1000, two to six times too often, and the gap grows with the "
     "sample size. With truly Student-t data the power depends on how heavy the tails are. For the "
     "along-track model (ν ≈ 1) the test catches the difference in 94% of samples of 50 and in every larger "
     "sample. For the radial model (ν ≈ 7) the power is only 4–7% up to n = 200 and 45% at n = 1000, so the "
     "Type II error stays large.")
rows = []
for truth, label in [("gaussian", "Type I"), ("student_t", "Power")]:
    for comp in ["radial", "along", "cross"]:
        s = sim[(sim.direction == comp) & (sim.truth == truth)].set_index("n").reject_rate
        rows.append([f"{label}, {comp}"] + [pct(s[n]) for n in (50, 200, 1000)])
table("Simulated rejection rates of the Wilcoxon test (1000 repetitions)",
      ["Truth / direction", "n = 50", "n = 200", "n = 1000"], rows, [1.2, 0.7, 0.7, 0.7])
figure(REPO / "report" / "figures" / "report_type_error.png",
       "Simulated rejection rate of the Wilcoxon test: Type I error (dashed) and power (solid).")

h2("Conjunctions and collision probability")
lam = pc.pc.sum()
body(f"The screening found {len(conj):,} satellite–debris conjunctions within 5 km in the seven days, "
     f"involving {conj.norad_id_1.nunique()} satellites and {conj.norad_id_2.nunique():,} debris objects. The "
     f"median miss distance is {conj.miss_distance_km.median():.2f} km, {(conj.miss_distance_km < 1).sum()} "
     f"passes are closer than 1 km, and the median relative speed is {conj.rel_speed_kms.median():.1f} km/s, "
     "typical of crossing orbits. Table VII summarises the exact Gaussian Pc with the learned scales. Most "
     "encounters are completely harmless: the median Pc is 2.1×10⁻¹⁸, only "
     f"{int((pc.pc >= 1e-6).sum())} conjunctions exceed 10⁻⁶, {int((pc.pc >= 1e-5).sum())} exceed 10⁻⁵, and "
     "just one exceeds the 10⁻⁴ alert level. All probabilities together give "
     f"Λ = {lam:.2e}, so by (18) the chance of at least one collision in the week is about 0.22%. Fig. 5 "
     f"shows the distribution of log₁₀ Pc; {int((pc.pc < 1e-60).sum()):,} conjunctions, "
     f"{int((pc.pc == 0).sum())} of them with a Pc too small for floating point, are below 10⁻⁶⁰ and not shown.")
figure(REPO / "report" / "figures" / "report_pc_hist.png",
       "Distribution of exact Gaussian Pc; dashed lines mark 10⁻⁶, 10⁻⁵ and 10⁻⁴.")
table("Exact Gaussian Pc over all conjunctions",
      ["Quantity", "Value"],
      [["Conjunctions", f"{len(pc):,}"],
       ["Median Pc", f"{pc.pc.median():.1e}"],
       ["Maximum Pc", f"{pc.pc.max():.2e}"],
       ["Pc ≥ 10⁻⁴ (alert)", int((pc.pc >= 1e-4).sum())],
       ["Pc ≥ 10⁻⁵", int((pc.pc >= 1e-5).sum())],
       ["Pc ≥ 10⁻⁶", int((pc.pc >= 1e-6).sum())],
       ["Λ = Σ Pc", f"{lam:.2e}"],
       ["P(N ≥ 1)", f"{1 - math.exp(-lam):.2e}"]],
      [1.6, 1.6])
top = pc.merge(conj[["conj_id", "miss_distance_km", "rel_speed_kms"]], on="conj_id").nlargest(5, "pc")
body("Table VIII lists the five riskiest encounters. The top one, conjunction 183, passes at 169 m and "
     "11.8 km/s with Pc = 4.3×10⁻⁴. Pc does not simply follow the miss distance: conjunction 1198 passes "
     "closer, at 118 m, yet its Pc is ten times smaller. How old the two TLEs are, and how the uncertainty "
     "ellipse is turned compared with the miss direction, matter just as much as the distance.")
table("Five highest-Pc conjunctions (Gaussian, exact)",
      ["ID", "Miss (m)", "Rel. speed (km/s)", "Pc"],
      [[int(r.conj_id), f"{1000 * r.miss_distance_km:.0f}", f"{r.rel_speed_kms:.2f}", f"{r.pc:.2e}"]
       for r in top.itertuples()], [0.6, 0.8, 1.0, 0.9])

h2("Validation of the sampling estimators")
r_all = pd.read_csv(PROC / "pc_results.csv")
wide = r_all.pivot_table(index="conj_id", columns=["method", "error_model"], values="pc")
wide.columns = ["_".join(c) for c in wide.columns]
conv = pd.read_csv(PROC / "mc_convergence.csv")
mc6 = conv[conv.n_samples == 1_000_000].set_index(["conj_id", "error_model"])
isr = r_all[r_all.method == "is"].set_index(["conj_id", "error_model"])
measurable = wide.exact_gaussian > 1e-12
max_rel = float((wide.is_gaussian[measurable] / wide.exact_gaussian[measurable] - 1).abs().max())
body(f"Before trusting the Student-t numbers we checked the samplers against the exact answer. For the "
     f"Gaussian model, importance sampling matches the exact Pc of Table VII to within {100 * max_rel:.2f}% "
     f"for all {int(measurable.sum()):,} conjunctions with Pc above 10⁻¹². Two independent implementations "
     "agreeing this closely is good evidence that both are right. Table IX compares every estimator for the "
     "three riskiest conjunctions. With 10⁵ samples, Monte Carlo for conjunction 183 came out about 2.8 "
     "standard errors above the exact value; the 10⁶-sample run brought it back, and its 95% Jeffreys "
     "intervals now contain the exact value in all three cases. For the Student-t model, importance sampling "
     "and Monte Carlo agree within their intervals. Fig. 6 shows the estimates settling as the sample size "
     "grows. Importance sampling is stable from 10³ samples, while plain Monte Carlo needs 10⁵–10⁶, because at "
     "Pc ≈ 10⁻⁴ only a few hundred of a million samples are hits. The median cost per conjunction was under "
     "3 ms for exact quadrature and Gaussian Monte Carlo, about 11 ms for Gaussian importance sampling and "
     "15–26 ms for the Student-t estimators. The complete pipeline, including the 10⁶-sample runs and the "
     "R simulation study, took about 13.5 hours on a laptop. Rerunning it from the cleaned data reproduced "
     "the exact and Monte Carlo results digit for digit and everything else to within rounding.")


def fmt(x):
    return f"{x:.2e}"


rows = []
for cid in (183, 911, 1198):
    g6, t6 = mc6.loc[(cid, "gaussian")], mc6.loc[(cid, "student_t")]
    rows.append([cid, "Gauss.", fmt(wide.loc[cid, "exact_gaussian"]), fmt(isr.loc[(cid, "gaussian"), "pc"]),
                 f"{fmt(g6.pc)} ({fmt(g6.ci_low)}–{fmt(g6.ci_high)})"])
    rows.append(["", "Stud.-t", "—", fmt(isr.loc[(cid, "student_t"), "pc"]),
                 f"{fmt(t6.pc)} ({fmt(t6.ci_low)}–{fmt(t6.ci_high)})"])
table("Pc of the three highest-risk conjunctions by estimator",
      ["ID", "Model", "Exact", "IS (10⁵)", "MC (10⁶), 95% interval"], rows, [0.4, 0.5, 0.6, 0.6, 1.3])
figure(REPO / "report" / "figures" / "report_convergence.png",
       "Convergence of Monte Carlo (dashed) and importance-sampling (solid) estimates with 95% intervals; "
       "dotted line: exact Gaussian Pc.", width=3.2)
sens = pd.read_csv(PROC / "pc_radius_sensitivity.csv").pivot_table(
    index=["conj_id", "radius_m"], columns="error_model", values="pc")
body("Table X shows how much the assumed object size matters. The density barely changes across a disk a "
     "few metres wide, so Pc grows almost exactly with the disk area: doubling R from 2.5 to 5 m multiplies "
     "Pc by about four in both models. The size assumption scales all probabilities by the same factor and "
     "leaves their order unchanged.")
table("Pc as a function of the hard-body radius R",
      ["ID", "R (m)", "Gaussian", "Student-t"],
      [[int(cid), f"{rad:g}", fmt(v.gaussian), fmt(v.student_t)] for (cid, rad), v in sens.iterrows()],
      [0.5, 0.6, 1.0, 1.0])

h2("Gaussian versus Student-t collision probability")
g_pc, t_pc = wide.exact_gaussian, wide.is_student_t
lam_t = float(t_pc.sum())
body("Switching from the Gaussian to the learned Student-t model changes which encounters are risky much "
     "more than it changes the total risk (Table XI). The riskiest conjunction, 183, drops from 4.2×10⁻⁴ to "
     "1.2×10⁻⁴. Conjunction 1198 doubles, from 4.5×10⁻⁵ to 9.4×10⁻⁵, and conjunction 1605 rises eightfold, "
     "from 7×10⁻⁶ to 5.8×10⁻⁵. Above 10⁻⁵ there are 38 conjunctions instead of 18, and only 15 of them are on "
     "both lists, while above 10⁻⁶ the count falls from 432 to 319. The fleet totals stay close: Λ goes from "
     f"2.18×10⁻³ to {lam_t:.2e}, and the chance of at least one collision from 0.22% to "
     f"{100 * (1 - math.exp(-lam_t)):.2f}%.")
rows = []
for th, lab in [(1e-4, "Pc ≥ 10⁻⁴"), (1e-5, "Pc ≥ 10⁻⁵"), (1e-6, "Pc ≥ 10⁻⁶")]:
    rows.append([lab, int((g_pc >= th).sum()), int((t_pc >= th).sum()), int(((g_pc >= th) & (t_pc >= th)).sum())])
rows.append(["Λ = Σ Pc", f"{g_pc.sum():.2e}", f"{lam_t:.2e}", ""])
rows.append(["P(N ≥ 1)", f"{1 - math.exp(-g_pc.sum()):.2e}", f"{1 - math.exp(-lam_t):.2e}", ""])
table("Risk under the Gaussian and Student-t models",
      ["Quantity", "Gaussian", "Student-t", "Both"], rows, [1.1, 0.75, 0.75, 0.6])
roc = pd.read_csv(PROC / "roc_summary.csv").iloc[0]
body(f"Fig. 7 is the ROC curve of the Gaussian Pc as a predictor of Student-t risk at 10⁻⁵. The AUC of "
     f"{roc.auc:.3f} says the Gaussian model orders encounters from risky to safe almost perfectly. Yet at "
     f"the same 10⁻⁵ threshold it flags only 15 of the {int(roc.student_t_alerts)} Student-t alerts, a "
     f"missed-alert rate of {100 * roc.missed_alert_rate:.1f}%, with a false-alarm rate of "
     f"{100 * roc.false_alarm_rate:.2f}% (3 of {int(roc.student_t_non_alerts):,}). The two models agree on the "
     "order of the encounters but not on how risky they are.")
figure(REPO / "figures" / "final_roc_curve.png",
       "ROC of the Gaussian Pc against Student-t labels (Pc ≥ 10⁻⁵).", width=3.0)

# ---------------- VI. Discussion -----------------
h1("Analysis & Discussion")
body("First, TLE position error grows almost linearly with TLE age in all three directions. That is what a "
     "small constant error in the mean motion would do along-track: a slightly wrong orbital period turns "
     "into a timing error that keeps accumulating. Badly modelled drag adds to it, which is why the "
     "along-track error is both the largest and the fastest growing. For anyone using public data, the age "
     "of the TLE at TCA matters as much as the miss distance. The same geometry gives a very different Pc "
     "with a fresh TLE than with a week-old one.")
body("Second, the errors are heavy-tailed, and this cuts both ways for collision probability. Usually the "
     "error is much smaller than a Gaussian fitted to the same data suggests, because the rare big errors "
     "inflate the Gaussian standard deviation; along-track the Gaussian scale is more than ten times the "
     "Student-t scale. So a Gaussian spreads the probability too thinly over ordinary encounters and, at the "
     "same time, underestimates the rare large error that turns a safe-looking pass into a collision. The "
     "held-out coverage shows the Student-t model getting new residuals right where the Gaussian does not. "
     "The same thing appears in Pc. Conjunction 183, whose risk comes from ordinary-sized errors, becomes "
     "less likely under the Student-t model, while 1198 and 1605, which need a larger error to collide, "
     "become two to eight times more likely, and twice as many encounters pass 10⁻⁵. A Gaussian alert rule "
     "therefore misses most of what the learned model flags, even though both models put the encounters in "
     "almost the same order (AUC 0.995). The Gaussian gets the ranking right but the level wrong, which is "
     "exactly the calibration problem we had already seen in the residuals.")
body("Third, a lesson about testing. With about 145,000 residuals per direction, any test of the Gaussian "
     "hypothesis will reject it, so a p-value tells us little about which model is better. The simulation "
     "showed a second problem: our Wilcoxon test on log-score differences rejects a true Gaussian far more "
     "often than 5%. The Student-t can always fit at least as well as the Gaussian, so the differences are "
     "not centred at zero even when H₀ is true. We therefore do not rely on its p-values. Our conclusions "
     "rest on the likelihood comparison, the goodness-of-fit distances and, above all, the out-of-sample "
     "calibration.")
body("The fleet results also show how rare high-risk encounters are: one conjunction out of 8,172 is above "
     "10⁻⁴, and the expected number of collisions in the week is 0.002. With public data a 10⁻⁴ threshold "
     "raises almost no alerts, which is why we compared the models at 10⁻⁵. Because the two models rank "
     "encounters so similarly, ranking by Pc is robust to the error model, but a fixed threshold is not. Any "
     "operational threshold should follow guidance such as [2] and be applied to Pc from a calibrated error "
     "model.")
body("The study has limitations. The newer TLE we use as the reference has its own error, so our residuals "
     "measure the combined error of two TLEs and slightly overstate that of one. Residuals that share a TLE "
     "are not independent; this makes formal p-values too optimistic and means our random train–test split "
     "has related residuals on both sides. Dropping residuals above 100 km cuts off the extreme tail. The Pc "
     "model assumes straight-line relative motion, a constant covariance during the encounter, independent "
     "errors for the two objects and a zero mean error, although the radial residuals have an offset. "
     "Finally, the power law was fitted for ages of 0.1–7 days, and anything outside that range is an "
     "extrapolation.")

# ---------------- VII. Conclusion -----------------
h1("Conclusion")
body("We estimated collision probabilities for satellite–debris encounters from public TLE data alone, by "
     "learning the position uncertainty the TLEs leave out from their own history. Across 145,467 residuals "
     "for 300 objects in the 700–900 km shell, the error grows almost linearly with TLE age, the along-track "
     "error is about three times the radial and cross-track errors, and the error is strongly heavy-tailed. A "
     "Student-t distribution with 7.1, 1.0 and 4.5 degrees of freedom fits far better than a Gaussian and "
     "gives correctly calibrated 95% intervals on held-out data. A simulation showed that a common rank test "
     "between the two models is too eager to reject, so we based the model choice on likelihood and "
     "calibration. For 8,172 conjunctions in one week, the learned Gaussian scales give a maximum Pc of "
     "4.3×10⁻⁴ and a 0.22% chance of at least one collision. Importance sampling, checked against the exact "
     "integral, carries Pc over to the heavy-tailed model. The total risk barely changes, but the risky "
     "encounters are reshuffled: twice as many exceed 10⁻⁵, and a Gaussian alert rule misses 23 of the 38 "
     "that the learned model flags.")
body("Next steps would be a better-matched importance-sampling proposal to tighten the Student-t "
     "intervals, an object-level or time-based train–test split and a simulation-calibrated test for "
     "stronger validation, including the radial offset in the Pc model, and longer TLE histories and other "
     "orbital shells.")

# ---------------- Acknowledgment -----------------
d.add_paragraph("Acknowledgment", style="Reference Head")
p = d.add_paragraph(style="Text")
r = p.add_run("We are grateful to Prof. Pallavi Mishra for her guidance throughout this project, and to "
              "Space-Track for free access to the TLE catalogue.")

# ---------------- References -----------------
d.add_paragraph("References", style="Reference Head")
refs = [
    "D. J. Kessler and B. G. Cour-Palais, “Collision frequency of artificial satellites: The creation of a debris belt,” J. Geophys. Res., vol. 83, no. A6, pp. 2637–2646, 1978.",
    "NASA, “NASA Spacecraft Conjunction Assessment and Collision Avoidance Best Practices Handbook,” NASA/SP-20205011318, Dec. 2020.",
    "F. R. Hoots and R. L. Roehrich, “Models for propagation of NORAD element sets,” Spacetrack Report No. 3, U.S. Air Force Aerospace Defense Command, 1980.",
    "D. A. Vallado, P. Crawford, R. Hujsak, and T. S. Kelso, “Revisiting Spacetrack Report #3,” in Proc. AIAA/AAS Astrodynamics Specialist Conf., Keystone, CO, 2006, AIAA 2006-6753.",
    "Space-Track.org, “Two-line element sets and satellite catalogue,” 18th Space Defense Squadron. [Online]. Available: https://www.space-track.org (accessed Oct. 2026).",
    "T. Flohrer, H. Krag, and H. Klinkrad, “Assessment and categorization of TLE orbit errors for the US SSN catalogue,” in Proc. Advanced Maui Optical and Space Surveillance Technologies Conf. (AMOS), 2008, pp. 513–524.",
    "V. P. Osweiler, “Covariance estimation and autocorrelation of NORAD two-line element sets,” M.S. thesis, Air Force Inst. of Technology, Wright-Patterson AFB, OH, 2006.",
    "J. Geul, E. Mooij, and R. Noomen, “TLE uncertainty estimation using robust weighted differencing,” Adv. Space Res., vol. 59, no. 10, pp. 2522–2535, 2017.",
    "J. L. Foster and H. S. Estes, “A parametric analysis of orbital debris collision probability and maneuver rate for space vehicles,” NASA JSC-25898, 1992.",
    "M. R. Akella and K. T. Alfriend, “Probability of collision between space objects,” J. Guid. Control Dyn., vol. 23, no. 5, pp. 769–772, 2000.",
    "S. Alfano, “A numerical implementation of spherical object collision probability,” J. Astronaut. Sci., vol. 53, no. 1, pp. 103–109, 2005.",
    "F. K. Chan, Spacecraft Collision Probability. El Segundo, CA: The Aerospace Press, 2008.",
    "K. L. Lange, R. J. A. Little, and J. M. G. Taylor, “Robust statistical modeling using the t distribution,” J. Amer. Statist. Assoc., vol. 84, no. 408, pp. 881–896, 1989.",
    "H. Akaike, “A new look at the statistical model identification,” IEEE Trans. Autom. Control, vol. 19, no. 6, pp. 716–723, 1974.",
    "F. Wilcoxon, “Individual comparisons by ranking methods,” Biometrics Bull., vol. 1, no. 6, pp. 80–83, 1945.",
    "R. Y. Rubinstein and D. P. Kroese, Simulation and the Monte Carlo Method, 3rd ed. Hoboken, NJ: Wiley, 2016.",
    "L. D. Brown, T. T. Cai, and A. DasGupta, “Interval estimation for a binomial proportion,” Statist. Sci., vol. 16, no. 2, pp. 101–133, 2001.",
    "J. L. Bentley, “Multidimensional binary search trees used for associative searching,” Commun. ACM, vol. 18, no. 9, pp. 509–517, 1975.",
]
for ref in refs:
    d.add_paragraph(ref, style="References")

d.save(OUT)
print("saved", OUT)
