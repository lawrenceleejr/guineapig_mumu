/* Common Block Declarations */

struct {
    double alp, bet, a, b, ba, bb, c, d, e, ep, sp;
} grvpar_;

#define grvpar_1 grvpar_

/* Table of constant values */

double fl(double x,double s)
{
  //    double sqrt(), log(), exp(double);
/* ****************************************************************** */
/* * Functional form of the distribution function of light partons. * */
/* * S = log( log(Q**2/.232^2) / log(.25/.232^2) ) (LO),            * */
/* * S = log( log(Q**2/.248^2) / log(.3/.248^2) )  (HO).            * */
/* ****************************************************************** */

    return (pow(x,grvpar_1.a) * (grvpar_1.ba + grvpar_1.bb * sqrt(x) 
	    + grvpar_1.c * pow(x,grvpar_1.b)) + pow(s,grvpar_1.alp) 
	    * exp(-grvpar_1.e + sqrt(-grvpar_1.ep * pow(s,grvpar_1.bet) *
	     log(x)))) * pow(1.0-x,grvpar_1.d) / x;
} /* fl */

double fh(double x,double s)
{
    /* System generated locals */

    /* Builtin functions */
  //    double sqrt(), log(), exp();

/* ************************************************* */
/* * Same as FL, but for the heavy s, c, b quarks. * */
/* ************************************************* */
/*      write (*,*) s-sp,ep**s**bet */
    if (s - grvpar_1.sp < 0.0) {
	return 0.0;
    }
    return ((s - grvpar_1.sp) * pow(x,grvpar_1.a) * (grvpar_1.ba + 
	    grvpar_1.bb * sqrt(x) + grvpar_1.c * pow(x,grvpar_1.b)) + 
	    pow(s-grvpar_1.sp,grvpar_1.alp) * exp(-grvpar_1.e + sqrt(
	    -grvpar_1.ep * pow(s,grvpar_1.bet) * log(x)))) *
	    pow(1.0-x,grvpar_1.d) / x;
} /* fh */

double grvul(double x,double q2)
{
    /* Builtin functions */
  //    double log(), sqrt();

    /* Local variables */
    static double s, s2;
    static double rts;

/* ********************************************************************** 
*/
/* * Leading order up-quark distributions. X is Bjorken-x, and Q2 is    * 
*/
/* * the factorization scale. Recall that alpha_em has been factored    * 
*/
/* * out. The parametrization is supposed to be valid for x > 10**-5,  * 
*/
/* * Q**2 < 10**6 GeV**2.                                               * 
*/
/* ********************************************************************** 
*/
    s = log(log(q2 * 18.58) / 1.536);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.alp = 1.717;
    grvpar_1.bet = .641;
    grvpar_1.a = .5 - s * .176;
    grvpar_1.b = 15. - rts * 5.687 - s2 * .552;
    grvpar_1.ba = rts * .046 + .235;
    grvpar_1.bb = .082 - s * .051 + s2 * .168;
    grvpar_1.c = s * .459;
    grvpar_1.d = .354 - s * .061;
    grvpar_1.e = s * 1.678 + 4.899;
    grvpar_1.ep = s * 1.389 + 2.046;
    return fl(x,s);
} /* grvul */

double grvdl(double x,double q2)
{
    /* Builtin functions */
  //    double log(), sqrt();

    /* Local variables */
    static double s, s2;
    static double rts;

/* ************************************ */
/* * Same as GRVUL, but for d-quarks. * */
/* ************************************ */
    s = log(log(q2 * 18.58) / 1.536);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.alp = 1.549;
    grvpar_1.bet = .782;
    grvpar_1.a = s * .026 + .496;
    grvpar_1.b = .685 - rts * .58 + s2 * .608;
    grvpar_1.ba = s * .302 + .233;
    grvpar_1.bb = s * -.818 + s2 * .198;
    grvpar_1.c = s * .154 + .114;
    grvpar_1.d = .405 - s * .195 + s2 * .046;
    grvpar_1.e = s * 1.226 + 4.807;
    grvpar_1.ep = s * .664 + 2.166;
    return fl(x,s);
} /* grvdl */

double grvgl(double x,double q2)
{
    /* Builtin functions */
  //    double log(), sqrt();

    /* Local variables */
    static double s, s2;
    static double rts;

/* ********************************** */
/* * Same as GRVUL, but for gluons. * */
/* ********************************** */
    s = log(log(q2 * 18.58) / 1.536);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.alp = .676;
    grvpar_1.bet = 1.089;
    grvpar_1.a = .462 - rts * .524;
    grvpar_1.b = 5.451 - s2 * .804;
    grvpar_1.ba = .535 - rts * .504 + s2 * .288;
    grvpar_1.bb = .364 - s * .52;
    grvpar_1.c = s2 * .115 - .323;
    grvpar_1.d = s * .79 + .233 - s2 * .139;
    grvpar_1.e = s * 1.968 + .893;
    grvpar_1.ep = s * .392 + 3.432;
    return fl(x,s);
} /* grvgl */

double grvsl(double x,double q2)
{
    /* Builtin functions */
  //    double log(), sqrt();

    /* Local variables */
    static double s, s2;
    static double rts;

/* ************************************ */
/* * Same as GRVUL, but for s-quarks. * */
/* ************************************ */
    s = log(log(q2 * 18.58) / 1.536);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.sp = 0.;
    grvpar_1.alp = 1.609;
    grvpar_1.bet = .962;
    grvpar_1.a = .47 - s2 * .099;
    grvpar_1.b = 3.246;
    grvpar_1.ba = .121 - rts * .068;
    grvpar_1.bb = s * .074 - .09;
    grvpar_1.c = s * .034 + .062;
    grvpar_1.d = s * .226 - s2 * .06;
    grvpar_1.e = s * 1.707 + 4.288;
    grvpar_1.ep = s * .656 + 2.122;
    return fh(x,s);
} /* grvsl */

double grvcl(double x,double q2)
{
    /* Builtin functions */
  //    double log();

    /* Local variables */
    static double s, s2;

/* ************************************ */
/* * Same as GRVUL, but for c-quarks. * */
/* ************************************ */
    s = log(log(q2 * 18.58) / 1.536);
    s2 = s * s;
    grvpar_1.sp = .888;
    grvpar_1.alp = .97;
    grvpar_1.bet = .545;
    grvpar_1.a = 1.254 - s * .251;
    grvpar_1.b = 3.932 - s2 * .327;
    grvpar_1.ba = s * .202 + .658;
    grvpar_1.bb = -.699;
    grvpar_1.c = .965;
    grvpar_1.d = s * .141 - s2 * .027;
    grvpar_1.e = s * .969 + 4.911;
    grvpar_1.ep = s * .952 + 2.796;
    return fh(x,s);
} /* grvcl_ */

double grvbl(double x,double q2)
{
    /* Builtin functions */
  //    double log();

    /* Local variables */
    static double s, s2;

/* ************************************ */
/* * Same as GRVUL, but for b-quarks. * */
/* ************************************ */
    s = log(log(q2 * 18.58) / 1.536);
    s2 = s * s;
    grvpar_1.sp = 1.351;
    grvpar_1.alp = 1.016;
    grvpar_1.bet = .338;
    grvpar_1.a = 1.961 - s * .37;
    grvpar_1.b = s * .119 + .923;
    grvpar_1.ba = s * .207 + .815;
    grvpar_1.bb = -2.275;
    grvpar_1.c = 1.48;
    grvpar_1.d = s * .173 - .223;
    grvpar_1.e = s * .623 + 5.426;
    grvpar_1.ep = s * .901 + 3.819;
    return fh(x,s);
} /* grvbl */

double grvuh(double x,double q2)
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    //    double log(), sqrt(), atan();

    /* Local variables */
    static double s, s2;
    static double bgamma, eq, pi, rts;

/* ********************************************************************** 
*/
/* * Higher order up-quark distributions. X is Bjorken-x, and Q2 is the * 
*/
/* * factorization scale. Recall that alpha_em has been factored out.   * 
*/
/* * The parametrization is supposed to be  valid for x > 10**-5,       * 
*/
/* * Q**2 < 10**6 GeV**2.                                               * 
*/
/* * WARNING: These distributions are in the DIS_gamma scheme; to trans-* 
*/
/* * late them into MSbar, a term proportional to B_gamma must be added.* 
*/
/* *--------------------------------------------------------------------* 
*/
/* * Modification by J. Zunft on 11/06/92: MS_bar scheme                * 
*/
/* ********************************************************************** 
*/
    s = log(log(q2 * 16.26) / 1.585);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.alp = .583;
    grvpar_1.bet = .688;
    grvpar_1.a = .449 - s * .025 - s2 * .071;
    grvpar_1.b = 5.06 - rts * 1.116;
    grvpar_1.ba = .103;
    grvpar_1.bb = s * .422 + .319;
    grvpar_1.c = s * 4.792 + 1.508 - s2 * 1.963;
    grvpar_1.d = rts * .222 + 1.075 - s2 * .193;
    grvpar_1.e = s * 1.131 + 4.147;
    grvpar_1.ep = s * .874 + 1.661;
    ret_val = fl(x,s);
    pi = atan(1.) * 4.;
    eq = .66666666666666663;
/* Computing 2nd power */
    bgamma = (((float)1. - x * (float)2. + x * x * (float)2.) * log(((
	    float)1. - x) / x) - (float)1. + x * (float)8. * ((float)1. - 
	    x)) * (float)4.;
/* Computing 2nd power */
    ret_val -= 3.0*eq*eq/(8.0*pi) * bgamma;
    return ret_val;
} /* grvuh_ */

double grvdh(double x,double q2)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Builtin functions */
    //    double log(), sqrt(), atan();

    /* Local variables */
    static double s, s2;
    static double bgamma, eq, pi, rts;

/* ************************************ */
/* * Same as GRVUH, but for d-quarks. * */
/* ************************************ */
    s = log(log(q2 * 16.26) / 1.585);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.alp = .591;
    grvpar_1.bet = .698;
    grvpar_1.a = .442 - s * .132 - s2 * .058;
    grvpar_1.b = 5.437 - rts * 1.916;
    grvpar_1.ba = .099;
    grvpar_1.bb = .311 - s * .059;
    grvpar_1.c = s * .078 + .8 - s2 * .1;
    grvpar_1.d = rts * .294 + .862 - s2 * .184;
    grvpar_1.e = s * 1.352 + 4.202;
    grvpar_1.ep = s * .99 + 1.841;
    ret_val = fl(x,s);
    pi = atan(1.) * 4.;
    eq = .33333333333333331;
/* Computing 2nd power */
    d__1 = x;
    bgamma = (((float)1. - x * (float)2. + d__1 * d__1 * (float)2.) * log(((
	    float)1. - x) / x) - (float)1. + x * (float)8. * ((float)1. - 
	    x)) * (float)4.;
/* Computing 2nd power */
    d__1 = eq;
    ret_val -= d__1 * d__1 * 3. / (pi * 8.) * bgamma;
    return ret_val;
} /* grvdh */

double grvgh(double x,double q2)
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    //    double log(), sqrt();

    /* Local variables */
    static double s, s2;
    static double rts;

/* **************************************************************** */
/* * Same as GRVUH, but for gluons. Note that DIS_gamma and MSbar * */
/* * scheme are identical for gluons (to the given order).        * */
/* **************************************************************** */
    s = log(log(q2 * 16.26) / 1.585);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.alp = 1.161;
    grvpar_1.bet = 1.591;
    grvpar_1.a = .53 - rts * .742 + s2 * .025;
    grvpar_1.b = 5.662;
    grvpar_1.ba = .533 - rts * .281 + s2 * .218;
    grvpar_1.bb = .025 - s * .518 + s2 * .156;
    grvpar_1.c = s2 * .209 - .282;
    grvpar_1.d = s * 1.058 + .107 - s2 * .218;
    grvpar_1.e = s * 2.704;
    grvpar_1.ep = 3.071 - s * .378;
    return fl(x,s);
} /* grvgh */

double grvsh(double x,double q2)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Builtin functions */
    //    double log(), sqrt(), atan();

    /* Local variables */
    static double s, s2;
    static double bgamma, eq, pi, rts;

/* ************************************ */
/* * Same as GRVUH, but for s-quarks. * */
/* ************************************ */
    s = log(log(q2 * 16.26) / 1.585);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.sp = 0.;
    grvpar_1.alp = .635;
    grvpar_1.bet = .456;
    grvpar_1.a = 1.77 - rts * .735 - s2 * .079;
    grvpar_1.b = 3.832;
    grvpar_1.ba = .084 - s * .023;
    grvpar_1.bb = .136;
    grvpar_1.c = 2.119 - s * .942 + s2 * .063;
    grvpar_1.d = s * .076 + 1.271 - s2 * .19;
    grvpar_1.e = s * .737 + 4.604;
    grvpar_1.ep = s * .976 + 1.641;
    ret_val = fh(x,s);
    pi = atan(1.) * 4.;
    eq = .33333333333333331;
/* Computing 2nd power */
    d__1 = x;
    bgamma = (((float)1. - x * (float)2. + d__1 * d__1 * (float)2.) * log(((
	    float)1. - x) / x) - (float)1. + x * (float)8. * ((float)1. - 
	    x)) * (float)4.;
/* Computing 2nd power */
    d__1 = eq;
    ret_val -= d__1 * d__1 * 3. / (pi * 8.) * bgamma;
    return ret_val;
} /* grvsh_ */

double grvch(double x,double q2)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Builtin functions */
    //    double log(), atan();

    /* Local variables */
    static double s, s2;
    static double bgamma, eq, pi;

/* ************************************ */
/* * Same as GRVUH, but for c-quarks. * */
/* ************************************ */
    s = log(log(q2 * 16.26) / 1.585);
    s2 = s * s;
    grvpar_1.sp = .82;
    grvpar_1.alp = .926;
    grvpar_1.bet = .152;
    grvpar_1.a = 1.142 - s * .175;
    grvpar_1.b = 3.276;
    grvpar_1.ba = s * .317 + .504;
    grvpar_1.bb = -.433;
    grvpar_1.c = 3.334;
    grvpar_1.d = s * .326 + .398 - s2 * .107;
    grvpar_1.e = s * .408 + 5.493;
    grvpar_1.ep = s * 1.277 + 2.426;
    ret_val = fh(x,s);
    pi = atan(1.) * 4.;
    eq = .66666666666666663;
/* Computing 2nd power */
    d__1 = x;
    bgamma = (((float)1. - x * (float)2. + d__1 * d__1 * (float)2.) * log(((
	    float)1. - x) / x) - (float)1. + x * (float)8. * ((float)1. - 
	    x)) * (float)4.;
/* Computing 2nd power */
    d__1 = eq;
    ret_val -= d__1 * d__1 * 3. / (pi * 8.) * bgamma;
    return ret_val;
} /* grvch_ */

double grvbh(double x,double q2)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Builtin functions */
    //    double log(), atan();

    /* Local variables */
    static double s, s2;
    static double bgamma, eq, pi;

/* ************************************ */
/* * Same as GRVUH, but for b-quarks. * */
/* ************************************ */
    s = log(log(q2 * 16.26) / 1.585);
    s2 = s * s;
    grvpar_1.sp = 1.297;
    grvpar_1.alp = .969;
    grvpar_1.bet = .266;
    grvpar_1.a = 1.953 - s * .391;
    grvpar_1.b = 1.657 - s * .161;
    grvpar_1.ba = s * .034 + 1.076;
    grvpar_1.bb = -2.015;
    grvpar_1.c = 1.662;
    grvpar_1.d = s * .016 + .353;
    grvpar_1.e = s * .249 + 5.713;
    grvpar_1.ep = s * .673 + 3.456;
    ret_val = fh(x,s);
    pi = atan(1.) * 4.;
    eq = .33333333333333331;
/* Computing 2nd power */
    d__1 = x;
    bgamma = (((float)1. - x * (float)2. + d__1 * d__1 * (float)2.) * log(((
	    float)1. - x) / x) - (float)1. + x * (float)8. * ((float)1. - 
	    x)) * (float)4.;
/* Computing 2nd power */
    d__1 = eq;
    ret_val -= d__1 * d__1 * 3. / (pi * 8.) * bgamma;
    return ret_val;
} /* grvbh */

double grvuh0(double x,double q2)
{
    /* Builtin functions */
  //    double log(), sqrt();

    /* Local variables */
    static double s, s2;
    static double rts;

/* ********************************************************************** 
*/
/* * The leading order part of the total HO up-quark distribution. This * 
*/
/* * part is needed since, strictly speaking, only this part should be  * 
*/
/* * multiplied with the HO part of a NLO cross section or Wilson co-   * 
*/
/* * efficient. For the meaning of the variables, see GRVUH.            * 
*/
/* ********************************************************************** 
*/
    s = log(log(q2 * 16.26) / 1.585);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.alp = 1.447;
    grvpar_1.bet = .848;
    grvpar_1.a = s * .2 + .527 - s2 * .107;
    grvpar_1.b = 7.106 - rts * .31 - s2 * .786;
    grvpar_1.ba = s * .533 + .197;
    grvpar_1.bb = .062 - s * .398 + s2 * .109;
    grvpar_1.c = s * .755 - s2 * .112;
    grvpar_1.d = .318 - s * .059;
    grvpar_1.e = s * 1.708 + 4.225;
    grvpar_1.ep = s * .866 + 1.752;
    return fl(x,s);
} /* grvuh0_ */

double grvdh0(double x,double q2)
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    //    double log(), sqrt();

    /* Local variables */
    static double s, s2;
    static double rts;

/* ************************************* */
/* * Same as GRVUH0, but for d-quarks. * */
/* ************************************* */
    s = log(log(q2 * 16.26) / 1.585);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.alp = 1.424;
    grvpar_1.bet = .77;
    grvpar_1.a = rts * .067 + .5 - s2 * .055;
    grvpar_1.b = .376 - rts * .453 + s2 * .405;
    grvpar_1.ba = s * .184 + .156;
    grvpar_1.bb = s * -.528 + s2 * .146;
    grvpar_1.c = s * .092 + .121;
    grvpar_1.d = .379 - s * .301 + s2 * .081;
    grvpar_1.e = s * 1.638 + 4.346;
    grvpar_1.ep = s * 1.016 + 1.645;
    ret_val = fl(x,s);
    return ret_val;
} /* grvdh0 */

double grvgh0(double x,double q2)
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    //    double log(), sqrt();

    /* Local variables */
    static double s, s2;
    static double rts;

/* ***************************************************************** */
/* * Same as GRVUH0, but for gluons. Note that DIS_gamma and MSbar * */
/* * scheme are identical for gluons (to the given order).         * */
/* ***************************************************************** */
    s = log(log(q2 * 16.26) / 1.585);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.alp = .661;
    grvpar_1.bet = .793;
    grvpar_1.a = .537 - rts * .6;
    grvpar_1.b = 6.389 - s2 * .953;
    grvpar_1.ba = .558 - rts * .383 + s2 * .261;
    grvpar_1.bb = s * -.305;
    grvpar_1.c = s2 * .078 - .222;
    grvpar_1.d = s * .978 + .153 - s2 * .209;
    grvpar_1.e = s * 1.772 + 1.429;
    grvpar_1.ep = s * .806 + 3.331;
    return fl(x,s);
} /* grvgh0 */

double grvsh0(double x,double q2)
{
    /* Builtin functions */
  //    double log(), sqrt();

    /* Local variables */
    static double s, s2;
    static double rts;

/* ************************************* */
/* * Same as GRVUH0, but for s-quarks. * */
/* ************************************* */
    s = log(log(q2 * 16.26) / 1.585);
    rts = sqrt(s);
    s2 = s * s;
    grvpar_1.sp = 0.;
    grvpar_1.alp = 1.578;
    grvpar_1.bet = .863;
    grvpar_1.a = s * .332 + .622 - s2 * .3;
    grvpar_1.b = 2.469;
    grvpar_1.ba = .211 - rts * .064 - s2 * .018;
    grvpar_1.bb = s * .122 - .215;
    grvpar_1.c = .153;
    grvpar_1.d = s * .253 - s2 * .081;
    grvpar_1.e = s * 2.014 + 3.99;
    grvpar_1.ep = s * .986 + 1.72;
    return fh(x,s);
} /* grvsh0 */

double grvch0(double x,double q2)
{
    /* Builtin functions */
  //    double log();

    /* Local variables */
    static double s, s2;

/* ************************************* */
/* * Same as GRVUH0, but for c-quarks. * */
/* ************************************* */
    s = log(log(q2 * 16.26) / 1.585);
    s2 = s * s;
    grvpar_1.sp = .82;
    grvpar_1.alp = .929;
    grvpar_1.bet = .381;
    grvpar_1.a = 1.228 - s * .231;
    grvpar_1.b = 3.806 - s2 * .337;
    grvpar_1.ba = s * .15 + .932;
    grvpar_1.bb = -.906;
    grvpar_1.c = 1.133;
    grvpar_1.d = s * .138 - s2 * .028;
    grvpar_1.e = s * .628 + 5.588;
    grvpar_1.ep = s * 1.054 + 2.665;
    return fh(x,s);
} /* grvch0 */

double grvbh0(double x,double q2)
{
    /* Builtin functions */
  //    double log();

    /* Local variables */
    static double s, s2;

/* ************************************* */
/* * Same as GRVUH0, but for b-quarks. * */
/* ************************************* */
    s = log(log(q2 * 16.26) / 1.585);
    s2 = s * s;
    grvpar_1.sp = 1.297;
    grvpar_1.alp = .97;
    grvpar_1.bet = .207;
    grvpar_1.a = 1.719 - s * .292;
    grvpar_1.b = s * .096 + .928;
    grvpar_1.ba = s * .178 + .845;
    grvpar_1.bb = -2.31;
    grvpar_1.c = 1.558;
    grvpar_1.d = s * .151 - .191;
    grvpar_1.e = s * .282 + 6.089;
    grvpar_1.ep = s * 1.062 + 3.379;
    return fh(x,s);
} /* grvbh0 */
