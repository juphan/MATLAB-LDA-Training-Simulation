/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mav_calc.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 04-Aug-2019 18:45:42
 */

/* Include Files */
#include <math.h>
#include "mav_calc.h"

/* Function Definitions */

/*
 * Arguments    : const float emg[100]
 * Return Type  : float
 */
float mav_calc(const float emg[100])
{
  float y;
  int k;
  float b_y[100];
  for (k = 0; k < 100; k++) {
    b_y[k] = fabsf(emg[k]);
  }

  y = b_y[0];
  for (k = 0; k < 99; k++) {
    y += b_y[k + 1];
  }

  y /= 100.0F;
  return y;
}

/*
 * File trailer for mav_calc.c
 *
 * [EOF]
 */
