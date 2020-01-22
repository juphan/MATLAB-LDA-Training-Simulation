/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: wav_calc.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 04-Aug-2019 18:57:54
 */

/* Include Files */
#include <math.h>
#include "wav_calc.h"

/* Function Definitions */

/*
 * Arguments    : const float emg[100]
 * Return Type  : float
 */
float wav_calc(const float emg[100])
{
  float y;
  int k;
  float b_y[99];
  for (k = 0; k < 99; k++) {
    b_y[k] = fabsf(emg[k] - emg[k + 1]);
  }

  y = b_y[0];
  for (k = 0; k < 98; k++) {
    y += b_y[k + 1];
  }

  y /= 99.0F;
  return y;
}

/*
 * File trailer for wav_calc.c
 *
 * [EOF]
 */
