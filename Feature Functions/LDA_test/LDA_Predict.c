/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: LDA_Predict.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 05-Aug-2019 18:19:58
 */

/* Include Files */
#include "LDA_Predict.h"

/* Function Definitions */

/*
 * Arguments    : const float features[4]
 *                const float Wg[20]
 *                const float Cg[5]
 * Return Type  : float
 */
float LDA_Predict(const float features[4], const float Wg[20], const float Cg[5])
{
  float prediction;
  int i0;
  float b_max;
  float tmp[5];
  int i1;
  for (i0 = 0; i0 < 5; i0++) {
    i1 = i0 << 2;
    tmp[i0] = (((features[0] * Wg[i1] + features[1] * Wg[1 + i1]) + features[2] *
                Wg[2 + i1]) + features[3] * Wg[3 + i1]) + Cg[i0];
  }

  b_max = tmp[0];
  prediction = 1.0F;
  if (tmp[1] > tmp[0]) {
    b_max = tmp[1];
    prediction = 2.0F;
  }

  if (tmp[2] > b_max) {
    b_max = tmp[2];
    prediction = 3.0F;
  }

  if (tmp[3] > b_max) {
    b_max = tmp[3];
    prediction = 4.0F;
  }

  if (tmp[4] > b_max) {
    prediction = 5.0F;
  }

  return prediction;
}

/*
 * File trailer for LDA_Predict.c
 *
 * [EOF]
 */
