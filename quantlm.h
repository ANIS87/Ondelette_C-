#ifndef QUANTLM_H
#define QUANTLM_H

/* Quantizes the signal x of size n with nq quantized values
 * Returns the quantized signal of x
 */

void quantlm(double* x,int n,int nq);

/* Quantizes the signal x of size n with nq quantized values
 * Returns the indexed quantifier of each value of x
 */
void quantlm_idx(double* x,int n,int nq);

#endif
