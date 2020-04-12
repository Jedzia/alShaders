#pragma once

#include <ai.h>
#include "alUtil.h"
#include <stdint.h>

inline bool refraction(const AtVector& I, const AtVector& N, float eta, AtVector& T) 
{
   // compute refracted direction
   // return value will be false if TIR occurs
   // NOTE: I is the incoming ray direction (points toward the surface, normalized)
   //       N is the surface normal (points toward the incoming ray origin, normalized)
   //       T is the outgoing refracted direction (points away from the surface)
   //AtVector T;
   float cosi = -AiV3Dot(I, N);
   // check which side of the surface we are on
   AtVector Nn; float neta;
   // if (cosi > 0) 
   // {
      // we are on the outside of the surface, going in
      neta = 1 / eta;
      Nn = N;
   // } 
   // else 
   // {
   //    // we are inside the surface,
   //    cosi = -cosi;
   //    neta = eta;
   //    Nn = -N;
   // }
   float arg = 1.0f - (neta * neta * (1.0f - cosi * cosi));
   if (arg >= 0) 
   {
      float dnp = sqrtf(arg);
      float nK = (neta * cosi) - dnp;
      T = I * neta + Nn * nK;
      return true;
   }
   else
   {
      T = I - 2 * cosi * -N;
      return false;
   }
   
}

inline float fast_log2(float x)
{
  union { float f; uint32_t i; } vx = { x };
  union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
  float y = vx.i;
  y *= 1.1920928955078125e-7f;

  return y - 124.22551499f
           - 1.498030302f * mx.f 
           - 1.72587999f / (0.3520887068f + mx.f);
}

#ifdef _WIN32
// thanks Windows!
#if (_MSC_VER < 1800) // VS2013 has erff
float erff(float x)
{
    // constants
    float a1 =  0.254829592;
    float a2 = -0.284496736;
    float a3 =  1.421413741;
    float a4 = -1.453152027;
    float a5 =  1.061405429;
    float p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    float t = 1.0/(1.0 + p*x);
    float y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}
#endif
#endif

inline float fast_ierf(float x)
{
   static const float invk = 0.30004578719350504f;
   static const float a = 0.020287853348211326f;
   static const float b = 0.07236892874789555f;
   static const float c = 0.9913030456864257f;
   static const float d = 0.8059775923760193f;

   float xsq = x * x;

   return invk * fast_log2 ((1.0f + x) / (1.0f - x)) 
       + x * (a - b * xsq) / (c - d * xsq);
}


struct MicrofacetTransmission
{
   inline AtVector2 sampleSlope(float cos_theta, float u1, float u2) const
   {
      static const float INV_SQRT_PI = 0.5641895835477563f;
      cos_theta = std::max(cos_theta, 1e-6f);
      const float cos_theta2 = SQR(cos_theta);
      const float tan_theta = sqrtf(1.0f - cos_theta2) / cos_theta;
      const float cot_theta = 1.0f / tan_theta;

      //float c = fast_erf(cot_theta);
      float c = erff(cot_theta);
      float k = tan_theta * INV_SQRT_PI;
      float y_approx = u1 * (1.0f + c + k * (1.0f - SQR(c)));
      float y_exact  = u1 * (1.0f + c + k * AiFastExp(-SQR(cot_theta)));
      float b = k > 0.0f ? (0.5f - sqrtf(k * (k - y_approx + 1.0f) + 0.25f)) / k : y_approx - 1.0f;

      float inv_erf = fast_ierf(b);
      float value = 1.0f + b + k * AiFastExp(-SQR(inv_erf)) - y_exact;

      AtVector2 slope;
      if (fabsf(value) > 1e-6f)
      {
         b -= value / (1.0f - inv_erf * tan_theta);
         inv_erf = fast_ierf(b);
         value = 1.0f + b + k * AiFastExp(-SQR(inv_erf)) - y_exact;
         b -= value / (1.0f - inv_erf * tan_theta);

         slope.x = fast_ierf(b);
      }
      else
      {
         slope.x = inv_erf;
      }

      slope.y = fast_ierf(2.0f * u2 - 1.0f);

      return slope;
   }

   inline AtVector sampleMicrofacetNormal(float u1, float u2) const
   {
      AtVector s_omega_o(AiV3Dot(U, omega_o), AiV3Dot(V, omega_o), AiV3Dot(N, omega_o));
      s_omega_o.x *= alpha_x;
      s_omega_o.y *= alpha_y;
      s_omega_o = AiV3Normalize(s_omega_o);

      float cos_theta = std::max(s_omega_o.z, 0.0f);
      float cos_phi = 1.0f;
      float sin_phi = 0.0f;

      if (cos_theta < -0.99999f)
      {
         float inv_norm = 1.0f / sqrtf(SQR(s_omega_o.x) + SQR(omega_o.y));
         cos_phi = s_omega_o.x * inv_norm;
         sin_phi = s_omega_o.y * inv_norm;
      }

      AtVector2 slope = sampleSlope(cos_theta, u1, u2);

      AtVector2 s = AtVector2(cos_phi * slope.x - sin_phi * slope.y, sin_phi * slope.x + cos_phi * slope.y);
      s.x *= alpha_x;
      s.y *= alpha_y;

      AtVector m;
      const float mlen = sqrtf(SQR(s.x) + SQR(s.y) + 1.0f);
      m.x = fabsf(s.x) < mlen ? -s.x / mlen : 1.0f;
      m.y = fabsf(s.y) < mlen ? -s.y / mlen : 1.0f;
      m.z = 1.0f / mlen;

      return U * m.x + V * m.y + N * m.z;
   }

   float D(const AtVector& M) const
   {
      float d = 0.0f;

      const float cos_theta = AiV3Dot(M, N);
      if (cos_theta > 1e-7f)
      {
         const float cos_theta2 = SQR(cos_theta);
         const float cos_theta4 = SQR(cos_theta2);
         const float sx = AiV3Dot(M, U) / (cos_theta * alpha_x);
         const float sy = AiV3Dot(M, V) / (cos_theta * alpha_y);
         const float tan_theta2 = SQR(sx) + SQR(sy);

         d = AiFastExp(-tan_theta2) / (AI_PI * alpha_x * alpha_y * cos_theta4);
         assert(AiIsFinite(d));
      }

      return d;
   }

   float lambda(const AtVector& omega) const
   {
      const float cos_theta2 = SQR(AiV3Dot(N, omega));
      const float cos_phi_2 = SQR(AiV3Dot(U, omega) * alpha_x);
      const float sin_phi_2 = SQR(AiV3Dot(V, omega) * alpha_y);
      const float a2 = cos_theta2 / (cos_phi_2 + sin_phi_2);
      const float a = sqrtf(a2);
      return a < 1.6f ? (1.0f - 1.259f * a + 0.396f * a2) / (3.535f * a + 2.181f * a2) : 0.0f;
   }

   inline static float G2(float lambda_o, float lambda_i)
   {
      return 1.0f / (lambda_i + lambda_o + 1.0f);
   }

   inline static float G1(float lambda_o)
   {
      return 1.0f / (lambda_o + 1.0f);
   }

   // if TIR, returned vector is zero
   inline bool sample(float u1, float u2, AtVector& T) const
   {
      const AtVector m = sampleMicrofacetNormal(u1, u2);
      return refraction(sg->Rd, m, eta, T);
   }

   inline AtRGB btdf(const AtVector& omega_i) const
   {
      AtRGB result = AI_RGB_BLACK;

      if (cos_N_o > 0 && AiV3Dot(omega_i, omega_o) <= 0.0f && !AiV3IsSmall(omega_i))
      {
         AtVector H = -(eta * omega_i + omega_o);
         if (eta < 1.0f) H = -H;
         const float inv_h2 = 1.0f / AiV3Dot(H, H);
         H = AiV3Normalize(H);

         const float cos_H_o = fabsf(AiV3Dot(H, omega_o));

         const float cos_H_i = AiV3Dot(H, omega_i);
         const float cos_theta = AiV3Dot(H, N);
         if (cos_theta > 0.0f)
         {
            const float d = D(H);
            const float lambda_o = lambda(omega_o);
            const float lambda_i = lambda(omega_i);
            const float g2 = G2(lambda_o, lambda_i);
            const float g1 = G1(lambda_o);

            result = (fabsf(cos_H_i * cos_H_o) * (eta*eta) * (g2 * d) * inv_h2) / fabsf(cos_N_o);
            assert(AiIsFinite(result));
            // cancelling out all the like terms from the pdf in the above leaves us:
            // result = g2;

         }

      }

      return result;
   }

   inline float pdf(const AtVector& omega_i) const
   {
      float result = 0.0f;

      if (cos_N_o > 0 && AiV3Dot(omega_i, omega_o) <= 0.0f && !AiV3IsSmall(omega_i))
      {
         AtVector H = -(eta * omega_i + omega_o);
         if (eta < 1.0f) H = -H;
         const float inv_h2 = 1.0f / AiV3Dot(H, H);
         H = AiV3Normalize(H);

         const float cos_H_o = AiV3Dot(H, omega_o);

         const float cos_H_i = AiV3Dot(H, omega_i);
         const float cos_theta = AiV3Dot(H, N);
         if (cos_theta > 0.0f)
         {
            const float d = D(H);
            const float lambda_o = lambda(omega_o);
            const float lambda_i = lambda(omega_i);
            const float g2 = G2(lambda_o, lambda_i);
            const float g1 = G1(lambda_o);

            result = (fabsf(cos_H_i * cos_H_o) * (eta*eta) * (g1 * d) * inv_h2) / fabsf(cos_N_o);
            // cancelling out all the like terms from the btdf in the above leaves us:
            // result = g1;
            assert(AiIsFinite(result));

         }

      }

      return result;
   }

   static MicrofacetTransmission* create(const AtShaderGlobals* sg, float rx, float ry, float eta, const AtVector& N, const AtVector& U, const AtVector& V)
   {
      MicrofacetTransmission* mt = reinterpret_cast<MicrofacetTransmission*>(AiShaderGlobalsQuickAlloc(sg, sizeof(MicrofacetTransmission)));
      mt->sg = sg;
      mt->alpha_x = std::max(0.001f, rx);
      mt->alpha_y = std::max(0.001f, ry);
      mt->eta = eta;
      mt->N = N;
      mt->U = U;
      mt->V = V;
      mt->omega_o = -sg->Rd;
      mt->cos_N_o = AiV3Dot(N, mt->omega_o);
      return mt;
   }

   static AtVector Sample(const void* brdf_data, float u1, float u2)
   {
      const MicrofacetTransmission* mt = reinterpret_cast<const MicrofacetTransmission*>(brdf_data);
      AtVector T;
      if (mt->sample(u1, u2, T))
      {
         return T;
      }
      else
      {
         return AtVector(0.0f, 0.0f, 0.0f);
      }
   }

   static AtRGB BTDF(const void* brdf_data, const AtVector* omega_i)
   {
      const MicrofacetTransmission* mt = reinterpret_cast<const MicrofacetTransmission*>(brdf_data);
      return mt->btdf(*omega_i);
   }

   static float PDF(const void* brdf_data, const AtVector* omega_i)
   {
      const MicrofacetTransmission* mt = reinterpret_cast<const MicrofacetTransmission*>(brdf_data);
      return mt->pdf(*omega_i);
   }

   const AtShaderGlobals* sg;
   float alpha_x;
   float alpha_y;
   float eta;
   AtVector U;
   AtVector V;
   AtVector N;
   AtVector omega_o;
   float cos_N_o;
};